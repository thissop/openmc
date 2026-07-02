#!/usr/bin/env python
"""dagmc_writer.py -- write a DAGMC .h5m for the conformal nested-shell wall
WITHOUT pymoab / MOAB python bindings, using only numpy + h5py.

Why this exists: every STL->DAGMC path (stl_to_h5m, stellarmesh, cad_to_dagmc)
imports pymoab, which conda-forge does not ship in importable form. But the DAGMC
.h5m is just an HDF5 file in MOAB's documented "mhdf" layout, and this geometry is
the easy case: nested CLOSED toroidal shells that share no edges, so there are no
curve/vertex geometry sets, only surfaces + volumes + senses. We build the boundary
surfaces ONCE (so adjacent shells share a surface by construction -> a merged model,
no coincident-surface welding needed), reusing the exact vertex/triangle routines
in stellarator_geometry.py so the mesh is byte-identical to the validated build.

Topology (6 material shells W/steel/Be/FLiBe/shield/coil; 'sol' is vacuum):
  7 boundary surfaces  S0..S6  at cumulative offsets from the LCFS
  6 volumes            V0..V5  (Vi bounded by Si inner, S(i+1) outer)
  6 material groups    NAME="mat:<layer>"  (run_conformal matches Material.name)
  implicit complement  = plasma interior (inside S0) + sol gap + exterior (outside S6)

Serialization is grounded in the MOAB H5M format spec:
  https://sigma.mcs.anl.gov/moab/h5m-file-format/
Structural self-consistency is checked by reading the file back with h5py
(_verify_structure). DAGMC's *semantic* acceptance is confirmed separately by
loading the file with openmc.DAGMCUniverse (run --load-test in the conda env).

Usage:
  python dagmc_writer.py --selftest                 # synthetic 2-shell, no repo deps
  python dagmc_writer.py --stem equil_precise_qa \
         --geom-dir spf_prototype/data/equil_precise_qa_geom \
         --out equil_precise_qa.h5m                 # real geometry
  python dagmc_writer.py --load-test equil_precise_qa.h5m   # DAGMC accept? (needs openmc)
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore", category=DeprecationWarning)  # quiet h5py/numpy noise

try:
    import h5py
except ImportError:
    sys.exit("h5py required: pip install h5py")

# ---- MOAB mhdf constants (from the format spec / moab-dev) -------------------
MESHSET_SET = 0x02          # MOAB entityset flag: unordered set (geometry sets)
MHDF_SPARSE = 1             # tag 'class' attribute: sparse storage
CAT_SIZE = 32               # CATEGORY / NAME tags are opaque char[32]
GEOM_DIMENSION = "GEOM_DIMENSION"
GLOBAL_ID = "GLOBAL_ID"
CATEGORY = "CATEGORY"
GEOM_SENSE_2 = "GEOM_SENSE_2"
NAME = "NAME"


# =============================================================================
# 1. SERIALIZER  --  assembled topology -> MOAB native .h5m
# =============================================================================
def _pad(s: str, n: int = CAT_SIZE) -> bytes:
    b = s.encode("ascii")
    if len(b) >= n:
        raise ValueError(f"tag string {s!r} too long for {n}-byte field")
    return b + b"\x00" * (n - len(b))


def write_h5m(path, node_xyz, tri_conn, surfaces, volumes, groups):
    """Write the MOAB native .h5m.

    node_xyz : (Nv,3) float64 vertex coordinates (cm)
    tri_conn : (Nt,3) int, 0-based indices into node_xyz (global)
    surfaces : list of dicts, in surface order:
                 {'tris': int array of GLOBAL 0-based triangle rows,
                  'gid' : int surface global id,
                  'sense': (fwd_vol, rev_vol)}  vol = 0-based volume index or None
    volumes  : list of dicts, in volume order:
                 {'gid': int, 'surfs': [inner_surf_idx, outer_surf_idx]}
    groups   : list of dicts, in group order:
                 {'gid': int, 'name': 'mat:W', 'vols': [vol_idx,...]}
    """
    node_xyz = np.ascontiguousarray(node_xyz, dtype="f8")
    tri_conn = np.ascontiguousarray(tri_conn, dtype="i8")
    Nv, Nt = len(node_xyz), len(tri_conn)
    nsurf, nvol, ngrp = len(surfaces), len(volumes), len(groups)

    # --- unified positive ID space, contiguous per entity block --------------
    NODE_START = 1
    TRI_START = NODE_START + Nv
    SET_START = TRI_START + Nt
    nsets = nsurf + nvol + ngrp
    MAX_ID = SET_START + nsets - 1

    surf_id = [SET_START + i for i in range(nsurf)]
    vol_id = [SET_START + nsurf + i for i in range(nvol)]
    grp_id = [SET_START + nsurf + nvol + i for i in range(ngrp)]

    def vhandle(v):  # volume index -> file id, or 0 for implicit complement
        return 0 if v is None else vol_id[v]

    # --- set meta table + flat contents/children/parents ---------------------
    # sets are written in order: [surfaces..., volumes..., groups...]
    # spec: each list column stores the index of the LAST entry (cumulative-1),
    # so an empty list repeats the previous end index and the very first is -1.
    contents, children, parents, setlist = [], [], [], []

    def row():
        return (len(contents) - 1, len(children) - 1, len(parents) - 1, MESHSET_SET)

    for si, s in enumerate(surfaces):                     # SURFACES
        contents.extend(int(TRI_START + t) for t in s["tris"])   # tris = contents
        par = [vol_id[vi] for vi, v in enumerate(volumes) if si in v["surfs"]]
        parents.extend(par)                               # bounding volumes = parents
        setlist.append(row())

    for v in volumes:                                     # VOLUMES
        children.extend(surf_id[si] for si in v["surfs"])  # bounding surfaces = children
        setlist.append(row())

    for g in groups:                                      # GROUPS
        contents.extend(vol_id[vi] for vi in g["vols"])   # member volumes = contents
        setlist.append(row())

    setlist = np.asarray(setlist, dtype="i8").reshape(nsets, 4)
    contents = np.asarray(contents, dtype="i8")
    children = np.asarray(children, dtype="i8")
    parents = np.asarray(parents, dtype="i8")

    # --- sparse tag tables ---------------------------------------------------
    geom_ids = surf_id + vol_id                     # dim tag: surfaces(2)+volumes(3)
    geom_dim = [2] * nsurf + [3] * nvol
    gid_ids = surf_id + vol_id + grp_id             # global id: everything
    gid_val = ([s["gid"] for s in surfaces] + [v["gid"] for v in volumes]
               + [g["gid"] for g in groups])
    cat_ids = surf_id + vol_id + grp_id
    cat_val = ([_pad("Surface")] * nsurf + [_pad("Volume")] * nvol
               + [_pad("Group")] * ngrp)
    sense_ids = surf_id                              # 2 volume handles per surface
    sense_val = np.array([[vhandle(s["sense"][0]), vhandle(s["sense"][1])]
                          for s in surfaces], dtype="u8").reshape(nsurf, 2)
    name_ids = grp_id
    name_val = [_pad(g["name"]) for g in groups]

    opaque = h5py.opaque_dtype(np.dtype(f"V{CAT_SIZE}"))  # MOAB char[32] opaque
    ELEM_TOPOS = ["Edge", "Tri", "Quad", "Polygon", "Tet", "Pyramid",
                  "Prism", "Knife", "Hex", "Polyhedron"]        # mhdf elemtypes enum
    elem_enum = h5py.enum_dtype({n: i for i, n in enumerate(ELEM_TOPOS)}, basetype="i4")

    with h5py.File(path, "w") as f:
        t = f.create_group("tstt")
        t.attrs.create("max_id", MAX_ID, dtype="u8")
        # history: list of variable-length strings (optional; write a stamp)
        vlen = h5py.string_dtype(encoding="ascii")
        f.create_dataset("tstt/history", data=np.array(
            ["dagmc_writer.py", "conformal nested-shell, pymoab-free"], dtype=object),
            dtype=vlen)

        # nodes
        d = f.create_dataset("tstt/nodes/coordinates", data=node_xyz)
        d.attrs.create("start_id", NODE_START, dtype="i8")
        f.create_group("tstt/nodes/tags")

        # committed enum of element topologies (spec requires /tstt/elemtypes)
        f["tstt/elemtypes"] = elem_enum

        # elements: 3-node triangles; group name = topology + conn length = "Tri3".
        # element_type attribute MUST be the elemtypes ENUM selecting "Tri" -- a plain
        # int attribute makes MOAB fail with "not an enumeration data type".
        conn = (tri_conn + NODE_START).astype("i8")
        d = f.create_dataset("tstt/elements/Tri3/connectivity", data=conn)
        d.attrs.create("start_id", TRI_START, dtype="i8")
        f["tstt/elements/Tri3"].attrs.create(
            "element_type", ELEM_TOPOS.index("Tri"), dtype=elem_enum)
        f.create_group("tstt/elements/Tri3/tags")

        # sets
        d = f.create_dataset("tstt/sets/list", data=setlist)
        d.attrs.create("start_id", SET_START, dtype="i8")
        f.create_dataset("tstt/sets/contents", data=contents)
        f.create_dataset("tstt/sets/children", data=children)
        f.create_dataset("tstt/sets/parents", data=parents)
        f.create_group("tstt/sets/tags")

        # tags
        tg = f.create_group("tstt/tags")

        def sparse_tag(name, ids, values, dtype, *, is_handle=False, opaque_t=False, comp=1):
            g = tg.create_group(name)
            g.attrs.create("class", MHDF_SPARSE, dtype="i4")  # MOAB requires this
            if is_handle:                                 # spec: values are entity IDs
                g.attrs.create("is_handle", 1, dtype="i4")
            g.create_dataset("id_list", data=np.asarray(ids, dtype="i8"))
            if opaque_t:
                g["type"] = opaque                            # committed opaque type
                arr = np.zeros(len(values), dtype=f"V{CAT_SIZE}")
                for i, b in enumerate(values):                # b is padded 32 bytes
                    arr[i] = np.void(b)
                g.create_dataset("values", data=arr)          # void -> HDF5 opaque
            elif comp > 1:
                # multi-component tag: 1-D dataspace of an HDF5 ARRAY datatype
                # (N values per entity), NOT a 2-D dataset -- MOAB requires this
                arr_dt = np.dtype((dtype, (comp,)))
                g["type"] = arr_dt
                ds = g.create_dataset("values", shape=(len(ids),), dtype=arr_dt)
                ds[...] = np.asarray(values, dtype=dtype).reshape(len(ids), comp)
            else:
                g["type"] = np.dtype(dtype)                   # committed scalar type
                g.create_dataset("values", data=np.asarray(values, dtype=dtype))

        sparse_tag(GEOM_DIMENSION, geom_ids, geom_dim, "i4")
        sparse_tag(GLOBAL_ID, gid_ids, gid_val, "i4")
        sparse_tag(CATEGORY, cat_ids, cat_val, None, opaque_t=True)
        sparse_tag(GEOM_SENSE_2, sense_ids, sense_val, "u8", is_handle=True, comp=2)
        sparse_tag(NAME, name_ids, name_val, None, opaque_t=True)

    return dict(nodes=Nv, tris=Nt, surfaces=nsurf, volumes=nvol, groups=ngrp,
                max_id=MAX_ID, path=str(path))


# =============================================================================
# 2. STRUCTURAL VERIFY  --  read the file back, check internal consistency
#    (this is NOT DAGMC acceptance; it catches malformed HDF5 / bad references)
# =============================================================================
def verify_structure(path):
    errs = []
    with h5py.File(path, "r") as f:
        node = f["tstt/nodes/coordinates"]
        tri = f["tstt/elements/Tri3/connectivity"]
        slist = f["tstt/sets/list"]
        contents = f["tstt/sets/contents"][:] if "tstt/sets/contents" in f else np.array([])
        children = f["tstt/sets/children"][:] if "tstt/sets/children" in f else np.array([])
        parents = f["tstt/sets/parents"][:] if "tstt/sets/parents" in f else np.array([])
        node_start = node.attrs["start_id"]
        tri_start = tri.attrs["start_id"]
        set_start = slist.attrs["start_id"]
        Nv, Nt, Ns = len(node), len(tri), len(slist)
        max_id = f["tstt"].attrs["max_id"]

        # id blocks must be contiguous and non-overlapping
        if tri_start != node_start + Nv:
            errs.append(f"tri start_id {tri_start} != {node_start}+{Nv}")
        if set_start != tri_start + Nt:
            errs.append(f"set start_id {set_start} != {tri_start}+{Nt}")
        if max_id != set_start + Ns - 1:
            errs.append(f"max_id {max_id} != {set_start}+{Ns}-1")

        # connectivity references valid node ids
        c = tri[:]
        if c.min() < node_start or c.max() > node_start + Nv - 1:
            errs.append(f"connectivity refs {c.min()}..{c.max()} outside node ids")

        # set list end-indices store the LAST-element index (len-1); monotonic
        if Ns:
            if slist[-1, 0] != len(contents) - 1:
                errs.append(f"contents end {slist[-1,0]} != len-1 {len(contents)-1}")
            if slist[-1, 1] != len(children) - 1:
                errs.append(f"children end {slist[-1,1]} != len-1 {len(children)-1}")
            if slist[-1, 2] != len(parents) - 1:
                errs.append(f"parents end {slist[-1,2]} != len-1 {len(parents)-1}")
            if np.any(np.diff(slist[:, 0]) < 0):
                errs.append("contents end-index not monotonic")

        # every referenced handle is a real entity id (or 0 in sense)
        all_ids = set(range(node_start, node_start + Nv + Nt + Ns))
        for nm, arr in [("contents", contents), ("children", children), ("parents", parents)]:
            bad = [int(x) for x in np.unique(arr) if int(x) not in all_ids]
            if bad:
                errs.append(f"{nm} references non-entity ids {bad[:5]}")

        # tags present, id_list/values same length
        for nm in (GEOM_DIMENSION, GLOBAL_ID, CATEGORY, GEOM_SENSE_2, NAME):
            g = f.get(f"tstt/tags/{nm}")
            if g is None:
                errs.append(f"missing tag {nm}")
                continue
            il, vv = g["id_list"], g["values"]
            if len(il) != len(vv):
                errs.append(f"tag {nm}: id_list {len(il)} != values {len(vv)}")
            bad = [int(x) for x in np.unique(il[:]) if int(x) not in all_ids]
            if bad:
                errs.append(f"tag {nm}: id_list refs non-entities {bad[:5]}")
        sense = f["tstt/tags/GEOM_SENSE_2/values"][:]
        senseflat = set(int(x) for x in np.unique(sense) if x != 0)
        vol_ids = set(range(set_start, set_start + Ns))
        if not senseflat <= vol_ids:
            errs.append(f"GEOM_SENSE_2 refs non-set ids {sorted(senseflat - vol_ids)[:5]}")

    return errs


# =============================================================================
# 3. TOPOLOGY ASSEMBLY  --  nested boundary surfaces -> surfaces/volumes/groups
# =============================================================================
def assemble(boundaries, mat_names):
    """boundaries : list of (V,T) per boundary surface, INNERMOST first, with
                    outward-oriented triangles (normals point radially out).
       mat_names  : material name per shell (len == len(boundaries)-1).
    Returns node_xyz, tri_conn, surfaces, volumes, groups for write_h5m.
    """
    nb = len(boundaries)
    assert nb == len(mat_names) + 1, "need one more boundary than materials"

    node_blocks, tri_blocks, surf_tris = [], [], []
    voff = 0
    toff = 0
    for (V, T) in boundaries:
        node_blocks.append(V)
        tri_blocks.append(T + voff)                       # global node indices
        surf_tris.append(np.arange(toff, toff + len(T)))  # global tri rows
        voff += len(V)
        toff += len(T)
    node_xyz = np.vstack(node_blocks)
    tri_conn = np.vstack(tri_blocks)

    surfaces, volumes, groups = [], [], []
    for j in range(nb):                                   # one surface per boundary
        # normals point radially OUTWARD. DAGMC sense = [forward, reverse] where
        # forward is the volume the normal points OUT of (the INNER shell j-1) and
        # reverse is the volume it points INTO (the OUTER shell j). None = implicit
        # complement (plasma void inside surface 0; exterior outside the last one).
        fwd = j - 1 if j - 1 >= 0 else None   # inner side  -> forward
        rev = j if j < nb - 1 else None       # outer side  -> reverse
        surfaces.append({"tris": surf_tris[j], "gid": j + 1, "sense": (fwd, rev)})
    for i in range(nb - 1):                               # one volume per shell
        volumes.append({"gid": i + 1, "surfs": [i, i + 1]})   # inner, outer
        groups.append({"gid": i + 1, "name": f"mat:{mat_names[i]}", "vols": [i]})
    return node_xyz, tri_conn, surfaces, volumes, groups


# =============================================================================
# 4. REAL GEOMETRY
# =============================================================================
def _read_stl_tris(path):
    """binary STL -> (Ntri,3,3) float64 array of triangle vertex coordinates."""
    import struct
    with open(path, "rb") as fh:
        fh.read(80)
        n = struct.unpack("<I", fh.read(4))[0]
        data = fh.read(n * 50)
    a = np.frombuffer(data, dtype=np.uint8).reshape(n, 50)
    return a[:, 12:48].copy().view("<f4").reshape(n, 3, 3).astype("f8")


def _dedup_surface(tri_xyz, decimals=4):
    """(m,3,3) explicit-vertex triangles -> (V (k,3) unique verts, conn (m,3) int).
    STL stores every triangle's 3 vertices explicitly; a structured surface repeats
    each vertex ~6x. Merge by rounded coordinate (mesh spacing >> 1e-4 cm here)."""
    flat = tri_xyz.reshape(-1, 3)
    _, first, inv = np.unique(np.round(flat, decimals), axis=0,
                              return_index=True, return_inverse=True)
    inv = inv.reshape(-1)
    V = np.zeros((first.size, 3))
    cnt = np.zeros(first.size)
    np.add.at(V, inv, flat)                       # average originals per unique vertex
    np.add.at(cnt, inv, 1)
    V /= cnt[:, None]
    return V, inv.reshape(-1, 3).astype("i8")


def _orient_outward(VT, sg):
    """Flip triangle winding if the closed surface's normals point inward, so every
    surface carries OUTWARD (radially out) normals -> the sense convention holds."""
    V, T = VT
    if sg.signed_volume(V, T) < 0:
        T = T[:, [0, 2, 1]]
    return V, T


def build_from_stls(geom_dir, out):
    """Build the merged DAGMC .h5m directly from the pipeline's shell STLs (ground
    truth). Each shell STL = [outer-surface tris | inner-surface tris]; the outer of
    shell k coincides with the inner of shell k+1, so we take each interface ONCE:
    boundary[0]   = shell0 inner half
    boundary[k]   = shell(k-1) outer half   (k = 1..N)  -> N+1 shared boundaries."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import stellarator_geometry as sg           # signed_volume only

    geom_dir = Path(geom_dir)
    man = json.loads((geom_dir / "manifest.json").read_text())
    layers = [L for L in man["layers"] if L["name"] not in {"sol"}]  # radial order
    mat_names = [L["name"] for L in layers]
    nb = len(layers) + 1

    boundaries = [None] * nb
    for k, L in enumerate(layers):
        tris = _read_stl_tris(geom_dir / L["stl"])
        if len(tris) % 2:
            raise SystemExit(f"{L['stl']}: odd triangle count {len(tris)}; not a "
                             "[outer|inner] shell as write_stl_shell emits")
        m = len(tris) // 2
        outer, inner = tris[:m], tris[m:]
        if k == 0:
            boundaries[0] = _orient_outward(_dedup_surface(inner), sg)
        boundaries[k + 1] = _orient_outward(_dedup_surface(outer), sg)

    node_xyz, tri_conn, surfaces, volumes, groups = assemble(boundaries, mat_names)

    # sanity: shell volume from these boundaries must match the STL shell volume
    print("[check] shell volume: assembled boundaries vs STL shell")
    worst = 0.0
    for i, L in enumerate(layers):
        Vin, Tin = boundaries[i]
        Vout, Tout = boundaries[i + 1]
        recon = abs(sg.signed_volume(Vout, Tout)) - abs(sg.signed_volume(Vin, Tin))
        ref = _stl_volume(geom_dir / L["stl"])
        rel = abs(recon - ref) / abs(ref)
        worst = max(worst, rel)
        print(f"  {L['name']:8s} boundaries={recon:13.4e}  STL={ref:13.4e}  "
              f"{rel:6.2%} {'OK' if rel < 0.01 else 'CHECK'}")
    if worst > 0.01:
        print(f"[check] WARNING: up to {worst:.1%} gap; dedup tolerance or the "
              "outer/inner split may be off. Inspect before trusting transport.")

    info = write_h5m(out, node_xyz, tri_conn, surfaces, volumes, groups)
    return info


# ---- reconstruction path (kept as fallback; see build_layers convention) -----
def build_from_geometry(stem, geom_dir, out):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import stellarator_geometry as sg   # reuse the validated mesh routines

    R, Z, nfp = sg.load_surface(stem)
    nt, nph = R.shape
    phi = 2.0 * np.pi * np.arange(nph) / nph          # full torus, periodic
    nR, nZ = sg.poloidal_outward_normals(R, Z)

    # cumulative offsets from the LCFS; 'sol' contributes its gap but is not a shell
    cum, acc = [], 0.0
    for _, thick in sg.DEFAULT_LAYERS:
        acc += thick
        cum.append(acc)                               # boundary AFTER each layer
    mat_names = [n for n, _ in sg.DEFAULT_LAYERS if n not in {"sol"}]
    # boundaries bounding the material shells: skip the LCFS itself (inside sol),
    # keep the boundary at the end of 'sol' (= W inner) through 'coil' outer
    keep = cum[len(sg.DEFAULT_LAYERS) - len(mat_names) - 1:]  # 7 boundaries

    boundaries = []
    for delta in keep:
        Ro, Zo = sg.offset_surface(R, Z, nR, nZ, delta)
        V = sg._verts(Ro, Zo, phi).reshape(-1, 3)
        T = sg._torus_triangles(nt, nph, flip=False)
        if sg.signed_volume(V, T) < 0:                # ensure OUTWARD normals
            T = sg._torus_triangles(nt, nph, flip=True)
        boundaries.append((V, T))

    node_xyz, tri_conn, surfaces, volumes, groups = assemble(boundaries, mat_names)

    # --- self-check vs the pipeline's own STL volumes (guards phi/offset) -----
    _selfcheck_volumes(geom_dir, boundaries, mat_names, sg)

    info = write_h5m(out, node_xyz, tri_conn, surfaces, volumes, groups)
    return info


def _stl_volume(path):
    """abs(signed volume) enclosed by a binary-STL closed mesh (fast, vectorized).
    For a shell-solid STL this is the annular shell volume directly."""
    import struct
    with open(path, "rb") as fh:
        fh.read(80)
        n = struct.unpack("<I", fh.read(4))[0]
        data = fh.read(n * 50)
    a = np.frombuffer(data, dtype=np.uint8).reshape(n, 50)
    tri = a[:, 12:48].copy().view("<f4").reshape(n, 3, 3).astype("f8")  # 3 verts/tri
    v0, v1, v2 = tri[:, 0], tri[:, 1], tri[:, 2]
    return abs(float(np.sum(np.einsum("ij,ij->i", v0, np.cross(v1, v2))) / 6.0))


def _selfcheck_volumes(geom_dir, boundaries, mat_names, sg, tol=0.02):
    """Cross-check reconstructed shell volumes against the pipeline's own STL
    shells (guaranteed on disk). If phi_axis or the offset convention were wrong,
    these disagree and we abort BEFORE trusting the mesh. Falls back to the
    manifest 'volume' field only where an STL is missing."""
    geom_dir = Path(geom_dir)
    man = geom_dir / "manifest.json"
    layers = {}
    if man.exists():
        layers = {L["name"]: L for L in json.loads(man.read_text())["layers"]}
    print("[selfcheck] shell volume: reconstructed vs on-disk STL")
    worst, compared = 0.0, 0
    for i, name in enumerate(mat_names):
        Vin, Tin = boundaries[i]
        Vout, Tout = boundaries[i + 1]
        recon = abs(sg.signed_volume(Vout, Tout)) - abs(sg.signed_volume(Vin, Tin))
        stl = layers.get(name, {}).get("stl")
        stl_path = geom_dir / stl if stl else geom_dir / f"{name}.stl"
        ref, src = None, ""
        if stl_path.exists():
            ref, src = _stl_volume(stl_path), stl_path.name
        elif layers.get(name, {}).get("volume"):
            ref, src = layers[name]["volume"], "manifest"
        if ref:
            rel = abs(recon - ref) / abs(ref)
            worst = max(worst, rel); compared += 1
            print(f"  {name:8s} recon={recon:13.4e}  {src:14s}={ref:13.4e}  "
                  f"{rel:6.1%} {'OK' if rel < tol else 'MISMATCH'}")
        else:
            print(f"  {name:8s} recon={recon:13.4e}  (no STL/manifest volume)")
    if compared == 0:
        print("[selfcheck] WARNING: nothing to compare against; phi/offset UNVERIFIED")
    elif worst > tol:
        raise SystemExit(
            f"[selfcheck] reconstructed shell volumes differ from the STLs by up to "
            f"{worst:.1%} (> {tol:.0%}). The phi grid or offset convention doesn't match "
            "build_layers, so this mesh is geometrically wrong even if DAGMC would load "
            "it. Send me the build_layers function (stellarator_geometry.py lines "
            "~174-245) and I'll correct the reconstruction.")


# =============================================================================
# 5. SELFTEST  --  synthetic 2-shell torus (no repo deps); write + verify
# =============================================================================
def _torus(Rmaj, rmin, nt=24, nph=48, flip=False):
    th = 2 * np.pi * np.arange(nt) / nt
    ph = 2 * np.pi * np.arange(nph) / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    X = (Rmaj + rmin * np.cos(TH)) * np.cos(PH)
    Y = (Rmaj + rmin * np.cos(TH)) * np.sin(PH)
    Zc = rmin * np.sin(TH)
    V = np.stack([X, Y, Zc], axis=-1).reshape(-1, 3)
    tris = []
    for i in range(nt):
        for j in range(nph):
            a = i * nph + j
            b = ((i + 1) % nt) * nph + j
            c = ((i + 1) % nt) * nph + (j + 1) % nph
            d = i * nph + (j + 1) % nph
            if not flip:
                tris += [(a, b, c), (a, c, d)]
            else:
                tris += [(a, c, b), (a, d, c)]
    return V, np.asarray(tris, dtype="i8")


def _signed_vol(V, T):
    tri = V[T]
    return float(np.sum(np.einsum("ij,ij->i", tri[:, 0],
                 np.cross(tri[:, 1], tri[:, 2]))) / 6.0)


def selftest(out="/tmp/selftest.h5m"):
    print("== selftest: synthetic 2-shell torus (inner=matA, outer=matB) ==")
    boundaries = []
    for rmin in (20.0, 25.0, 40.0):          # 3 boundaries -> 2 shells
        V, T = _torus(100.0, rmin)
        if _signed_vol(V, T) < 0:
            V, T = _torus(100.0, rmin, flip=True)
        boundaries.append((V, T))
    node, tri, surfaces, volumes, groups = assemble(boundaries, ["matA", "matB"])
    info = write_h5m(out, node, tri, surfaces, volumes, groups)
    print(f"  wrote {out}: {info}")
    errs = verify_structure(out)
    if errs:
        print("  STRUCTURE ERRORS:")
        for e in errs:
            print("   -", e)
        return 1
    print("  structure OK: ids contiguous, refs valid, senses reference volumes,",
          "all 5 tags present and length-consistent")
    return 0


def _write_shell_stl(inner, outer, path):
    """Mimic stellarator_geometry.write_stl_shell: [outer tris | inner tris],
    outer normals out, inner normals in."""
    import struct
    (Vo, To), (Vi, Ti) = outer, inner
    V = np.vstack([Vo, Vi])
    T = np.vstack([To, Ti + len(Vo)])
    tri = V[T]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= (np.linalg.norm(n, axis=1, keepdims=True) + 1e-30)
    with open(path, "wb") as f:
        f.write(b"\0" * 80)
        f.write(struct.pack("<I", len(T)))
        for k in range(len(T)):
            f.write(struct.pack("<3f", *n[k]))
            for v in tri[k]:
                f.write(struct.pack("<3f", *v))
            f.write(struct.pack("<H", 0))


def selftest_stl(tmp="/tmp/dw_stl"):
    """End-to-end synthetic: write 2 shell STLs + manifest, build via build_from_stls,
    verify structure. Exercises the real code path with no repo dependency except a
    tiny stub stellarator_geometry providing signed_volume."""
    import types
    d = Path(tmp); d.mkdir(exist_ok=True)
    # boundaries at r=20,25,40 -> shells matA (20->25), matB (25->40)
    tori = {}
    for r in (20.0, 25.0, 40.0):
        V, T = _torus(100.0, r)
        if _signed_vol(V, T) < 0:
            V, T = _torus(100.0, r, flip=True)
        tori[r] = (V, T)
    inv = lambda VT: (VT[0], VT[1][:, [0, 2, 1]])          # inward-oriented copy
    _write_shell_stl(inv(tori[20.0]), tori[25.0], d / "matA.stl")
    _write_shell_stl(inv(tori[25.0]), tori[40.0], d / "matB.stl")
    (d / "manifest.json").write_text(json.dumps({"layers": [
        {"name": "matA", "stl": "matA.stl", "watertight": True, "simple_cross_section": True},
        {"name": "matB", "stl": "matB.stl", "watertight": True, "simple_cross_section": True},
    ]}))
    # stub stellarator_geometry so build_from_stls imports cleanly here
    stub = types.ModuleType("stellarator_geometry")
    stub.signed_volume = _signed_vol
    sys.modules["stellarator_geometry"] = stub

    print("== selftest-stl: build_from_stls on synthetic shell STLs ==")
    info = build_from_stls(str(d), str(d / "out.h5m"))
    print("  wrote:", info)
    errs = verify_structure(str(d / "out.h5m"))
    if errs:
        print("  STRUCTURE ERRORS:")
        for e in errs:
            print("   -", e)
        return 1
    print("  structure OK")
    return 0


# =============================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true", help="serializer structure test")
    ap.add_argument("--selftest-stl", action="store_true",
                    help="end-to-end STL->h5m test on synthetic shells")
    ap.add_argument("--load-test", metavar="H5M",
                    help="load with openmc.DAGMCUniverse (needs conda openmc)")
    ap.add_argument("--stem", default="equil_precise_qa")
    ap.add_argument("--geom-dir")
    ap.add_argument("--out", default="equil_precise_qa.h5m")
    ap.add_argument("--reconstruct", action="store_true",
                    help="use LCFS reconstruction instead of the STLs (needs build_layers match)")
    a = ap.parse_args()

    if a.selftest:
        return selftest()
    if a.selftest_stl:
        return selftest_stl()

    if a.load_test:
        import openmc
        u = openmc.DAGMCUniverse(a.load_test)
        print("DAGMC LOADED. material_names:", u.material_names)
        try:
            print("bounding_box:", u.bounding_box)
        except Exception as e:
            print("bounding_box unavailable:", e)
        return 0

    geom_dir = a.geom_dir or f"spf_prototype/data/{a.stem}_geom"
    if a.reconstruct:
        info = build_from_geometry(a.stem, geom_dir, a.out)
    else:
        info = build_from_stls(geom_dir, a.out)         # default: from the STLs
    print("[dagmc_writer] wrote", info)
    errs = verify_structure(a.out)
    if errs:
        print("STRUCTURE ERRORS (fix before load-test):")
        for e in errs:
            print("  -", e)
        return 1
    print("[dagmc_writer] structure OK. Now confirm DAGMC accepts it:")
    print(f"    python dagmc_writer.py --load-test {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())