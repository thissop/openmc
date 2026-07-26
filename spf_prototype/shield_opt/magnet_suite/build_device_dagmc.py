#!/usr/bin/env python
"""STAGE dagmc (pstl env): device wout + coils -> watertight DAGMC h5m via ParaStell.

Device-agnostic generalization of pstl_test/corrected/build_corrected.py. Same validated
literature-anchored 8-layer radial build (FW 3.2 / mult BE / breeder 50 / backwall 4 /
shield 40 / gap 2 / VV 25 / TS 3 cm ~= 1.29 m), but PARAMETERIZED on:
  --wout   : VMEC wout .nc (reactor-scaled; from equil_device.py)
  --coils  : MAKEGRID coils file (reactor-scaled to the SAME factor as the wout)
  --nfp    : field periods -> one-period toroidal span 360/nfp, repeat = nfp-1
  --delta-npz : optional {delta_placed|delta_uniform} field for the shield<->breeder trade
                (the adaptive-shielding kill-shot); absent -> uniform baseline.

Material tags are byte-identical to build_corrected so coil_run_v3.py's material set matches.
Idempotent: skips if dagmc_<outname>.h5m exists (unless --force).

Usage: python build_device_dagmc.py --wout W --coils C --nfp N --outname TAG --export-dir D \
                                     [--be-cm 2.0] [--delta-npz F --variant placed|uniform]
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np


def write_status(export_dir, outname, status, **extra):
    """Per-stage GO/NO-GO record the orchestrator classifies on (never crash -> always a status)."""
    rec = dict(stage="dagmc", outname=outname, status=status, **extra)
    (Path(export_dir) / f"dagmc_{outname}_build_status.json").write_text(json.dumps(rec, indent=2))
    print(f"[dagmc] STATUS={status} " + " ".join(f"{k}={v}" for k, v in extra.items()), flush=True)
    return rec


def _count_subsolids(cq_shape):
    """Number of disconnected TopoDS SOLIDs inside one CadQuery object (a split/pinched layer
    is a compound of >1 solid). Returns >=1; 1 for a normal connected solid."""
    try:
        s = cq_shape.val() if hasattr(cq_shape, "val") and not hasattr(cq_shape, "Solids") else cq_shape
        subs = s.Solids()
        return max(1, len(subs))
    except Exception:
        return 1


def export_dagmc_tag_safe(stellarator, filename, export_dir,
                          min_mesh_size=5.0, max_mesh_size=20.0, algorithm=1):
    """GUARD + RECOVER replacement for stellarator.export_cad_to_dagmc.

    ParaStell builds one material tag PER in-vessel layer (+1 per coil), but gmsh tessellates a
    geometrically DISCONNECTED layer (thin conformal offset that pinches into 2 lobes for a
    high-aspect plasma) into 2+ volumes -> len(material_tags) < n_volumes -> cad_to_dagmc raises
    'number of material_tags provided is N but number of triangle sets is M'. This is a BUILD-time
    crash (upstream of the transport watertightness gate).

    Discipline (pre-flight-then-commit): run gmsh, count volumes, and
      * if n_tags == n_volumes -> commit (no split; normal export);
      * else RECOVER -> expand the tag list by each ORIGINAL solid's disconnected-sub-solid count,
        in build order (gmsh numbers a compound's sub-solids consecutively, and ParaStell appends
        the coils LAST so the trailing volumes stay 'magnets'). Commit ONLY if the expanded list
        reconciles exactly (sum == n_volumes) AND the coil-cell count is preserved;
      * else REJECT cleanly (status='dagmc_tag_mismatch') -- never let cad_to_dagmc raise.

    Returns a status dict.
    """
    import cad_to_dagmc

    export_path = Path(export_dir) / (Path(filename).with_suffix(".h5m").name)

    # rebuild the per-original-solid list + tags exactly as build_cad_to_dagmc_model did
    ivb_solids, ivb_tags = stellarator.invessel_build.extract_solids_and_mat_tags()
    magnet_solids = list(stellarator.magnet_set.all_coil_solids)
    mtag = stellarator.magnet_set.mat_tag
    magnet_tag = (mtag[0] if isinstance(mtag, (list, tuple)) else mtag)
    solids = list(ivb_solids) + magnet_solids
    per_solid_tags = list(ivb_tags) + [magnet_tag] * len(magnet_solids)
    n_short = len(per_solid_tags)

    # tessellate (mirror ParaStell.export_cad_to_dagmc up to vertices_to_h5m)
    gmsh_obj = cad_to_dagmc.init_gmsh()
    _, volumes = cad_to_dagmc.get_volumes(gmsh_obj, stellarator._geometry, method="in memory")
    cad_to_dagmc.set_sizes_for_mesh(gmsh_obj, min_mesh_size=min_mesh_size,
                                    max_mesh_size=max_mesh_size, mesh_algorithm=algorithm)
    gmsh_obj.model.mesh.generate(dim=2)
    vertices, tris_by_solid = cad_to_dagmc.mesh_to_vertices_and_triangles(volumes)
    gmsh_obj.finalize()
    n_vol = len(tris_by_solid)

    if n_short == n_vol:
        tags = per_solid_tags                      # no split -> commit as-is
        recovered = False
    else:
        # RECOVER: expand by each solid's disconnected sub-solid count, in order.
        counts = [_count_subsolids(s) for s in solids]
        expanded = []
        for tag, c in zip(per_solid_tags, counts):
            expanded.extend([tag] * c)
        n_magnet_pieces = sum(counts[len(ivb_solids):])
        if len(expanded) == n_vol and n_magnet_pieces == len(magnet_solids):
            tags = expanded
            recovered = True
            print(f"[dagmc] RECOVERED tag/solid split: {n_short} tags -> {len(expanded)} "
                  f"(n_vol={n_vol}); split solids: "
                  f"{[(t, c) for t, c in zip(per_solid_tags, counts) if c > 1]}", flush=True)
        else:
            # sub-solid counting did not reconcile with gmsh -> clean REJECT (never crash/mis-tag)
            return write_status(export_dir, Path(filename).name.replace("dagmc_", ""),
                                "dagmc_tag_mismatch", n_tags=n_short, n_volumes=n_vol,
                                sum_subsolids=int(sum(counts)),
                                n_magnet_pieces=int(n_magnet_pieces),
                                n_coils=len(magnet_solids),
                                note="sub-solid recovery did not reconcile; rejected (no crash)")

    cad_to_dagmc.vertices_to_h5m(vertices, tris_by_solid, tags, h5m_filename=export_path)
    outname = Path(filename).name.replace("dagmc_", "")
    return write_status(export_dir, outname, "go", n_tags=len(tags), n_volumes=n_vol,
                        n_coils=len(magnet_solids), recovered_split=recovered)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wout", required=True)
    ap.add_argument("--coils", required=True)
    ap.add_argument("--nfp", type=int, required=True)
    ap.add_argument("--outname", required=True)
    ap.add_argument("--export-dir", required=True)
    ap.add_argument("--be-cm", type=float, default=2.0)
    ap.add_argument("--wall-s", type=float, default=1.08)
    ap.add_argument("--n-tor", type=int, default=13, help="toroidal grid pts over one period")
    ap.add_argument("--n-pol", type=int, default=19)
    ap.add_argument("--coil-thickness", type=float, nargs=2, default=(40.0, 50.0),
                    help="magnet cross-section (cm) width height for construct_magnets_from_filaments")
    ap.add_argument("--delta-npz", default=None,
                    help="npz with delta_placed/delta_uniform on the (n_tor,n_pol) grid")
    ap.add_argument("--variant", choices=["placed", "uniform"], default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    export_dir = args.export_dir
    os.makedirs(export_dir, exist_ok=True)
    out_h5m = Path(export_dir) / f"dagmc_{args.outname}.h5m"
    if out_h5m.exists() and not args.force:
        print(f"[dagmc] {out_h5m.name} exists, skipping (use --force)", flush=True)
        return

    import parastell.parastell as ps

    nfp = args.nfp
    period_deg = 360.0 / nfp
    toroidal_angles = list(np.linspace(0.0, period_deg, args.n_tor))
    poloidal_angles = list(np.linspace(0.0, 360.0, args.n_pol))
    repeat = nfp - 1
    ones = np.ones((len(toroidal_angles), len(poloidal_angles)))

    # corrected baseline thicknesses (cm) -- identical to build_corrected.py
    T_FW, T_MULT, T_BR0, T_BW, T_SH0, T_GAP, T_VV, T_TS = \
        3.2, args.be_cm, 50.0, 4.0, 40.0, 2.0, 25.0, 3.0

    if args.delta_npz:
        assert args.variant, "--variant required with --delta-npz"
        dd = np.load(args.delta_npz)
        delta = dd["delta_placed"] if args.variant == "placed" else dd["delta_uniform"]
        assert delta.shape == ones.shape, (delta.shape, ones.shape)
    else:
        delta = np.zeros_like(ones)
    t_shield = T_SH0 + delta
    t_breeder = T_BR0 - delta
    assert np.all(t_breeder >= 15.0 - 1e-9), f"breeder below floor: {t_breeder.min():.3f}"
    assert np.allclose(t_shield + t_breeder, T_SH0 + T_BR0), "envelope not conserved"

    total = T_FW + T_MULT + T_BR0 + T_BW + T_SH0 + T_GAP + T_VV + T_TS
    print(f"[dagmc] outname={args.outname} nfp={nfp} period={period_deg:.1f}deg repeat={repeat}",
          flush=True)
    print(f"[dagmc] standoff (uniform envelope) = {total:.1f} cm; "
          f"delta max={delta.max():.2f} mean={delta.mean():.3f}", flush=True)

    radial_build_dict = {
        "first_wall":    {"thickness_matrix": ones * T_FW},
        "multiplier":    {"thickness_matrix": ones * T_MULT},
        "breeder":       {"thickness_matrix": t_breeder},
        "back_wall":     {"thickness_matrix": ones * T_BW},
        "shield":        {"thickness_matrix": t_shield},
        "gap":           {"thickness_matrix": ones * T_GAP},
        "vacuum_vessel": {"thickness_matrix": ones * T_VV, "mat_tag": "vac_vessel"},
        "thermal_shield":{"thickness_matrix": ones * T_TS},
    }

    stellarator = ps.Stellarator(args.wout)
    print("[dagmc] construct_invessel_build ...", flush=True)
    stellarator.construct_invessel_build(
        toroidal_angles, poloidal_angles, args.wall_s, radial_build_dict, repeat=repeat)
    stellarator.export_invessel_build_step(export_dir=export_dir)

    print("[dagmc] construct_magnets_from_filaments ...", flush=True)
    stellarator.construct_magnets_from_filaments(
        args.coils, args.coil_thickness[0], args.coil_thickness[1], 360.0, sample_mod=1)
    stellarator.export_magnets_step(export_dir=export_dir)
    print("[dagmc] n coil solids:", len(stellarator.magnet_set.all_coil_solids), flush=True)

    print("[dagmc] build_cad_to_dagmc_model ...", flush=True)
    stellarator.build_cad_to_dagmc_model()
    print("[dagmc] material tags (ParaStell short list):",
          len(stellarator._material_tags), flush=True)

    # GUARDED export: tag/solid-count mismatch (split layers) recovers or rejects cleanly,
    # and any unexpected error becomes a caught status -- never a mid-build crash / MPI_ABORT.
    try:
        st = export_dagmc_tag_safe(stellarator, filename=f"dagmc_{args.outname}",
                                   export_dir=export_dir)
    except Exception as e:
        st = write_status(export_dir, args.outname, f"error:dagmc:{type(e).__name__}",
                          message=str(e)[:200])
    if st["status"] != "go":
        print(f"[dagmc] NO-GO ({st['status']}); no DAGMC written. Device rejected cleanly.",
              flush=True)
        raise SystemExit(7)

    np.savez(Path(export_dir) / f"delta_{args.outname}.npz",
             delta=delta, t_shield=t_shield, t_breeder=t_breeder,
             variant=str(args.variant), be_cm=args.be_cm, nfp=nfp,
             toroidal_angles=np.array(toroidal_angles),
             poloidal_angles=np.array(poloidal_angles))
    print(f"[dagmc] DONE build_{args.outname} (recovered_split={st.get('recovered_split')})",
          flush=True)


if __name__ == "__main__":
    main()
