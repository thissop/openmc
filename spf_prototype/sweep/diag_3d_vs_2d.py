#!/usr/bin/env python
"""Diagnose the 3D leak_gate vs 2D per-phi discrepancy for the buffer offset.

For each adjacent boundary pair (computed IN MEMORY at a given pol_res, identical to
what the STL holds), compute:
  - 2D per-phi: fraction of inner vertices outside the outer polygon at their own phi
  - 3D trimesh.contains: fraction of inner vertices outside the outer 3D mesh
  - for the 3D-flagged vertices, their 2D signed distance to the outer polygon (cm):
    strongly negative => well inside 2D => 3D 'outside' is a faceting/coplanar
    artifact; near/positive => a real clearance/fold problem.
"""
import sys
from pathlib import Path
import numpy as np
from matplotlib.path import Path as MplPath
from shapely.geometry import Polygon

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "python"))
import stellarator_geometry as sg
import validate_offset2 as V
import trimesh


def boundaries(stem, scale, pol_res):
    R, Z, _ = sg.load_surface(stem); R, Z = R * scale, Z * scale
    return [V.sg2.offset_boundary(R, Z, d, pol_res)
            for d in V._material_boundary_dists(V.BASELINE)], R.shape[1]


def mesh_of(RZ, phi_axis):
    Ro, Zo = RZ
    Vv = sg._verts(Ro, Zo, phi_axis).reshape(-1, 3)
    T = sg._torus_triangles(Ro.shape[0], Ro.shape[1], flip=False)
    m = trimesh.Trimesh(vertices=Vv, faces=T, process=False)
    if m.volume < 0:
        T = sg._torus_triangles(Ro.shape[0], Ro.shape[1], flip=True)
        m = trimesh.Trimesh(vertices=Vv, faces=T, process=False)
    return m


def main(stem, scale, pol_res):
    B, nph = boundaries(stem, scale, pol_res)
    phi_axis = np.linspace(0, 2 * np.pi, nph, endpoint=False)
    print(f"{stem} scale={scale:.2f} pol_res={pol_res} nph={nph}")
    for i in range(len(B) - 1):
        Ri, Zi = B[i]; Ro, Zo = B[i + 1]
        nt = Ri.shape[0]
        # 2D per phi
        out2d = 0; tot = 0
        signed_min = []  # signed distance of inner verts to outer polygon (neg=inside)
        for j in range(nph):
            poly2 = Polygon(np.column_stack([Ro[:, j], Zo[:, j]]))
            for k in range(nt):
                p = (Ri[k, j], Zi[k, j])
                d = poly2.exterior.distance(_pt(p))
                inside = poly2.contains(_pt(p))
                sd = -d if inside else d
                signed_min.append(sd)
                if not inside:
                    out2d += 1
                tot += 1
        signed_min = np.array(signed_min)
        # 3D
        mo = mesh_of((Ro, Zo), phi_axis)
        Vi = sg._verts(Ri, Zi, phi_axis).reshape(-1, 3)
        ins3d = mo.contains(Vi)
        out3d = int((~ins3d).mean() * len(Vi))
        # of the 3D-flagged, what is their 2D signed distance?
        flagged_sd = signed_min[~ins3d]
        di, do = V._material_boundary_dists(V.BASELINE)[i], V._material_boundary_dists(V.BASELINE)[i + 1]
        msg = (f"pair d=[{di:5.1f},{do:5.1f}] 2D_out={out2d/tot:.4f} "
               f"3D_out={out3d/len(Vi):.4f} watertight={mo.is_watertight}")
        if flagged_sd.size:
            msg += (f" | 3D-flagged 2D-signed-dist cm: "
                    f"min={flagged_sd.min():.3f} max={flagged_sd.max():.3f} "
                    f"median={np.median(flagged_sd):.3f} "
                    f"(neg=inside2D)")
        print(msg, flush=True)


def _pt(xy):
    from shapely.geometry import Point
    return Point(xy[0], xy[1])


if __name__ == "__main__":
    import validate_offset2 as V
    stem = sys.argv[1]
    pr = int(sys.argv[2]) if len(sys.argv) > 2 else 192
    main(stem, V.scale_for(stem), pr)
