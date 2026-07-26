#!/usr/bin/env python
"""STAGE fetch (any env with numpy): QUASR device ID -> reactor-scaled MAKEGRID coils file
   + coil-centroids npz, the inputs ParaStell (magnets) and the adjoint (filament source) need.

Reuses sweep/quasr_loader.py (pure-JSON simsopt-serial parse; no simsopt import). The coil
polylines come out at the device's NATIVE QUASR scale; we multiply by the SAME reactor-scale
factor equil_device.py applied to the boundary (read from <equil_out>/<label>_reactor_scale.txt)
so coils and the wout share one coordinate frame. MAKEGRID units are the wout units (m); the
DAGMC build is in cm, so ParaStell's construct_magnets_from_filaments handles the m->cm at
read time exactly as for QH/QA (coils files are stored x100). We therefore write the MAKEGRID
in CENTIMETERS (x100), matching pstl_test/coils_qh convention.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_THIS = Path(__file__).resolve()
# quasr_loader lives in sweep/; it does `import quasr_geom` which lives in python/ -- add both.
for base in (_THIS.parents[2],                                     # local repo spf_prototype/
             Path.home() / "src/GitHub/openmc-spf/spf_prototype",
             Path("/ginsburg/astro/users/tjk2147/spf_work/spf_pp")):
    for sub in ("sweep", "python"):
        cand = base / sub
        if cand.exists():
            sys.path.insert(0, str(cand))


def write_makegrid(coils, nfp, path, unit_scale=100.0):
    """coils = [(poly (Ni,3) meters, current A)] -> MAKEGRID file in cm (x100)."""
    with open(path, "w") as f:
        f.write(f"periods {int(nfp)}\n")
        f.write("begin filament\nmirror NIL\n")
        for k, (poly, I) in enumerate(coils):
            P = np.asarray(poly, float) * unit_scale
            for (x, y, z) in P:
                f.write(f"{x: .10E} {y: .10E} {z: .10E} {I: .10E}\n")
            x0, y0, z0 = P[0]
            f.write(f"{x0: .10E} {y0: .10E} {z0: .10E} {0.0: .10E} {k+1} coil{k+1}\n")
        f.write("end\n")
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ID", type=int)
    ap.add_argument("outdir")
    ap.add_argument("--reactor-a", type=float, default=1.704,
                    help="target minor radius (m); scale = reactor_a / catalogue minor_radius")
    ap.add_argument("--reactor-scale", type=float, default=None,
                    help="explicit coil scale factor (overrides reactor-a)")
    ap.add_argument("--n-samples", type=int, default=256)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    # Only the serial JSON (coils) + catalogue row (nfp, minor_radius) are needed -- NOT the
    # boundary/quasr_geom -- so fetch is self-contained and runs anywhere with numpy+internet.
    import quasr_loader as ql
    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    label = f"quasr{int(a.ID)}"
    coils_path = outdir / f"{label}.coils"
    cen_path = outdir / f"{label}_coil_centroids.npz"
    scale_file = outdir / f"{label}_reactor_scale.txt"
    if coils_path.exists() and cen_path.exists() and scale_file.exists() and not a.force:
        print(f"[coils] {coils_path.name} exists, skipping", flush=True)
        return

    raw_coils = ql.load_coils(a.ID, n_samples=a.n_samples)   # [(poly,current)], serial JSON only
    row = ql._catalogue_row(a.ID) or {}
    nfp = int(float(row["nfp"])) if row.get("nfp") else None
    a_native = float(row["minor_radius"]) if row.get("minor_radius") else None
    qs_class = ql._symmetry_class(row)
    if nfp is None:
        raise SystemExit(f"[coils] device {a.ID}: no nfp in catalogue -- cannot proceed")

    # reactor scale factor (SHARED with equil via the scale file): f = reactor_a / a_native.
    # The catalogue minor_radius is the single authoritative a_native so coils and the
    # DESC-scaled wout share one frame.
    f = a.reactor_scale
    if f is None:
        if a_native and a_native > 0 and a.reactor_a > 0:
            f = a.reactor_a / a_native
        else:
            f = 1.0
            print("[coils] WARNING: no catalogue minor_radius; scale=1.0 (native).", flush=True)
    scale_file.write_text(
        f"a_native {a_native}\nreactor_a {a.reactor_a}\ncoil_scale_factor {f}\n")

    coils = [(np.asarray(poly, float) * f, float(cur)) for poly, cur in raw_coils]
    Lnat = ql.total_filament_length(raw_coils)
    print(f"[coils] device {a.ID}: nfp {nfp} class {qs_class} "
          f"ncoils {len(coils)} native_len {Lnat:.4f} scale_factor {f:.4f}", flush=True)

    write_makegrid(coils, nfp, coils_path)          # writes in cm (x100)
    centroids = np.array([c.mean(axis=0) * 100.0 for c, _ in coils])  # cm, reactor-scaled
    np.savez(cen_path, centroids=centroids, currents=np.array([I for _, I in coils]),
             nfp=nfp, scale_factor=f, symmetry_class=qs_class,
             aspect=float(row["aspect_ratio"]) if row.get("aspect_ratio") else None,
             minor_radius=a_native)
    print(f"[coils] wrote {coils_path.name} ({len(coils)} coils, cm) + {cen_path.name}", flush=True)
    print("DONE_COILS", flush=True)


if __name__ == "__main__":
    main()
