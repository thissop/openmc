#!/usr/bin/env python
"""Load public precise-QA (Wechsung 2022) and precise-QH (Wiedman/Buller/Landreman 2023)
coil sets into one normalized representation for the shield optimizer (brief S15.3).

Formats (see shield_opt/COIL_DATA_FORMATS.md, verified against simsopt 1.10.6):
  QA  : output/<cfg>/xmin.txt = simsopt BiotSavart.x DOF vector (399 floats, order 16).
        Stores only 4 base CurveXYZFourier; expand via coils_via_symmetries(nfp=2, stellsym)
        -> 16 physical coils. Requires the fixed create_curves recipe (R0=1.1,R1=0.6,order16).
        Experiment scale (R~0.24-2.23 m) -> pass scale to reach reactor size.
  QH  : coils.curves_22_7_21 = MAKEGRID Cartesian polylines, ALL 20 coils explicit
        (nfp=4, no stellsym doubling -> 20, NOT 40). Reactor scale (R~5.1-22.6 m).

QA uses simsopt (importable) for the safe symmetry expansion; QH is pure numpy.

LICENSES: QA repo has NO license (treat all-rights-reserved; cite PNAS, do not redistribute
raw files). QH is CC-BY-4.0 (attribution).
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RAW = HERE / "data" / "raw"
PROC = HERE / "data" / "processed"

QA_XMIN = {
    "QA24": RAW / "CoilsForPreciseQS" / "qs_xmin" / "xmin_QA24.txt",
    "QAWell24": RAW / "CoilsForPreciseQS" / "qs_xmin" / "xmin_QAWell24.txt",
    "QA18": RAW / "CoilsForPreciseQS" / "qs_xmin" / "xmin_QA18.txt",
    "QAWell18": RAW / "CoilsForPreciseQS" / "qs_xmin" / "xmin_QAWell18.txt",
}
QH_COILS = RAW / "QH_unzipped" / "zenodo" / "configurations" / "LandremanPaulQH_coils" / "coils.curves_22_7_21"

# Reactor scaling: fusion-relevant SPF plasmas must be reactor scale. Common reference =
# plasma minor radius a = 1.7 m (Wiedman/Buller/Landreman QH reactor scaling: a=1.7 m,
# R0~14 m, <B>=5.86 T). Native minor radii (a, meters) fix each config's scale factor;
# aspect ratios are preserved (QA A=6 -> R0=10.2 m; QH A=8 -> R0=13.7 m).
REACTOR_A_M = 1.7
NATIVE_A = {
    # from input.LandremanPaul2021_QA VMEC (R0=1.010, a=0.1683, aspect=6.000)
    "QA": 0.1683,
    # Wiedman QH is distributed already at reactor scale (a~1.7 m); factor ~1.0
    "QH": 1.70,
}


@dataclass
class CoilSet:
    name: str
    centerlines: list[np.ndarray]        # each (Ni,3), meters (already scaled)
    currents: np.ndarray                 # (n_coils,) amps
    nfp: int
    symmetry_class: str                  # 'QA' | 'QH'
    scale: float = 1.0                   # native->output multiplier applied
    meta: dict = field(default_factory=dict)

    @property
    def n_coils(self) -> int:
        return len(self.centerlines)

    def closure_gap(self) -> float:
        """Max |first - last| over coils, in output meters (0 for a closed loop)."""
        return max(float(np.linalg.norm(c[0] - c[-1])) for c in self.centerlines)

    def extent(self):
        P = np.vstack(self.centerlines)
        R = np.hypot(P[:, 0], P[:, 1])
        return dict(R_min=float(R.min()), R_max=float(R.max()),
                    Z_min=float(P[:, 2].min()), Z_max=float(P[:, 2].max()),
                    R0=float(0.5 * (R.min() + R.max())))

    def save_npz(self, path):
        # object array since coils may have unequal point counts
        arr = np.empty(self.n_coils, dtype=object)
        for i, c in enumerate(self.centerlines):
            arr[i] = c
        np.savez(path, name=self.name, nfp=self.nfp, symmetry_class=self.symmetry_class,
                 scale=self.scale, currents=self.currents, centerlines=arr,
                 meta=np.array(self.meta, dtype=object))
        return path

    def to_makegrid(self, path):
        """Write a MAKEGRID coils. file (the format ParaStell / MAKEGRID consume)."""
        with open(path, "w") as f:
            f.write(f"periods {self.nfp}\n")
            f.write("begin filament\n")
            f.write("mirror NIL\n")
            for k, (c, I) in enumerate(zip(self.centerlines, self.currents)):
                for j, (x, y, z) in enumerate(c):
                    f.write(f"{x: .10E} {y: .10E} {z: .10E} {I: .10E}\n")
                # closing terminator row: duplicate first point, current 0, group+label
                x0, y0, z0 = c[0]
                f.write(f"{x0: .10E} {y0: .10E} {z0: .10E} {0.0: .10E} {k+1} coil{k+1}\n")
            f.write("end\n")
        return path


# --------------------------------------------------------------------------- #
# QA  (Wechsung) -- simsopt path
# --------------------------------------------------------------------------- #
def load_qa(config="QA24", n_points=200, reactor_a_m=REACTOR_A_M, scale=None) -> CoilSet:
    from simsopt.geo import create_equally_spaced_curves
    from simsopt.field import Current, coils_via_symmetries, BiotSavart
    from simsopt.field.coil import ScaledCurrent

    if scale is None:                       # default: scale to reactor minor radius
        scale = reactor_a_m / NATIVE_A["QA"]
    xmin = np.loadtxt(QA_XMIN[config])
    assert xmin.size == 399, f"expected 399 dofs, got {xmin.size}"
    base = create_equally_spaced_curves(4, 2, stellsym=True, R0=1.1, R1=0.6,
                                        order=16, numquadpoints=n_points)
    base_currents = []
    for i in range(4):
        c = Current(1.0)
        if i == 0:
            c.fix_all()                       # coil-0 current fixed (not in xmin)
        base_currents.append(ScaledCurrent(c, 1e5))
    coils = coils_via_symmetries(base, base_currents, 2, True)   # 16 physical
    bs = BiotSavart(coils)
    bs.x = xmin
    centerlines = [np.asarray(c.curve.gamma()) * scale for c in coils]
    currents = np.array([c.current.get_value() for c in coils])
    cs = CoilSet(name=f"Wechsung_{config}", centerlines=centerlines, currents=currents,
                 nfp=2, symmetry_class="QA", scale=scale,
                 meta=dict(source="CoilsForPreciseQS", license="NONE-cite-PNAS",
                           native_a_m=NATIVE_A["QA"], reactor_a_m=scale * NATIVE_A["QA"],
                           order=16))
    assert cs.n_coils == 16
    return cs


# --------------------------------------------------------------------------- #
# QH  (Wiedman) -- pure-numpy MAKEGRID parser
# --------------------------------------------------------------------------- #
def load_qh(path=QH_COILS, reactor_a_m=REACTOR_A_M, scale=None) -> CoilSet:
    if scale is None:                       # QH file already ~reactor scale (a~1.7 m)
        scale = reactor_a_m / NATIVE_A["QH"]
    lines = Path(path).read_text().splitlines()
    # header: 'periods N' / 'begin filament' / 'mirror NIL'
    nfp = int(lines[0].split()[1])
    coils, currents, cur = [], [], []
    for ln in lines[3:]:
        p = ln.split()
        if not p or p[0] == "end":
            break
        x, y, z, I = float(p[0]), float(p[1]), float(p[2]), float(p[3])
        cur.append((x, y, z))
        if len(p) > 4:                        # terminator row (group+label, current 0)
            arr = np.asarray(cur[:-1]) * scale   # drop duplicate closing point
            coils.append(arr)
            currents.append(I_prev)
            cur = []
        else:
            I_prev = I
    currents = np.asarray(currents)
    cs = CoilSet(name="Wiedman_LandremanPaulQH", centerlines=coils, currents=currents,
                 nfp=nfp, symmetry_class="QH", scale=scale,
                 meta=dict(source="Zenodo-10211349", license="CC-BY-4.0",
                           native_scale="reactor ~14 m",
                           note="20 physical coils (nfp x5), NOT 40"))
    return cs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("which", choices=["qa", "qh", "both"], default="both", nargs="?")
    ap.add_argument("--qa_config", default="QA24")
    ap.add_argument("--reactor_a", type=float, default=REACTOR_A_M,
                    help="target plasma minor radius (m); reactor scale, default 1.7")
    ap.add_argument("--save", action="store_true", help="write processed npz + makegrid")
    args = ap.parse_args()

    PROC.mkdir(parents=True, exist_ok=True)
    sets = []
    if args.which in ("qa", "both"):
        sets.append(load_qa(args.qa_config, reactor_a_m=args.reactor_a))
    if args.which in ("qh", "both"):
        sets.append(load_qh(reactor_a_m=args.reactor_a))

    for cs in sets:
        e = cs.extent()
        print(f"\n{cs.name}: {cs.n_coils} coils  nfp={cs.nfp}  class={cs.symmetry_class}")
        print(f"  |I| range: {np.abs(cs.currents).min():.4e} .. {np.abs(cs.currents).max():.4e} A")
        print(f"  R: {e['R_min']:.3f}..{e['R_max']:.3f} m  Z: {e['Z_min']:.3f}..{e['Z_max']:.3f} m"
              f"  (R0~{e['R0']:.2f} m)")
        print(f"  closure gap: {cs.closure_gap():.2e} m")
        if args.save:
            npz = PROC / f"{cs.name}.npz"; mg = PROC / f"{cs.name}.coils"
            cs.save_npz(npz); cs.to_makegrid(mg)
            print(f"  saved {npz.name}, {mg.name}")


if __name__ == "__main__":
    main()
