"""Three-angle poloidal cross-sections of the conformal stellarator shield with a
NON-CONSTANT (poloidally varying) radial offset. Uses the non-folding morphological
offset (stellarator_geometry_offset2): compute nested, phase-aligned shells, then
interpolate between them per poloidal vertex so the shield standoff varies smoothly
without shell-crossing. Grayscale / smplotlib serif, legend-only, for the poster.
Honest: the varying band illustrates the tool's per-location offset capability."""
import sys
from pathlib import Path
import numpy as np
import smplotlib  # noqa: F401
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon, Patch

plt.rcParams.update({"axes.labelsize": 20, "xtick.labelsize": 15,
                     "ytick.labelsize": 15, "axes.titlesize": 22})

REPO = Path("/Users/tkiker/Documents/GitHub/openmc/spf_prototype")
sys.path.insert(0, str(REPO / "python"))
import stellarator_geometry as sg
import stellarator_geometry_offset2 as sg2

OUT = "/private/tmp/claude-501/-Users-tkiker-Documents-GitHub-openmc/638db22c-edea-421d-88fe-e7f5d5ecfd10/scratchpad/poster_figs/better_cross_sections.png"
SCALE = 10.0
M = 100.0

R, Z, nfp = sg.load_surface("equil_precise_qa")
R = R * SCALE; Z = Z * SCALE
nt, nph = R.shape
phi_axis = np.linspace(0.0, 2 * np.pi, nph, endpoint=False)
period = 2 * np.pi / nfp
cols = [int(np.argmin(np.abs(phi_axis - a))) for a in (0.0, period / 4, period / 2)]

POL = 240
d_fw, d_lo, d_hi = 40.0, 110.0, 300.0           # scaled-cm offsets
shells = sg2.offset_boundaries_shared(R, Z, [d_fw, d_lo, d_hi], POL)
Rfw, Zfw = shells[d_fw]
Rlo, Zlo = shells[d_lo]
Rhi, Zhi = shells[d_hi]

# per-poloidal-vertex blend between the lo and hi shield shells -> non-constant offset
th = np.linspace(0, 2 * np.pi, POL, endpoint=False)
frac = 0.5 * (1.0 - np.cos(th))                 # 0 at theta=0, 1 at theta=pi
frac = frac.reshape(POL, 1)
Rsh = Rlo + frac * (Rhi - Rlo)
Zsh = Zlo + frac * (Zhi - Zlo)

def closed(a):
    return np.append(a, a[:1])

fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.7))
G_SH, G_FW, G_PL = "0.80", "0.45", "0.96"
for ax, j, a in zip(axes, cols, (0.0, period / 4, period / 2)):
    for (Rb, Zb, fc, z) in [(Rsh, Zsh, G_SH, 1), (Rfw, Zfw, G_FW, 2), (R, Z, G_PL, 3)]:
        rr = closed(Rb[:, j]) / M
        zz = closed(Zb[:, j]) / M
        ax.add_patch(MplPolygon(np.column_stack([rr, zz]), closed=True,
                                facecolor=fc, edgecolor="0.2", linewidth=0.9, zorder=z))
    ax.set_aspect("equal")
    allR = np.concatenate([Rsh[:, j], R[:, j]]) / M
    allZ = np.concatenate([Zsh[:, j], Z[:, j]]) / M
    mr = 0.08 * (allR.max() - allR.min()); mz = 0.08 * (allZ.max() - allZ.min())
    ax.set_xlim(allR.min() - mr, allR.max() + mr)
    ax.set_ylim(allZ.min() - mz, allZ.max() + mz)
    ax.set_xlabel("R  [m]")
    ax.set_title(rf"$\phi = {np.degrees(a):.0f}^\circ$")
axes[0].set_ylabel("Z  [m]")
fig.tight_layout()
fig.savefig(OUT, dpi=300, bbox_inches="tight")
print("wrote", OUT, " nfp", nfp, " cols", cols)
