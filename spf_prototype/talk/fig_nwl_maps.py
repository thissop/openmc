"""Talk fig: side-by-side (theta,phi) NWL maps for QUASR 886079, unpolarized (left)
vs optimally polarized (right). Restyle of paper/src/fig_conformal_map.py with subplot
titles showing the peaking factors. NWL is exactly linear in a2: NWL(a2)=unpol+a2*(par-unpol)."""
from pathlib import Path
import numpy as np

try:
    import smplotlib  # noqa: F401  (registers the house serif/full-box/inward-tick style)
    STYLE = "smplotlib"
except Exception as e:  # pragma: no cover
    STYLE = f"plain-matplotlib (smplotlib import failed: {e})"
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "paper" / "data" / "conformal_maps" / "quasr886079_conformalmap.npz"
FIGS = HERE / "figs"
FIGS.mkdir(exist_ok=True)

f = np.load(DATA)
unpol, par, dA = f["unpol"], f["par"], f["dA"]
d = dA / dA.sum()


def PF(m):
    return m.max() / (m * d).sum()


a2 = np.linspace(-1, 1, 401)
pf = np.array([PF(unpol + t * (par - unpol)) for t in a2])
opt = a2[pf.argmin()]
pf_unpol = PF(unpol)
pf_opt = pf.min()
reduction = (pf_unpol - pf_opt) / pf_unpol

m_un = unpol / (unpol * d).sum()
m_op = unpol + opt * (par - unpol)
m_op = m_op / (m_op * d).sum()

vmax = max(m_un.max(), m_op.max())
fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.2), sharey=True)
ext = [0, 360, 0, 360]
panels = [(axes[0], m_un, f"Unpolarized (P.F. = {pf_unpol:.2f})"),
          (axes[1], m_op, f"Optimally Polarized (P.F. = {pf_opt:.2f})")]
for ax, m, ttl in panels:
    im = ax.imshow(m, origin="lower", extent=ext, aspect="auto",
                   cmap="gray_r", vmin=0, vmax=vmax)
    ax.set_xlabel("Toroidal Angle $\\phi$ [deg]")
    ax.set_title(ttl, fontsize=12)
axes[0].set_ylabel("Poloidal Angle $\\theta$ [deg]")
cb = fig.colorbar(im, ax=axes, fraction=0.046, pad=0.02)
cb.set_label("NWL / Mean(NWL)")

for e in ("pdf", "png"):
    fig.savefig(FIGS / f"nwl_maps.{e}", dpi=200, bbox_inches="tight")
plt.close(fig)

print(f"[style] {STYLE}")
print(f"pf_unpol   = {pf_unpol:.6f}")
print(f"pf_opt     = {pf_opt:.6f}")
print(f"opt (a2)   = {opt:+.6f}")
print(f"reduction  = {reduction:.6f}  ({reduction*100:.2f}%)")
print(f"wrote {FIGS/'nwl_maps.pdf'} and .png")
