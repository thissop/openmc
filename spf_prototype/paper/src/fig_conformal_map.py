"""Fig: NWL/mean(NWL) on the DAGMC conformal first wall (theta,phi), unpolarized (left)
vs peaking-optimal polarization (right), for QUASR 886079. Real free-streaming transport.
Grayscale, minimalist. NWL is exactly linear in a2: NWL(a2)=unpol+a2*(par-unpol)."""
import numpy as np
import matplotlib.pyplot as plt
from _style import DATA, finish

f = np.load(DATA / "quasr886079_conformalmap.npz")
unpol, par, dA = f["unpol"], f["par"], f["dA"]
d = dA / dA.sum()


def PF(m):
    return m.max() / (m * d).sum()


a2 = np.linspace(-1, 1, 401)
pf = np.array([PF(unpol + t * (par - unpol)) for t in a2])
opt = a2[pf.argmin()]
m_un = unpol / (unpol * d).sum()
m_op = (unpol + opt * (par - unpol)); m_op = m_op / (m_op * d).sum()

vmax = max(m_un.max(), m_op.max())
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), sharey=True)
ext = [0, 360, 0, 360]
for ax, m, ttl in [(axes[0], m_un, f"unpolarized  (PF={PF(unpol):.2f})"),
                   (axes[1], m_op, f"optimal a2={opt:+.2f}  (PF={pf.min():.2f})")]:
    im = ax.imshow(m, origin="lower", extent=ext, aspect="auto",
                   cmap="gray_r", vmin=0, vmax=vmax)
    ax.set_title(ttl, fontsize=9)
    ax.set_xlabel("toroidal phi [deg]")
axes[0].set_ylabel("poloidal theta [deg]")
cb = fig.colorbar(im, ax=axes, fraction=0.046, pad=0.02)
cb.set_label("NWL / mean(NWL)")
fig.suptitle("Conformal-wall neutron load: QUASR 886079", fontsize=10)
for ext_ in ("pdf", "png"):
    fig.savefig(f"{DATA.parent}/figures/fig_conformal_map.{ext_}", bbox_inches="tight")
plt.close(fig)
print(f"wrote figures/fig_conformal_map.pdf  (optimal a2={opt:+.2f}, PF {PF(unpol):.3f}->{pf.min():.3f})")
