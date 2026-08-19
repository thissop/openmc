"""Teaching figures for the ParaStell presentation deck (smplotlib style, matches the check-in deck).
Two original explanatory diagrams (NOT reproductions of the paper's CAD figures):
  1. radial_build.pdf   -- the in-vessel radial build as a cumulative offset stack (Eq. 3), using
                           MY 8-layer QH build (mirrors ARIES-CS / the paper's WISTELL-D build).
  2. param_study.pdf    -- SCHEMATIC of the paper's Fig. 12: magnet nuclear heating over the
                           (total-thickness %, breeder %) design grid, TBR=1.05 feasibility contour,
                           optimum star. Illustrative trend only (not the paper's exact 99-config data).
"""
import numpy as np
import smplotlib  # noqa: F401  (sets the global style on import)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG = "/Users/tkiker/Documents/GitHub/openmc/spf_prototype/talk/figs"
KEY = "#1a599f"; WARN = "#b2182b"; GREY = "#7a7a7a"

# ---- Figure 1: radial build stack (cumulative offset from the first wall) ----
layers = [("First wall", 3.2, "#3b4a5a"), ("Multiplier (Be)", 2.0, "#6b8fb5"),
          ("Breeder (FLiBe)", 50.0, KEY), ("Back wall", 4.0, "#6b8fb5"),
          ("Shield (WC)", 40.0, "#9c3a3a"), ("Gap", 2.0, "#c9c9c9"),
          ("Vac. vessel", 25.0, "#4a4a4a"), ("Thermal shield", 3.0, "#8a8a8a")]
fig, ax = plt.subplots(figsize=(9.6, 2.7))
x = 0.0
for name, t, c in layers:
    ax.barh(0, t, left=x, height=0.6, color=c, edgecolor="k", linewidth=1.0)
    if t >= 3.5:
        ax.text(x + t/2, 0, f"{name}\n{t:.0f} cm", ha="center", va="center",
                color="white" if c not in ("#c9c9c9", "#8a8a8a") else "k", fontsize=8.5)
    x += t
ax.annotate("", xy=(x, 0.55), xytext=(0, 0.55),
            arrowprops=dict(arrowstyle="<->", color="k", lw=1.2))
ax.text(x/2, 0.72, f"cumulative offset  $o_l(\\phi,\\theta)=\\sum_{{i=1}}^{{l}} t_i$   "
        f"= {x:.1f} cm  (first wall $\\to$ magnet)", ha="center", va="bottom", fontsize=9)
ax.text(-1, 0, "plasma /\nSOL", ha="right", va="center", fontsize=8.5, color=GREY)
ax.set_xlim(-14, x + 2); ax.set_ylim(-0.5, 1.0); ax.axis("off")
ax.set_title("In-vessel radial build: each layer is a flux-surface offset (Eq. 3)", fontsize=11)
fig.tight_layout(); fig.savefig(f"{FIG}/radial_build.pdf"); plt.close(fig)

# ---- Figure 2: schematic parametric study (paper Fig. 12) ----
thick = np.linspace(0, 100, 11)      # total FWBS thickness %  (0..100)
breed = np.linspace(10, 90, 9)       # breeder %  (10..90)
TH, BR = np.meshgrid(thick, breed)
# heating DOWN with more thickness and less breeder (more HTS) -- the paper's qualitative trend
heat = 152.0 * (1.0 + 2.6*(1 - TH/100)) * (0.55 + 0.55*(BR/90))
# TBR UP with more breeder and more thickness; feasibility TBR>=1.05
tbr = 0.45 + 0.62*(BR/90) + 0.28*(TH/100)
fig, ax = plt.subplots(figsize=(7.4, 4.6))
pc = ax.pcolormesh(TH, BR, heat, cmap="magma_r", shading="auto")
cs = ax.contour(TH, BR, tbr, levels=[1.05], colors=[KEY], linewidths=2.2)
ax.clabel(cs, fmt="TBR = 1.05", fontsize=9)
ax.plot(100, 60, marker="*", markersize=22, color="white", markeredgecolor="k",
        markeredgewidth=1.2, zorder=5)
ax.text(100, 60, "  optimum\n  TBR 1.10\n  152 kW", ha="right", va="center", color="white", fontsize=8.5)
ax.set_xlabel("total FWBS thickness used (%)"); ax.set_ylabel("fraction to breeder (%)")
ax.set_title("Paper Fig. 12 (schematic): magnet nuclear heating vs. the design grid")
cb = fig.colorbar(pc, ax=ax); cb.set_label("total magnet nuclear heating (kW)")
fig.text(0.01, 0.01, "schematic trend, not the paper's exact 99-config data", fontsize=7, color=GREY)
fig.tight_layout(); fig.savefig(f"{FIG}/param_study.pdf"); plt.close(fig)

print("WROTE radial_build.pdf, param_study.pdf")
