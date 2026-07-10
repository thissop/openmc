"""Shared minimalist smplotlib style + paths. ASCII-only labels (smplotlib font
renders em-dash and angle brackets as tofu; use 'mean', '-', 'perp')."""
from pathlib import Path
import smplotlib  # noqa: F401  (registers the style on import)
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
FIGS = HERE.parent / "figures"
FIGS.mkdir(exist_ok=True)

# Embodies paper/PLOTTING_STYLE.md: full SM box, inward ticks all sides, big readable
# labels/ticks, frameless legends. Do NOT strip spines or add titles/grids.
plt.rcParams.update({
    "savefig.dpi": 300,
    "axes.grid": False,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "legend.fontsize": 9,
    "legend.frameon": False,
})
BLACK, BLUE, RED = "k", "#2166ac", "#b2182b"   # Rule 2 palette
GRAY_SEQ, DIVERGENT = "gray_r", "RdBu_r"        # Rule 5 colormaps


def finish(fig, name):
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGS / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote figures/{name}.pdf")
