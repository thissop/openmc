"""Shared minimalist smplotlib style + paths. ASCII-only labels (smplotlib font
renders em-dash and angle brackets as tofu; use 'mean', '-', 'perp')."""
from pathlib import Path
import smplotlib  # noqa: F401  (registers the style on import)
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
FIGS = HERE.parent / "figures"
FIGS.mkdir(exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 200,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "font.size": 10,
})


def finish(fig, name):
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGS / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote figures/{name}.pdf")
