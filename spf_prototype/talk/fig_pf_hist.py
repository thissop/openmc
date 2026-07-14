"""Talk fig: histogram of the spin-polarization peaking-factor benefit across the 15
QUASR devices. Point: the achievable peaking-factor reduction is weak/minor for most
stellarator designs. Y = (PF_unpol - min_a2 PF)/PF_unpol, plotted in percent."""
from pathlib import Path
import numpy as np

try:
    import smplotlib  # noqa: F401  (registers the house serif/full-box/inward-tick style)
    STYLE = "smplotlib"
except Exception as e:  # pragma: no cover
    STYLE = f"plain-matplotlib (smplotlib import failed: {e})"
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CSV = HERE.parent / "paper" / "data" / "predictor_n15.csv"
FIGS = HERE / "figs"
FIGS.mkdir(exist_ok=True)

# Columns: ID,pf_unpol,a2_opt,Y,PF_geo,standoff_peak_ratio,spike_area,shaping_amp,nematic_Q
# Two leading '#' comment lines then a header row -> skip 3 rows, index by position.
data = np.genfromtxt(CSV, delimiter=",", comments="#", skip_header=3,
                     usecols=(3, 7))
Y = data[:, 0] * 100.0           # percent
shaping = data[:, 1]
N = Y.size
median = np.median(Y)

fig, ax = plt.subplots(figsize=(5.2, 3.6))
bins = np.arange(0, 18 + 1e-9, 2)
ax.hist(Y, bins=bins, histtype="step", color="k")
ax.axvline(median, color="k", ls="--", lw=1.0, label="Median")
ax.set_xlabel("Peaking-Factor Reduction [%]")
ax.set_ylabel("Number Of Devices")
ax.legend(frameon=False)

fig.tight_layout()
for e in ("pdf", "png"):
    fig.savefig(FIGS / f"pf_hist.{e}", dpi=200, bbox_inches="tight")
plt.close(fig)

print(f"[style] {STYLE}")
print(f"N devices      = {N}")
print(f"Y%  min        = {Y.min():.4f}")
print(f"Y%  median     = {median:.4f}")
print(f"Y%  mean       = {Y.mean():.4f}")
print(f"Y%  max        = {Y.max():.4f}")
print(f"shaping_amp min= {shaping.min():.4f}")
print(f"shaping_amp max= {shaping.max():.4f}")
print(f"wrote {FIGS/'pf_hist.pdf'} and .png")
