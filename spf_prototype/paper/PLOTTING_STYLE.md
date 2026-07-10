# Plotting style — the five rules

*A prompt for any agent making figures for me. Read it before writing the first line of plot
code. Terse on purpose (after Orwell). When in doubt, obey Rule 0: smplotlib default, and
don't fight it.*

**Rule 0 — Always `import smplotlib`, and let its SM/AASTeX look stand.** It is the house
style: serif font, full box frame, inward ticks. Do not restyle away from it.

1. **One plot, one story.** Never put more on a graph than the single point it exists to
   make. If you feel you need many series to say it, you need more than one graph, not a
   bigger legend. Three flux-surface cross-sections per field period, not ten.

2. **Black and white until colour earns its place.** Separate lines by *style* — solid,
   then dashed, then dotted — all black. Bring in blue and red **only for a genuine two-way
   contrast** (two regimes, perp vs. par, ±), never for decoration. If one series must lead,
   make it bold black and let a second be dashed black; do not mix grey with black. Almost
   every plot you will ever make needs only black, white, blue, and red.

3. **Make the axes shout and the ink whisper.** Large, readable axis labels and tick labels
   — they must stand out on a printed page. **Every axis and legend label is Title Case**
   ("Toroidal Angle", "Fusion Rate", "A (Perpendicular)") — never lowercase, never sentence
   case. Keep the full box with inward ticks. No grid unless it carries information. No
   chartjunk. (Units, math symbols, and proper names keep their own case: `[%]`, `$s$`,
   `OpenMC`, `PF$_\mathrm{geo}$`.)

4. **Never title a figure, and never draw in the panel what the caption can say.** Kill
   titles. Kill any in-plot annotation that belongs in the LaTeX caption (a stray "χ²" label
   is caption text). Keep legend entries to the shortest phrase that names the series, put
   the legend *inside* the axes where it covers no data (centre it if that is the clear
   spot), and leave it frameless.

5. **Match the mark to the data.** Histograms: black outline, no fill. 2-D maps: grayscale
   (`gray_r`) for magnitudes, blue–white–red divergent (`RdBu_r`) **only** for signed
   quantities with a meaningful zero. 3-D: wireframe mesh, not a shaded solid.

> And above all: break any of these sooner than make something ugly or unreadable.

---

## Defaults (copy this stub)

```python
import smplotlib                      # Rule 0 — house style; do not override the SM look
import matplotlib.pyplot as plt

plt.rcParams.update({
    "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,  # Rule 3: big & readable
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,        # full box with ticks on all sides
    "legend.fontsize": 9, "legend.frameon": False, # Rule 4: concise, frameless
})
# DO NOT: ax.set_title(...), ax.grid(True), ax.spines[...].set_visible(False)

BLACK, BLUE, RED = "k", "#2166ac", "#b2182b"       # Rule 2 palette (blue/red = RdBu endpoints)
GRAY_SEQ, DIVERGENT = "gray_r", "RdBu_r"           # Rule 5 colormaps

# Save vector for the paper AND a png for preview:
def save(fig, name):
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"figures/{name}.{ext}")
```

## Sizing
Target the journal column, not the screen. Nuclear Fusion single column ≈ 3.4 in wide,
double ≈ 7 in. Pick the width first, then trust Rule 3's font sizes to stay legible at that
width. Do not shrink fonts to fit more in — cut content instead (Rule 1).

## Recipes

**Line plot (≤3 series, one story).** All black; solid / dashed / dotted; frameless legend
with 1–3-word labels placed off the data.

**Stellarator cross-sections.** Exactly three toroidal cuts per field period: black solid,
black dashed, black dotted. Legend centred inside the nested curves (it overlaps nothing
there). No title; label axes `R` and `Z`.

**Histogram.** `ax.hist(x, histtype="step", color="k")` — black outline, no fill.

**2-D map.** `ax.imshow(M, cmap="gray_r", origin="lower")` for a magnitude; switch to
`cmap="RdBu_r"` with a symmetric `vmin=-v, vmax=+v` only when `M` is signed. One colorbar,
labelled; no title.

**Two-way contrast.** The only routine colour case: one thing blue (`#2166ac`), the other
red (`#b2182b`), on white. (This is the two-phase / CMB look.)

**3-D.** `ax.plot_wireframe(...)` in black on white — wireframe, never a shaded surface.

## Never
- a lowercase or sentence-case axis/legend label (always Title Case)
- a figure title; an in-panel label that the caption should carry
- grey lines beside black ones; colour used for decoration rather than a real contrast
- a boxed/opaque legend, or a legend that sits on top of the data
- stripped spines / open axes (keep the full SM box)
- rainbow or `jet`; a divergent map for unsigned data
- tiny axis or tick labels to cram more onto one panel
