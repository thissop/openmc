# ParaStell Geometry Demo — Fusion 360 add-in

An educational Fusion 360 add-in that builds, in one click, a CAD demonstration of the geometry ideas
behind ParaStell: **sweep one profile vs. loft many *changing* poloidal sections**, a **parametric
radial build** (offsets along the local normal), **swept coils**, and a **uniform vs. non-uniform**
shield build where breeder thickness is locally traded for shield at a **fixed outer envelope**.

## 1. Folder / file structure

```
fusion_demo/
  README.md
  ParaStellDemo/
    ParaStellDemo.manifest      # add-in manifest (id, type=addin, supportedOS)
    ParaStellDemo.py            # the whole add-in (run/stop + command + geometry)
```
The folder name and the two file basenames **must match** (`ParaStellDemo/ParaStellDemo.py` +
`ParaStellDemo.manifest`) — that is how Fusion discovers an add-in.

## 2. Install & run on macOS

1. Fusion 360 → top menu **Utilities → ADD-INS** (or press **Shift+S**) → **Add-Ins** tab.
2. Click the small **green “+”** next to “My Add-Ins” and select the `ParaStellDemo` **folder**.
   - Or copy the `ParaStellDemo` folder into Fusion's add-ins directory:
     `~/Library/Application Support/Autodesk/Autodesk Fusion 360/API/AddIns/`
     then click the refresh icon; it appears in the list.
3. Select **ParaStellDemo** → **Run** (tick “Run on Startup” if you want it always loaded).
4. Have a **Design** document active. In the **SOLID** tab, the **ADD-INS** panel (right end)
   now has **“ParaStell Geometry Demo”**. Click it.
5. A dialog appears with sensible defaults — just click **OK** to generate everything, or tweak first.
6. A summary box reports which stages succeeded. Everything lands under the component
   **`ParaStell_Geometry_Demo`** in the browser.

To remove it: Add-Ins tab → select → **Stop** (the add-in cleans up its toolbar button in `stop()`).

## 3. What gets built (browser layout)

- **Axis_Centerline**, **GuideRail_theta_***  — the toroidal magnetic-axis-like curve and four
  fixed-θ guide rails (teaching curves; the loft does **not** depend on them).
- **Simple_Swept_Plasma_Demo** — one bean profile swept along the centerline → *constant* section.
- **Uniform_Radial_Build** (sub-component) — `Plasma`, `SOL`, `First_Wall`, `Breeder`, `Shield`,
  `Vacuum_Vessel` as nested lofts, constant offsets everywhere.
- **Nonuniform_Radial_Build** (sub-component) — same, but breeder is thinned and shield thickened by
  the same amount in a smooth (φ,θ) region → same outer envelope, locally more shield.
- **Coil_00 … Coil_07** — 8 swept coil rings (2 per field period for nfp=4).

Hide/show bodies to narrate; each has a clear name and (best-effort) color.

## 4. Suggested talk track (hide/show)

1. Show **Simple_Swept_Plasma_Demo**: *“one profile + one path = a sweep; the section never changes.”*
2. Show **Plasma** in `Uniform_Radial_Build` next to it: *“a real equilibrium gives a **different**
   poloidal section at every toroidal angle — a sweep can’t make this.”*
3. Show the **Plasma_phi_***`sketches`: *“these are the VMEC-like sections ParaStell samples.”*
4. Show the lofted layers: *“each wall point is pushed out along its local normal by the cumulative
   layer thickness — that’s the ParaStell radial build — then lofted.”*
5. Show **Uniform** vs **Nonuniform**: *“same total breeder+shield, same outer envelope — but here we
   locally sacrifice breeder thickness to add shield exactly where the magnets need protection.”*
6. Show **Coil_0x**: *“coils are a circular conductor swept along a coil centerline.”*

## 5. Concept mapping (paper term → Fusion feature in this add-in)

| Paper / pipeline concept | In this add-in |
|---|---|
| VMEC-like equilibrium (φ-dependent poloidal profile) | `local_uv(theta, phi)` → the `Plasma_phi_***` closed-spline **sketches** (shape changes with φ) |
| ParaStell radial offset `o_l = Σ t_i` along the FW normal | `outward_normals()` + cumulative `layer_thicknesses()` → the `SOL/First_Wall/…` profiles |
| CadQuery **loft** across poloidal profiles | the **Loft** features → `Plasma`, layer solids |
| Coil filament **sweep** of a cross-section | `build_coils()` → **Sweep** of a circle along each `Coil_xx_centerline` |
| Single-profile sweep (the *contrast* case) | `Simple_Swept_Plasma_Demo` (**Sweep** of one profile along the axis) |
| Fixed-envelope breeder↔shield trade | `layer_thicknesses(..., nonuniform=True)`: `t_breeder−δ`, `t_shield+δ`, sum constant |

## 6. Two-file workflow you asked for

The add-in builds **both** `Uniform_Radial_Build` and `Nonuniform_Radial_Build` in one document.
For your two separate Fusion files:
- Generate once. **File → Save As** “…_uniform”, hide/delete the Nonuniform component.
- **Save As** again “…_nonuniform”, hide/delete the Uniform component.
- To loft between *many* close sections in the perturbation region (smoother demo), re-run with a
  higher **N_phi** (e.g. 24–36) and a wider window; or run only the Non-uniform stage (uncheck others).

## 7. Known Fusion-API limitations / risks (all guarded, non-fatal)

- **Closed full-torus loft** can occasionally fail to self-close on organic sections. If so, untick
  **“Close full torus”** in the dialog to build an open (one-turn) loft, or lower N_theta.
- **Boolean shells**: the ordered cut cascade is wrapped in try/except. If a cut fails, you get the
  **nested solids** instead (rendered translucent so you can still see inside) — geometry still valid.
- **Coil sweeps**: each coil is independent and guarded; a failed coil is skipped and reported, the
  rest still build. Coils here are planar “picture-frame” rings (robust), not true non-planar
  modular coils.
- **Appearance/color** is fully best-effort (library appearance names vary by install); if it can’t
  set a color it silently keeps the default — geometry is never blocked by appearance.
- **Performance**: default N_phi=12, N_theta=48 builds ~100+ sketches per radial build; expect tens of
  seconds. Raising N_phi/N_theta increases quality and time.
- Everything uses **documented** API calls only (`setByAngle`, `modelToSketchSpace`,
  `sketchFittedSplines.add`, `loftFeatures`, `sweepFeatures`, `createPath`, `combineFeatures`); no
  invented methods. Fusion’s internal unit is **cm** — the dialog takes mm/deg and converts for you.

## 8. If something errors

The execute handler shows a message box with a per-stage report (and a full traceback if the whole
command fails). Paste that back and it’s an easy fix — usually a single flag (close-torus off, or a
lower resolution) gets a stubborn loft to cooperate.
