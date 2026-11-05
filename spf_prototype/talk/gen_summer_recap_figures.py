"""Generate the evidence figures and engineering illustrations for the 2026 summer recap.

Every quantitative panel has one declared job in the talk. Engineering views use
the same monochrome technical language as the plots: visible edges, restrained
perspective, and no decorative lighting.
"""
import smplotlib  # noqa: F401 -- plotting protocol: first import statement

from pathlib import Path
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Patch, Polygon, Rectangle
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
import numpy as np
from PIL import Image, ImageEnhance, ImageOps
from scipy.io import netcdf_file
from shapely.geometry import Polygon as ShapelyPolygon


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIGS = HERE / "figs" / "summer_recap_2026"
DATA = ROOT / "data"
PAPER_DATA = ROOT / "paper" / "data"
RESPONSE_DATA = ROOT / "shield_opt" / "patch_worth" / "analysis_bfd_20260824"
PYTHON = ROOT / "python"
WOUT_QH = HERE / "geom" / "wout_qh.nc"
COILS_QH = HERE / "geom" / "coils_qh"
PAPER_SCREENSHOTS = HERE / "figs" / "paper_screenshots"
FIGS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(PYTHON))

from stellarator_source import ConformalStellaratorSampler  # noqa: E402
import spf_mirror as spf  # noqa: E402


plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 300,
    "axes.grid": False,
    "axes.edgecolor": "black",
    "axes.linewidth": 0.9,
    "axes.labelsize": 10,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "lines.linewidth": 1.2,
    "text.color": "black",
})


def finish(fig, name, pad=0.02, tight=True):
    if tight:
        fig.tight_layout(pad=0.55)
    fig.savefig(FIGS / f"{name}.pdf", bbox_inches="tight", pad_inches=pad)
    fig.savefig(FIGS / f"{name}.png", bbox_inches="tight", pad_inches=pad)
    plt.close(fig)


def load_vmec_surface(path=WOUT_QH, n_theta=121, n_phi=241, phi_start=0.0,
                      phi_stop=2.0 * np.pi, radial_index=-1):
    """Return one QH VMEC flux surface on a regular grid in metres."""
    handle = netcdf_file(str(path), mmap=False)
    xm = np.asarray(handle.variables["xm"].data).copy()
    xn = np.asarray(handle.variables["xn"].data).copy()
    rmnc = np.asarray(handle.variables["rmnc"].data[radial_index]).copy()
    zmns = np.asarray(handle.variables["zmns"].data[radial_index]).copy()
    rmns = (np.asarray(handle.variables["rmns"].data[radial_index]).copy()
            if "rmns" in handle.variables else None)
    zmnc = (np.asarray(handle.variables["zmnc"].data[radial_index]).copy()
            if "zmnc" in handle.variables else None)
    nfp = int(handle.variables["nfp"].data)
    handle.close()

    theta = np.linspace(0.0, 2.0 * np.pi, n_theta)
    phi = np.linspace(phi_start, phi_stop, n_phi)
    th, ph = np.meshgrid(theta, phi, indexing="ij")
    phase = xm[None, None, :] * th[..., None] - xn[None, None, :] * ph[..., None]
    radius = np.sum(rmnc[None, None, :] * np.cos(phase), axis=-1)
    height = np.sum(zmns[None, None, :] * np.sin(phase), axis=-1)
    if rmns is not None:
        radius += np.sum(rmns[None, None, :] * np.sin(phase), axis=-1)
    if zmnc is not None:
        height += np.sum(zmnc[None, None, :] * np.cos(phase), axis=-1)
    x = radius * np.cos(ph)
    y = radius * np.sin(ph)
    return dict(theta=theta, phi=phi, R=radius, Z=height, X=x, Y=y, nfp=nfp)


def moderated_stellarator_surface(n_theta=81, n_phi=181, blend=0.62):
    """Return the QH VMEC surface with its three-dimensional shaping moderated.

    Blending the real boundary toward its toroidal mean retains the field-period
    structure and rotating cross-section without the visually extreme torsion of
    the full-strength reactor surface.
    """
    surface = load_vmec_surface(n_theta=n_theta, n_phi=n_phi)
    radius_mean = surface["R"].mean(axis=1, keepdims=True)
    height_mean = surface["Z"].mean(axis=1, keepdims=True)
    radius = radius_mean + blend * (surface["R"] - radius_mean)
    height = height_mean + blend * (surface["Z"] - height_mean)
    phi = surface["phi"][None, :]
    x = radius * np.cos(phi)
    y = radius * np.sin(phi)
    return x, y, height


def parse_makegrid(path=COILS_QH):
    """Read full-device MAKEGRID filaments as Cartesian polylines."""
    coils = []
    current = []
    for line in path.read_text().splitlines():
        words = line.split()
        if len(words) < 4 or words[0] in {"periods", "begin", "mirror", "end"}:
            continue
        try:
            point = [float(words[0]), float(words[1]), float(words[2])]
            amperes = float(words[3])
        except ValueError:
            continue
        current.append(point)
        if amperes == 0.0 or len(words) > 4:
            coils.append(np.asarray(current))
            current = []
    return coils


def equalize_3d(ax, arrays, z_scale=1.0, zoom=1.25):
    points = np.vstack([np.column_stack(a) for a in arrays])
    low = points.min(axis=0)
    high = points.max(axis=0)
    center = 0.5 * (low + high)
    span = high - low
    ax.set_xlim(center[0] - span[0] / 2, center[0] + span[0] / 2)
    ax.set_ylim(center[1] - span[1] / 2, center[1] + span[1] / 2)
    ax.set_zlim(center[2] - z_scale * span[2] / 2,
                center[2] + z_scale * span[2] / 2)
    try:
        ax.set_box_aspect((span[0], span[1], z_scale * span[2]), zoom=zoom)
    except TypeError:
        ax.set_box_aspect((span[0], span[1], z_scale * span[2]))


def qh_reactor_bw():
    """Actual QH LCFS and coil filaments: the opening engineering object."""
    surface = load_vmec_surface(phi_start=np.deg2rad(42.0), phi_stop=2.0 * np.pi)
    coils = parse_makegrid()
    fig = plt.figure(figsize=(9.0, 5.2))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_wireframe(surface["X"], surface["Y"], surface["Z"],
                      rstride=6, cstride=6, color="black", linewidth=0.40,
                      alpha=0.82)
    for coil in coils:
        ax.plot(coil[:, 0], coil[:, 1], coil[:, 2], color="0.30",
                linewidth=0.72, alpha=0.88)
    equalize_3d(ax, [(surface["X"].ravel(), surface["Y"].ravel(),
                      surface["Z"].ravel())] +
                    [(coil[:, 0], coil[:, 1], coil[:, 2]) for coil in coils],
                z_scale=1.18, zoom=1.42)
    ax.view_init(elev=24, azim=-58)
    ax.set_axis_off()
    fig.subplots_adjust(0, 0, 1, 1)
    finish(fig, "qh_reactor_bw", pad=0.0, tight=False)


def qh_plasma_bw():
    """Actual QH last closed flux surface, without coils, for the title slide."""
    surface = load_vmec_surface(phi_start=np.deg2rad(34.0), phi_stop=2.0 * np.pi)
    fig = plt.figure(figsize=(9.0, 5.0))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_wireframe(surface["X"], surface["Y"], surface["Z"],
                      rstride=6, cstride=6, color="black", linewidth=0.45,
                      alpha=0.84)
    equalize_3d(ax, [(surface["X"].ravel(), surface["Y"].ravel(),
                      surface["Z"].ravel())], z_scale=1.20, zoom=1.48)
    ax.view_init(elev=24, azim=-58)
    ax.set_axis_off()
    fig.subplots_adjust(0, 0, 1, 1)
    finish(fig, "qh_plasma_bw", pad=0.0, tight=False)


def title_stellarator_bw():
    """Moderately shaped plasma-only stellarator mesh for title and summary."""
    x, y, z = moderated_stellarator_surface(n_theta=81, n_phi=201)
    fig = plt.figure(figsize=(8.5, 4.8))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_wireframe(x, y, z, rstride=6, cstride=9, color="black",
                      linewidth=0.43, alpha=0.84)
    equalize_3d(ax, [(x.ravel(), y.ravel(), z.ravel())],
                z_scale=1.18, zoom=1.08)
    ax.view_init(elev=24, azim=-58)
    ax.set_axis_off()
    fig.subplots_adjust(0, 0, 1, 1)
    finish(fig, "title_stellarator_bw", pad=0.0, tight=False)


def tokamak_vs_stellarator_bw():
    """Axisymmetric torus and moderated QH stellarator in a common view."""
    theta = np.linspace(0, 2 * np.pi, 81)
    phi = np.linspace(0, 2 * np.pi, 161)
    th, ph = np.meshgrid(theta, phi, indexing="ij")
    radius = 10.9 + 3.1 * np.cos(th + 0.18 * np.sin(th))
    height = 4.9 * np.sin(th)
    tokamak = (radius * np.cos(ph), radius * np.sin(ph), height)
    stellarator = moderated_stellarator_surface(n_theta=81, n_phi=181)

    fig = plt.figure(figsize=(9.0, 4.0))
    axes = [fig.add_subplot(121, projection="3d"),
            fig.add_subplot(122, projection="3d")]
    for ax, arrays in zip(axes, [tokamak, stellarator], strict=True):
        ax.plot_wireframe(*arrays, rstride=6, cstride=8, color="black",
                          linewidth=0.38, alpha=0.82)
        equalize_3d(ax, [(arrays[0].ravel(), arrays[1].ravel(),
                          arrays[2].ravel())], z_scale=1.20, zoom=1.10)
        ax.view_init(elev=24, azim=-58)
        ax.set_axis_off()
    fig.subplots_adjust(0, 0, 1, 1, wspace=0.02)
    finish(fig, "tokamak_vs_stellarator_bw", pad=0.0, tight=False)


def vmec_coordinates_bw():
    """Nested surfaces from the actual QH VMEC equilibrium."""
    indices = [40, 80, 120, 160, 200]
    surfaces = [load_vmec_surface(n_theta=81, n_phi=161,
                                  radial_index=index) for index in indices]
    fig = plt.figure(figsize=(8.6, 4.6))
    ax = fig.add_subplot(111, projection="3d")
    shades = np.linspace(0.72, 0.05, len(surfaces))
    for surface, shade in zip(surfaces, shades, strict=True):
        ax.plot_wireframe(surface["X"], surface["Y"], surface["Z"],
                          rstride=10, cstride=12, color=str(shade),
                          linewidth=0.42, alpha=0.82)
    outer = surfaces[-1]
    equalize_3d(ax, [(outer["X"].ravel(), outer["Y"].ravel(),
                      outer["Z"].ravel())], z_scale=1.22, zoom=1.43)
    ax.view_init(elev=25, azim=-57)
    ax.set_axis_off()
    fig.subplots_adjust(0, 0, 1, 1)
    finish(fig, "vmec_coordinates_bw", pad=0.0, tight=False)


def source_births_bw():
    """VMEC birth sites and local field directions on the reactor-scale QH."""
    fluxmap = dict(np.load(DATA / "qh_freestream_fluxmap.npz"))
    sampler = ConformalStellaratorSampler(fluxmap, S_of_rho=lambda rho: 1 - rho ** 2)
    samples = sampler.sample(18_000, rng=np.random.default_rng(20260824))
    surface = load_vmec_surface()

    fig = plt.figure(figsize=(8.2, 4.8))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_wireframe(surface["X"], surface["Y"], surface["Z"],
                      rstride=8, cstride=8, color="0.55", linewidth=0.30,
                      alpha=0.65)
    keep = slice(None, None, 14)
    ax.scatter(samples["x"][keep], samples["y"][keep], samples["z"][keep],
               s=2.0, color="black", alpha=0.28, depthshade=False,
               rasterized=True)
    arrow_ids = [900, 5100, 11_500]
    ax.quiver(samples["x"][arrow_ids], samples["y"][arrow_ids],
              samples["z"][arrow_ids],
              samples["bhat"][arrow_ids, 0], samples["bhat"][arrow_ids, 1],
              samples["bhat"][arrow_ids, 2], length=2.5, normalize=True,
              color="black", linewidth=1.0, arrow_length_ratio=0.22)
    equalize_3d(ax, [(surface["X"].ravel(), surface["Y"].ravel(),
                      surface["Z"].ravel())], z_scale=1.20, zoom=1.42)
    ax.view_init(elev=25, azim=-57)
    ax.set_axis_off()
    fig.subplots_adjust(0, 0, 1, 1)
    finish(fig, "source_births_bw", pad=0.0, tight=False)


def birth_density():
    """The sampler's physical-space and radial birth distributions."""
    fluxmap = dict(np.load(DATA / "qh_freestream_fluxmap.npz"))
    sampler = ConformalStellaratorSampler(fluxmap, S_of_rho=lambda rho: 1 - rho ** 2)
    samples = sampler.sample(90_000, rng=np.random.default_rng(803097))

    fig, axes = plt.subplots(1, 2, figsize=(8.3, 3.35),
                             gridspec_kw={"width_ratios": [1.05, 1]})
    keep = slice(None, None, 5)
    axes[0].scatter(samples["R"][keep], samples["z"][keep], s=0.55,
                    color="black", alpha=0.12, rasterized=True)
    axes[0].set_aspect("equal")
    axes[0].set_xlabel("Major Radius (m)")
    axes[0].set_ylabel("Vertical Position (m)")

    rho = fluxmap["rho"]
    drho = rho[1] - rho[0]
    edges = np.concatenate([[rho[0] - drho / 2], rho + drho / 2])
    observed, _ = np.histogram(samples["rho"], bins=edges, density=True)
    expected = (1 - rho ** 2) * fluxmap["sqrtg"].sum(axis=(1, 2))
    expected = expected / np.sum(expected * drho)
    axes[1].step(rho, observed, where="mid", color="black", label="Sampled")
    axes[1].plot(rho, expected, color="0.55", linestyle="--",
                 label=r"$S(\rho)\langle\sqrt{g}\rangle$")
    axes[1].set_xlabel(r"Normalized Flux Radius $\rho$")
    axes[1].set_ylabel("Probability Density")
    axes[1].legend(loc="upper right")
    finish(fig, "birth_density")


def spf_geometry_bw():
    """How the polarized angular kernels steer emission relative to B."""
    theta = np.linspace(0, np.pi, 45)
    phi = np.linspace(0, 2 * np.pi, 73)
    th, ph = np.meshgrid(theta, phi, indexing="ij")
    mu = np.cos(th)
    kernels = [
        ("Unpolarized", np.full_like(mu, 0.50)),
        ("Perpendicular Mode", 0.75 * (1.0 - mu ** 2)),
        ("Parallel Mode", 0.25 + 0.75 * mu ** 2),
    ]
    fig = plt.figure(figsize=(9.2, 3.5))
    for index, (title, density) in enumerate(kernels, start=1):
        ax = fig.add_subplot(1, 3, index, projection="3d")
        radius = np.sqrt(np.maximum(density, 0.0))
        x = radius * np.sin(th) * np.cos(ph)
        y = radius * np.sin(th) * np.sin(ph)
        z = radius * np.cos(th)
        ax.plot_wireframe(x, y, z, rstride=4, cstride=6, color="black",
                          linewidth=0.38, alpha=0.82)
        ax.plot([0, 0], [0, 0], [-1.02, 1.02], color="black", linewidth=1.0)
        ax.quiver(0, 0, 0.62, 0, 0, 0.28, color="black", linewidth=1.0,
                  arrow_length_ratio=0.28)
        label_height = 1.24 if index == 3 else 1.04
        ax.text(0.05, 0, label_height, r"$\widehat{\mathbf{B}}$", fontsize=10)
        equalize_3d(ax, [(x.ravel(), y.ravel(), z.ravel())],
                    z_scale=1.0, zoom=1.12)
        if index == 3:
            low, _ = ax.get_zlim()
            ax.set_zlim(low, 1.34)
        ax.view_init(elev=17, azim=-56)
        ax.set_title(title, fontsize=14, pad=18)
        ax.set_axis_off()
    fig.subplots_adjust(0, 0, 1, 0.84, wspace=0.01)
    finish(fig, "spf_geometry_bw", pad=0.0, tight=False)


def sampler_validation():
    """Monte Carlo direction samples against the two analytic SPF kernels."""
    rng_a = spf.Rng(20260824)
    rng_b = spf.Rng(20260825)
    mode_a = spf.make_mode_weights(1.0, 0.0, 0.0)
    mode_b = spf.make_mode_weights(0.0, 1.0, 0.0)
    n = 45_000
    xa = np.array([spf.sample_local_direction(mode_a, rng_a)[2] for _ in range(n)])
    xb = np.array([spf.sample_local_direction(mode_b, rng_b)[2] for _ in range(n)])
    edges = np.linspace(-1, 1, 35)
    centers = (edges[:-1] + edges[1:]) / 2
    x = np.linspace(-1, 1, 400)

    fig, ax = plt.subplots(figsize=(5.3, 3.45))
    ha, _ = np.histogram(xa, bins=edges, density=True)
    hb, _ = np.histogram(xb, bins=edges, density=True)
    ax.plot(centers, ha, "o", ms=3.0, mfc="white", mec="black",
            label="Mode A Sample")
    ax.plot(centers, hb, "s", ms=2.8, mfc="white", mec="0.45",
            label="Mode B Sample")
    ax.plot(x, 0.75 * (1 - x ** 2), color="black",
            label=r"$\frac{3}{4}(1-\mu^2)$")
    ax.plot(x, 0.25 + 0.75 * x ** 2, color="0.45", linestyle="--",
            label=r"$\frac{1}{4}+\frac{3}{4}\mu^2$")
    ax.set_xlabel(r"Pitch Cosine $\mu=\widehat{\mathbf{u}}\!\cdot\!\widehat{\mathbf{B}}$")
    ax.set_ylabel("Probability Density")
    ax.set_xlim(-1, 1)
    ax.legend(ncol=2, loc="upper center")
    finish(fig, "sampler_validation")


def transport_validation():
    """Direct residuals for the conformal-wall transport cross-check."""
    with open(PAPER_DATA / "conformal_xcheck.csv") as handle:
        rows = list(csv.DictReader(row for row in handle if not row.startswith("#")))
    regions = ["Inboard", "Outboard", "Global Mean"]
    a_error = np.array([float(row["A_rel_pct"]) for row in rows])
    b_error = np.array([float(row["B_rel_pct"]) for row in rows])

    fig, ax = plt.subplots(figsize=(6.8, 3.5))
    x = np.arange(len(regions))
    width = 0.34
    ax.bar(x - width / 2, a_error, width, facecolor="white",
           edgecolor="black", linewidth=0.9, label="Mode A")
    ax.bar(x + width / 2, b_error, width, facecolor="white",
           edgecolor="black", linewidth=0.9, hatch="////", label="Mode B")
    ax.set_xticks(x)
    ax.set_xticklabels(regions)
    ax.set_ylabel("OpenMC Difference (percent)")
    ax.set_ylim(0, 2.8)
    ax.legend(loc="upper right", ncol=2)
    finish(fig, "transport_validation")


def rate_steering():
    """Tokamak comparison in which rate and load changes must be read together."""
    with open(PAPER_DATA / "rate_steering.csv") as handle:
        rows = {row["mode"]: row for row in csv.DictReader(
            line for line in handle if not line.startswith("#"))}
    modes = ["nonpol", "A", "B", "C"]
    labels = ["Unpolarized", "A", "B", "C"]
    metrics = [("fusion_rate", "Fusion Rate"),
               ("center_stack_load", "Center-Stack Load"),
               ("outboard_load", "Outboard Load")]
    hatches = ["", "////", "...."]
    x = np.arange(len(modes))
    width = 0.25

    fig, ax = plt.subplots(figsize=(6.1, 3.55))
    for index, ((key, label), hatch) in enumerate(zip(metrics, hatches)):
        values = [float(rows[mode][key]) for mode in modes]
        ax.bar(x + (index - 1) * width, values, width, label=label,
               facecolor="white", edgecolor="black", linewidth=0.7,
               hatch=hatch)
    ax.axhline(1.0, color="black", linestyle=":", linewidth=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Relative To Unpolarized Source")
    ax.set_ylim(0, 1.82)
    ax.legend(loc="upper right")
    finish(fig, "rate_steering")


def conformal_map():
    """The one wall-map comparison that establishes the modest stellarator effect."""
    data = np.load(PAPER_DATA / "conformal_maps" / "quasr886079_conformalmap.npz")
    unpolarized = data["unpol"]
    parallel = data["par"]
    area = data["dA"] / np.sum(data["dA"])
    maps = [unpolarized / np.sum(unpolarized * area),
            parallel / np.sum(parallel * area)]
    extent = [0, 360, 0, 360]
    norm = Normalize(vmin=min(values.min() for values in maps),
                     vmax=max(values.max() for values in maps))

    fig = plt.figure(figsize=(8.7, 3.95))
    grid = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.060], wspace=0.24)
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    color_axis = fig.add_subplot(grid[0, 2])
    for ax, values in zip(axes, maps):
        image = ax.imshow(values, origin="lower", extent=extent, aspect="auto",
                          cmap="gray_r", norm=norm, interpolation="nearest")
        ax.set_xlabel(r"Toroidal Angle $\phi$ (degrees)", fontsize=13)
        ax.tick_params(labelsize=11.5)
        ax.set_xticks([0, 120, 240, 360])
        ax.set_yticks([0, 120, 240, 360])
    axes[0].set_title("Unpolarized", fontsize=14, pad=9)
    axes[1].set_title("Parallel Mode", fontsize=14, pad=9)
    axes[0].set_ylabel(r"Poloidal Angle $\vartheta$ (degrees)", fontsize=13)
    axes[1].tick_params(labelleft=False)
    bar = fig.colorbar(image, cax=color_axis)
    bar.set_label("Wall Load / Mean Wall Load", fontsize=12.5)
    color_axis.tick_params(labelsize=11.0)
    fig.subplots_adjust(left=0.10, right=0.93, bottom=0.19, top=0.86)
    finish(fig, "conformal_map", tight=False)


def predictor_scan():
    """The powered check of whether simple geometry predicts SPF benefit."""
    with open(PAPER_DATA / "predictor_n15.csv") as handle:
        rows = list(csv.DictReader(row for row in handle if not row.startswith("#")))
    benefit = np.array([float(row["Y"]) * 100 for row in rows])
    headroom = np.array([float(row["PF_geo"]) for row in rows])
    shape = np.array([float(row["nematic_Q"]) for row in rows])

    fig, axes = plt.subplots(1, 2, figsize=(7.9, 3.35), sharey=True)
    axes[0].scatter(headroom, benefit, s=56, facecolor="white", edgecolor="black")
    axes[0].set_xlabel(r"Headroom Proxy $\mathrm{PF}_{\mathrm{geo}}$", fontsize=13.5)
    axes[0].set_ylabel("Peaking Reduction (percent)", fontsize=13.5)
    axes[1].scatter(shape, benefit, marker="s", s=52, facecolor="white",
                    edgecolor="0.40")
    axes[1].set_xlabel(r"Nematic Shape Metric $Q$", fontsize=13.5)
    for ax in axes:
        ax.axhline(0, color="0.75", linewidth=0.8)
        ax.set_ylim(-1, 18)
        ax.tick_params(labelsize=12)
    finish(fig, "predictor_scan")


def crop_to_content(image, cutoff=248, margin=8):
    gray = np.asarray(image.convert("L"))
    mask = gray < cutoff
    if not np.any(mask):
        return image
    y, x = np.where(mask)
    box = (max(0, x.min() - margin), max(0, y.min() - margin),
           min(image.width, x.max() + margin), min(image.height, y.max() + margin))
    return image.crop(box)


def config_montage_bw():
    """Crop the actual QUASR coil/plasma render from five solved configurations."""
    ids = [1960314, 59509, 803097, 886079, 932746]
    panels = []
    for device in ids:
        path = ROOT / "figs" / "vmec_devices" / f"device_{device}.png"
        image = Image.open(path).convert("L")
        crop = image.crop((25, 200, int(image.width * 0.17), image.height - 40))
        crop = crop_to_content(ImageOps.autocontrast(crop), cutoff=246, margin=10)
        scale = min(505 / crop.width, 400 / crop.height)
        crop = crop.resize((max(1, round(crop.width * scale)),
                            max(1, round(crop.height * scale))),
                           Image.Resampling.LANCZOS)
        canvas = Image.new("L", (520, 420), 255)
        canvas.paste(crop, ((520 - crop.width) // 2, (420 - crop.height) // 2))
        panels.append(canvas)
    montage = Image.new("L", (520 * len(panels), 420), 255)
    for index, panel in enumerate(panels):
        montage.paste(panel, (index * 520, 0))
    montage.save(FIGS / "config_montage_bw.png", dpi=(300, 300))


def paper_screenshots_bw():
    """Convert the four user-selected literature screenshots to monochrome."""
    sources = {
        "bae_screenshot_bw": PAPER_SCREENSHOTS / "bae_spherical_tokamak_spf.png",
        "schwartz_screenshot_bw": PAPER_SCREENSHOTS / "analytic_schwartz.png",
        "parastell_screenshot_bw": PAPER_SCREENSHOTS / "parastell.png",
        "aries_cs_screenshot_bw": PAPER_SCREENSHOTS / "aries_cs.png",
    }
    for name, path in sources.items():
        image = Image.open(path).convert("L")
        image = ImageEnhance.Contrast(ImageOps.autocontrast(image)).enhance(1.04)
        image = ImageOps.expand(image, border=3, fill=0)
        image.save(FIGS / f"{name}.png", dpi=(300, 300))


def literature_crops_bw():
    """Legible title-page fragments; these are citations, not decoration."""
    crops = {
        "parastell_cover_bw": (HERE / "figs" / "paper_parastell.png",
                                (300, 130, 1210, 650)),
        "helios_cover_bw": (HERE / "figs" / "paper_helios.png",
                             (65, 20, 1190, 510)),
    }
    for name, (path, box) in crops.items():
        image = Image.open(path).convert("L").crop(box)
        image = ImageEnhance.Contrast(ImageOps.autocontrast(image)).enhance(1.08)
        image = ImageOps.expand(image, border=2, fill=0)
        image.save(FIGS / f"{name}.png", dpi=(300, 300))


def qh_cross_section(phi_value=0.0, n_theta=241):
    surface = load_vmec_surface(n_theta=n_theta, n_phi=2,
                                phi_start=phi_value, phi_stop=phi_value)
    radius = surface["R"][:, 0]
    height = surface["Z"][:, 0]
    dr = 0.5 * (np.roll(radius, -1) - np.roll(radius, 1))
    dz = 0.5 * (np.roll(height, -1) - np.roll(height, 1))
    nr, nz = dz.copy(), -dr.copy()
    mag = np.hypot(nr, nz)
    nr, nz = nr / mag, nz / mag
    outward = (radius - radius.mean()) * nr + (height - height.mean()) * nz
    sign = np.where(outward >= 0, 1.0, -1.0)
    return radius, height, nr * sign, nz * sign


def radial_build_bw():
    """Actual QH poloidal section with the Step-1 layer dimensions."""
    radius, height, nr, nz = qh_cross_section()
    layers = [
        ("First Wall", 0.032, "", "0.88"),
        ("Be Multiplier", 0.020, "..", "0.94"),
        ("FLiBe Breeder", 0.500, "", "0.78"),
        ("Back Wall", 0.040, "xx", "0.92"),
        ("WC Shield", 0.400, "////", "0.98"),
        ("Gap", 0.020, "", "1.00"),
        ("Vacuum Vessel", 0.250, "\\\\", "0.90"),
        ("Thermal Shield", 0.030, "++", "0.96"),
    ]
    cumulative = np.cumsum([0.0] + [item[1] for item in layers])
    plasma = ShapelyPolygon(np.column_stack([radius, height])).buffer(0)

    fig, ax = plt.subplots(figsize=(5.4, 5.3))
    for index in range(len(layers) - 1, -1, -1):
        label, _, hatch, shade = layers[index]
        offset = cumulative[index + 1]
        shell = plasma.buffer(float(offset), join_style="round", resolution=20)
        rr, zz = np.asarray(shell.exterior.coords).T
        ax.add_patch(Polygon(np.column_stack([rr, zz]), closed=True,
                             facecolor=shade, edgecolor="black", linewidth=0.65,
                             hatch=hatch, label=label))
    rr, zz = np.asarray(plasma.exterior.coords).T
    ax.add_patch(Polygon(np.column_stack([rr, zz]), closed=True,
                         facecolor="white", edgecolor="black", linewidth=0.9))
    ax.set_aspect("equal")
    ax.autoscale_view()
    ax.set_xlabel("Major Radius (m)", fontsize=13.5)
    ax.set_ylabel("Vertical Position (m)", fontsize=13.5)
    ax.tick_params(labelsize=12)
    finish(fig, "radial_build_bw", pad=0.03)


def parastell_iteration_bw():
    """Three actual QH offset surfaces representing a regenerated in-vessel build."""
    x, y, z, _ = shield_surface_geometry(n_theta=49, n_phi=121,
                                         shield_mid_offset=0.79)
    plasma = load_vmec_surface(n_theta=49, n_phi=121)
    breeder_x, breeder_y, breeder_z, _ = shield_surface_geometry(
        n_theta=49, n_phi=121, shield_mid_offset=0.32)
    fig = plt.figure(figsize=(8.6, 4.8))
    ax = fig.add_subplot(111, projection="3d")
    surfaces = [
        (plasma["X"], plasma["Y"], plasma["Z"], "0.68", 0.32),
        (breeder_x, breeder_y, breeder_z, "0.35", 0.40),
        (x, y, z, "0.02", 0.52),
    ]
    for sx, sy, sz, shade, width in surfaces:
        ax.plot_wireframe(sx, sy, sz, rstride=8, cstride=8,
                          color=shade, linewidth=width, alpha=0.82)
    equalize_3d(ax, [(x.ravel(), y.ravel(), z.ravel())],
                z_scale=1.20, zoom=1.40)
    ax.view_init(elev=25, azim=-57)
    ax.set_axis_off()
    fig.subplots_adjust(0, 0, 1, 1)
    finish(fig, "parastell_iteration_bw", pad=0.0, tight=False)


def periodic_boundaries_bw():
    """Full torus and rotational periodic wedge in top view."""
    surface = load_vmec_surface(n_theta=49, n_phi=193)
    fig, axes = plt.subplots(1, 2, figsize=(7.7, 3.25))
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_axis_off()
    for row in range(0, surface["X"].shape[0], 8):
        axes[0].plot(surface["X"][row], surface["Y"][row], color="black",
                     linewidth=0.55)
    axes[0].set_title("Full Torus", fontsize=12, pad=4)

    wedge = surface["phi"] <= (np.pi / 2 + 1e-10)
    for row in range(0, surface["X"].shape[0], 8):
        axes[1].plot(surface["X"][row, wedge], surface["Y"][row, wedge],
                     color="black", linewidth=0.55)
    axes[1].plot([0, 19], [0, 0], color="black", linewidth=0.9)
    axes[1].plot([0, 0], [0, 19], color="black", linewidth=0.9)
    axes[1].set_title("One Rotational Period", fontsize=12, pad=4)
    for ax in axes:
        ax.set_xlim(-20, 20)
        ax.set_ylim(-20, 20)
    finish(fig, "periodic_boundaries_bw", pad=0.0)


def ssr_comparison_bw():
    """Full transport versus first-wall surface-source replay across all coils."""
    directory = Path("/Users/tkiker/.claude/jobs/876e30be/tmp")
    full = np.load(directory / "coil_real_master.npz")
    replay = np.load(directory / "coil_real_pristine.npz")
    cells = full["coil_cells"]
    full_shape = full["coil_fast_flux"] / full["coil_fast_flux"].max()
    replay_shape = replay["coil_fast_flux"] / replay["coil_fast_flux"].max()

    fig, ax = plt.subplots(figsize=(7.2, 3.05))
    x = np.arange(len(cells))
    width = 0.42
    ax.bar(x - width / 2, full_shape, width, facecolor="white",
           edgecolor="black", linewidth=0.85, label="Full Transport")
    ax.bar(x + width / 2, replay_shape, width, facecolor="white",
           edgecolor="black", linewidth=0.85, hatch="////", label="SSR Replay")
    ax.set_xticks(x)
    ax.set_xticklabels([int(cell) for cell in cells], fontsize=7)
    ax.set_xlabel("Coil Cell")
    ax.set_ylabel("Load / Peak Load")
    ax.legend(loc="upper left")
    finish(fig, "ssr_comparison_bw")


def shield_surface_geometry(n_theta=33, n_phi=65, shield_mid_offset=0.79):
    """Approximate the Step-1 shield mid-surface by poloidal-normal offset."""
    base = load_vmec_surface(n_theta=n_theta, n_phi=n_phi)
    radius = base["R"]
    height = base["Z"]
    dr = 0.5 * (np.roll(radius, -1, axis=0) - np.roll(radius, 1, axis=0))
    dz = 0.5 * (np.roll(height, -1, axis=0) - np.roll(height, 1, axis=0))
    nr, nz = dz.copy(), -dr.copy()
    mag = np.hypot(nr, nz)
    nr, nz = nr / mag, nz / mag
    outward = ((radius - radius.mean(axis=0, keepdims=True)) * nr +
               (height - height.mean(axis=0, keepdims=True)) * nz)
    sign = np.where(outward >= 0, 1.0, -1.0)
    radius = radius + shield_mid_offset * nr * sign
    height = height + shield_mid_offset * nz * sign
    ph = np.broadcast_to(base["phi"][None, :], radius.shape)
    return radius * np.cos(ph), radius * np.sin(ph), height, base["phi"]


def patch_geometry_bw():
    """Actual repeated QH regions for the reducer and backfire interventions."""
    x, y, z, phi = shield_surface_geometry()
    coils = parse_makegrid()
    ntheta, nphi = x.shape
    theta_bins = np.linspace(0, ntheta - 1, 5, dtype=int)
    toroidal_sector = np.floor(((phi % (np.pi / 2)) / (np.pi / 2)) * 4).astype(int)
    toroidal_sector = np.clip(toroidal_sector, 0, 3)

    fig = plt.figure(figsize=(8.5, 5.3))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(x, y, z, color="white", edgecolor="none", alpha=0.10,
                    shade=False)
    for coil in coils:
        ax.plot(coil[:, 0], coil[:, 1], coil[:, 2], color="0.55",
                linewidth=0.52, alpha=0.75)

    segments = []
    for theta_index in theta_bins:
        theta_index = min(theta_index, ntheta - 1)
        for j in range(nphi - 1):
            segments.append([(x[theta_index, j], y[theta_index, j], z[theta_index, j]),
                             (x[theta_index, j + 1], y[theta_index, j + 1],
                              z[theta_index, j + 1])])
    for j in range(nphi):
        if j == 0 or toroidal_sector[j] != toroidal_sector[j - 1]:
            for i in range(ntheta - 1):
                segments.append([(x[i, j], y[i, j], z[i, j]),
                                 (x[i + 1, j], y[i + 1, j], z[i + 1, j])])
    ax.add_collection3d(Line3DCollection(segments, colors="0.45", linewidths=0.38))

    def add_patch(toroidal_index, poloidal_index, shade):
        faces = []
        i0, i1 = theta_bins[poloidal_index], theta_bins[poloidal_index + 1]
        for i in range(i0, i1):
            for j in range(nphi - 1):
                if toroidal_sector[j] != toroidal_index:
                    continue
                faces.append([(x[i, j], y[i, j], z[i, j]),
                              (x[i, j + 1], y[i, j + 1], z[i, j + 1]),
                              (x[i + 1, j + 1], y[i + 1, j + 1], z[i + 1, j + 1]),
                              (x[i + 1, j], y[i + 1, j], z[i + 1, j])])
        ax.add_collection3d(Poly3DCollection(faces, facecolor=shade,
                                             edgecolor="black", linewidth=0.22,
                                             alpha=0.88))

    add_patch(3, 3, "0.18")
    add_patch(2, 3, "0.72")
    equalize_3d(ax, [(x.ravel(), y.ravel(), z.ravel())] +
                    [(coil[:, 0], coil[:, 1], coil[:, 2]) for coil in coils],
                z_scale=1.28, zoom=1.40)
    ax.view_init(elev=31, azim=-56)
    ax.set_axis_off()
    fig.subplots_adjust(0, 0, 1, 1)
    finish(fig, "patch_geometry_bw", pad=0.0, tight=False)


def response_matrix_bw():
    """Full 12-by-20 finite-intervention response matrix, including signs."""
    with open(RESPONSE_DATA / "patch_response_vectors.csv") as handle:
        rows = list(csv.DictReader(handle))
    coils = np.arange(11, 31)
    patches = [row["patch"] for row in rows]
    delta = np.array([[float(row[f"dR_c{coil}"]) for coil in coils] for row in rows])
    baseline_peak = 3.4521079534924136e-2
    magnitude = np.abs(delta) / baseline_peak

    fig, ax = plt.subplots(figsize=(8.3, 4.2))
    image = ax.imshow(magnitude, cmap="Greys", origin="upper", aspect="auto",
                      vmin=0.0, vmax=np.percentile(magnitude, 98),
                      interpolation="nearest")
    for row in range(delta.shape[0]):
        for column in range(delta.shape[1]):
            if delta[row, column] > 0:
                ax.add_patch(Rectangle((column - 0.5, row - 0.5), 1, 1,
                                       facecolor="none", edgecolor="black",
                                       linewidth=0.0, hatch="////"))
    hot_column = int(np.where(coils == 20)[0][0])
    ax.add_patch(Rectangle((hot_column - 0.5, -0.5), 1, len(patches),
                           facecolor="none", edgecolor="black", linewidth=1.2))
    ax.set_xticks(np.arange(len(coils)))
    ax.set_xticklabels(coils, fontsize=7.3)
    ax.set_yticks(np.arange(len(patches)))
    ax.set_yticklabels([rf"$({tag[1]},{tag[3]})$" for tag in patches], fontsize=7.3)
    ax.set_xlabel("Coil Cell")
    ax.set_ylabel(r"Shield Patch $(i_{\phi},j_{\vartheta})$")
    bar = fig.colorbar(image, ax=ax, fraction=0.030, pad=0.025)
    bar.set_label(r"$|\Delta R_j|/R_{\max,0}$")
    ax.legend(handles=[Patch(facecolor="white", edgecolor="black", hatch="////",
                             label="Load Increase")], loc="upper left",
              bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    finish(fig, "response_matrix_bw", pad=0.03)


def current_response_curves():
    """Show the spread of coil responses and the reducer/backfire extremes."""
    with open(RESPONSE_DATA / "baseline_responses.csv") as handle:
        baseline_rows = list(csv.DictReader(handle))
    baseline = np.array([float(row["fast_flux"]) for row in baseline_rows])
    coils = np.array([int(float(row["coil_cell"])) for row in baseline_rows])
    with open(RESPONSE_DATA / "patch_response_vectors.csv") as handle:
        rows = list(csv.DictReader(handle))
    peak = baseline.max()

    fig, ax = plt.subplots(figsize=(8.6, 3.65))
    for row in rows:
        response = baseline + np.array([float(row[f"dR_c{coil}"]) for coil in coils])
        ax.plot(coils, response / peak, color="0.74", linewidth=0.85, zorder=1)

    styles = {
        "t3p3": dict(color="black", linewidth=1.8, marker="o",
                     markersize=4.2, label=r"Reducer $(i_\phi,j_\vartheta)=(3,3)$"),
        "t2p3": dict(color="black", linewidth=1.8, linestyle="--", marker="s",
                     markerfacecolor="white", markersize=4.2,
                     label=r"Backfire $(i_\phi,j_\vartheta)=(2,3)$"),
    }
    for tag, style in styles.items():
        row = next(item for item in rows if item["patch"] == tag)
        response = baseline + np.array([float(row[f"dR_c{coil}"]) for coil in coils])
        ax.plot(coils, response / peak, zorder=3, **style)

    ax.axhline(1.0, color="0.35", linewidth=0.9, linestyle=":")
    ax.set_xlim(coils.min() - 0.5, coils.max() + 0.5)
    ax.set_xlabel("Coil Cell", fontsize=13)
    ax.set_ylabel("Response / Baseline Peak", fontsize=13)
    ax.tick_params(labelsize=11.5)
    ax.legend(loc="upper left", fontsize=10.5)
    finish(fig, "current_response_curves", pad=0.03)


def main():
    title_stellarator_bw()
    qh_plasma_bw()
    qh_reactor_bw()
    tokamak_vs_stellarator_bw()
    vmec_coordinates_bw()
    source_births_bw()
    birth_density()
    spf_geometry_bw()
    sampler_validation()
    transport_validation()
    rate_steering()
    conformal_map()
    predictor_scan()
    config_montage_bw()
    paper_screenshots_bw()
    literature_crops_bw()
    radial_build_bw()
    parastell_iteration_bw()
    periodic_boundaries_bw()
    patch_geometry_bw()
    response_matrix_bw()
    current_response_curves()
    print(f"Wrote deck figures and engineering views to {FIGS}")


if __name__ == "__main__":
    main()
