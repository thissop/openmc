"""Uncollided NWL 'culprit map' for QUASR device 59509.

Attribute the peak first-wall neutron loading back to the plasma source regions that
CAUSE it, for unpolarized vs spin-polarized (SPF) emission. Pure analytic point-kernel
free-streaming (NO Monte Carlo, NO ParaStell). Reuses the source-voxel loader, conformal
wall construction (LCFS + 0.30a outward) and angular factors of
python/conformal_freestream_xcheck.py.

Forward uncollided NWL on wall patch k:
    NWL_k = sum_v S_v g(theta_vk) max(0, n_k.u_vk) / |w_k - r_v|^2
with u_vk = (w_k - r_v)/|w_k - r_v|, theta_vk = angle(u_vk, Bhat_v),
S_v = |sqrt(g)|_v (volume element ~ reactivity support, as in the existing pipeline).
Angular factor g: unpol g=1; A-mode g=(3/2) sin^2 (perp to B); B/C g=(1/4 + 3/4 cos^2)
(parallel to B).

Hot region H = wall patches carrying the top 5% of the wall AREA, ranked by NWL_k
(per mode). Culprit field:
    C_v = sum_{k in H} S_v g(theta_vk) max(0, n_k.u_vk)/|w_k - r_v|^2 dA_k
normalized to sum 1 (fraction of the peak-region load each voxel drives).
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SPF = HERE.parent
sys.path.insert(0, str(SPF / "python"))
sys.path.insert(0, str(SPF / "sweep"))
import stellarator_geometry as sg  # noqa: E402

CM = 100.0
GAP_FRAC = 0.30
DATA = SPF / "data"
OUTFIG = HERE / "figs"
OUTDATA = HERE / "data" / "processed"
OUTFIG.mkdir(parents=True, exist_ok=True)
OUTDATA.mkdir(parents=True, exist_ok=True)

# angular factors g(cos_theta), cos_theta = u_hat . Bhat_v
GFUN = {
    "unpol": lambda c2: np.ones_like(c2),
    "A":     lambda c2: 1.5 * (1.0 - c2),          # (3/2) sin^2 theta  (perp to B)
    "BC":    lambda c2: 0.25 + 0.75 * c2,          # 1/4 + 3/4 cos^2 theta (par to B)
}


def bhat_cart(BR, Bphi, BZ, phi):
    bx = BR * np.cos(phi) - Bphi * np.sin(phi)
    by = BR * np.sin(phi) + Bphi * np.cos(phi)
    bz = BZ
    n = np.sqrt(bx * bx + by * by + bz * bz) + 1e-30
    return np.stack([bx / n, by / n, bz / n], axis=-1)


def build(ID):
    f = dict(np.load(DATA / f"quasr{ID}_vmec_fluxmap.npz"))
    R, Z, phi, sqrtg = f["R"], f["Z"], f["phi"], f["sqrtg"]
    BR, Bphi, BZ = f["BR"], f["Bphi"], f["BZ"]
    nrho, nt, nz = R.shape
    # --- source voxels (cm), weight S_v = |sqrt(g)| ---
    Rc, Zc = R * CM, Z * CM
    xv = (Rc * np.cos(phi)).ravel(); yv = (Rc * np.sin(phi)).ravel(); zv = Zc.ravel()
    Sv = np.stack([xv, yv, zv], axis=1)                      # (Nv,3)
    wv = np.abs(sqrtg).ravel()                               # S_v (unnormalized)
    Bv = bhat_cart(BR, Bphi, BZ, phi).reshape(-1, 3)         # (Nv,3)
    # --- conformal wall: LCFS offset outward by 0.30 a (matches conformal_freestream) ---
    Rl, Zl = R[-1] * CM, Z[-1] * CM
    nRl, nZl = sg.poloidal_outward_normals(Rl, Zl)
    a = 0.5 * (Rl.max() - Rl.min())
    Rw, Zw = sg.offset_surface(Rl, Zl, nRl, nZl, GAP_FRAC * a)
    phiw = np.linspace(0.0, 2 * np.pi, nz, endpoint=False)
    Pg = sg._verts(Rw, Zw, phiw)                             # (nt,nz,3)
    nRw, nZw = sg.poloidal_outward_normals(Rw, Zw)
    cph = np.cos(phiw)[None, :]; sph = np.sin(phiw)[None, :]
    Nw = np.stack([nRw * cph, nRw * sph, np.broadcast_to(nZw, nRw.shape)], axis=-1)
    Nw = Nw / (np.linalg.norm(Nw, axis=-1, keepdims=True) + 1e-30)
    Pw = Pg.reshape(-1, 3); Nwf = Nw.reshape(-1, 3)
    # wall cell areas (metric on the conformal grid), same recipe as predictor_test
    dth, dph = 2 * np.pi / nt, 2 * np.pi / nz
    Pt = (np.roll(Pg, -1, 0) - np.roll(Pg, 1, 0)) / (2 * dth)
    Pp = (np.roll(Pg, -1, 1) - np.roll(Pg, 1, 1)) / (2 * dph)
    dA = (np.linalg.norm(np.cross(Pt, Pp), axis=-1) * dth * dph).ravel()  # (Nw,)
    return dict(Sv=Sv, wv=wv, Bv=Bv, Pw=Pw, Nw=Nwf, dA=dA,
                nrho=nrho, nt=nt, nz=nz,
                Rlcfs=Rl, Zlcfs=Zl, phiw=phiw,
                theta=f["theta"], zeta=f["zeta"])


def forward_nwl(G, modes, chunk=128):
    """NWL_k per mode: sum_v S_v g cos_inc / d2, front-facing only. Returns dict mode->(Nw,)."""
    Sv, wv, Bv, Pw, Nw = G["Sv"], G["wv"], G["Bv"], G["Pw"], G["Nw"]
    out = {m: np.empty(Pw.shape[0]) for m in modes}
    for i0 in range(0, Pw.shape[0], chunk):
        W = Pw[i0:i0 + chunk]; NW = Nw[i0:i0 + chunk]
        d = W[:, None, :] - Sv[None, :, :]
        d2 = (d * d).sum(-1)
        dist = np.sqrt(d2)
        u = d / dist[..., None]
        cos_inc = np.clip(np.einsum('kj,kvj->kv', NW, u), 0, None)
        c2 = np.einsum('vj,kvj->kv', Bv, u) ** 2
        geo = wv[None, :] * cos_inc / d2                     # S_v cos_inc / d2
        for m in modes:
            out[m][i0:i0 + chunk] = (geo * GFUN[m](c2)).sum(1)
    return out


def hot_mask(nwl, dA, area_frac=0.05):
    """Top `area_frac` of wall AREA ranked by NWL. Returns boolean mask over patches."""
    order = np.argsort(nwl)[::-1]
    cum = np.cumsum(dA[order]) / dA.sum()
    ncut = np.searchsorted(cum, area_frac) + 1
    H = np.zeros(nwl.shape[0], dtype=bool)
    H[order[:ncut]] = True
    return H


def attribute(G, mode, H, chunk=64):
    """C_v = sum_{k in H} S_v g cos_inc / d2 * dA_k, over hot patches only. (Nv,) normalized."""
    Sv, wv, Bv, Pw, Nw, dA = G["Sv"], G["wv"], G["Bv"], G["Pw"], G["Nw"], G["dA"]
    idx = np.nonzero(H)[0]
    C = np.zeros(Sv.shape[0])
    for i0 in range(0, idx.size, chunk):
        kk = idx[i0:i0 + chunk]
        W = Pw[kk]; NW = Nw[kk]; dAk = dA[kk]
        d = W[:, None, :] - Sv[None, :, :]
        d2 = (d * d).sum(-1)
        dist = np.sqrt(d2)
        u = d / dist[..., None]
        cos_inc = np.clip(np.einsum('kj,kvj->kv', NW, u), 0, None)
        c2 = np.einsum('vj,kvj->kv', Bv, u) ** 2
        geo = wv[None, :] * cos_inc / d2
        C += (geo * GFUN[mode](c2) * dAk[:, None]).sum(0)
    return C / (C.sum() + 1e-300)


def participation(C, wv, level):
    """Fraction of plasma VOLUME (weight=S_v=|sqrt g|) supplying `level` of the peak-region
    load, i.e. take voxels in descending C until cumulative C reaches `level`."""
    order = np.argsort(C)[::-1]
    cumC = np.cumsum(C[order])
    ncut = np.searchsorted(cumC, level) + 1
    vol = wv / wv.sum()
    return float(vol[order[:ncut]].sum())


def main(ID=59509):
    G = build(ID)
    nrho, nt, nz = G["nrho"], G["nt"], G["nz"]
    modes = ["unpol", "A", "BC"]
    nwl = forward_nwl(G, modes)
    dA = G["dA"]
    awm = {m: (nwl[m] * dA).sum() / dA.sum() for m in modes}   # area-weighted mean
    print(f"=== culprit map: device {ID}  (Nv={G['Sv'].shape[0]}, Nw={dA.size}) ===")
    res = {}
    for m in modes:
        H = hot_mask(nwl[m], dA, 0.05)
        peak = nwl[m].max() / awm[m]
        C = attribute(G, m, H)
        Cmap = C.reshape(nrho, nt, nz).sum(0)                 # project onto (theta,zeta)
        p50 = participation(C, G["wv"], 0.50)
        p90 = participation(C, G["wv"], 0.90)
        res[m] = dict(C=C, Cmap=Cmap, H=H, peak=peak, p50=p50, p90=p90)
        harea = dA[H].sum() / dA.sum()
        print(f"  {m:6s}  peak NWL/mean = {peak:5.2f}   hot-area = {100*harea:4.1f}%   "
              f"participation: 50% <- {100*p50:4.1f}% vol   90% <- {100*p90:4.1f}% vol")
    # concentration: is SPF shifting the culprits? correlation & centroid shift on (th,zeta)
    cu = res["unpol"]["Cmap"]; cA = res["A"]["Cmap"]
    corr = float(np.corrcoef(cu.ravel(), cA.ravel())[0, 1])
    print(f"  culprit-map corr(unpol, A) = {corr:.3f}   "
          f"(1=identical culprits, <1 => SPF redistributes)")
    # save
    np.savez(OUTDATA / f"culprit_{ID}.npz",
             C_unpol=res["unpol"]["C"], C_A=res["A"]["C"], C_BC=res["BC"]["C"],
             Cmap_unpol=res["unpol"]["Cmap"], Cmap_A=res["A"]["Cmap"], Cmap_BC=res["BC"]["Cmap"],
             H_unpol=res["unpol"]["H"], H_A=res["A"]["H"], H_BC=res["BC"]["H"],
             nwl_unpol=nwl["unpol"], nwl_A=nwl["A"], nwl_BC=nwl["BC"],
             Sv=G["Sv"], wv=G["wv"], dA=dA,
             theta=G["theta"], zeta=G["zeta"], nrho=nrho, nt=nt, nz=nz,
             p50=np.array([res[m]["p50"] for m in modes]),
             p90=np.array([res[m]["p90"] for m in modes]),
             peak=np.array([res[m]["peak"] for m in modes]), modes=np.array(modes))
    print(f"  saved {OUTDATA / f'culprit_{ID}.npz'}")
    plot(ID, G, res, modes)
    return res


def plot(ID, G, res, modes):
    import smplotlib  # noqa: F401  (house style)
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "savefig.dpi": 300, "savefig.bbox": "tight",
        "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "legend.fontsize": 9, "legend.frameon": False,
    })
    th = np.degrees(G["theta"]); ze = np.degrees(G["zeta"])
    ext = [ze.min(), ze.max(), th.min(), th.max()]
    cu, cA = res["unpol"]["Cmap"], res["A"]["Cmap"]
    vmax = max(cu.max(), cA.max())
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.6), constrained_layout=True)
    for a, m, ttl in [(ax[0], "unpol", "Unpolarized"), (ax[1], "A", "SPF (A, Perp)")]:
        im = a.imshow(res[m]["Cmap"], origin="lower", extent=ext, aspect="auto",
                      cmap="gray_r", vmin=0, vmax=vmax)
        a.set_xlabel("Toroidal Angle Zeta [Deg]")
        a.set_ylabel("Poloidal Angle Theta [Deg]")
        a.text(0.03, 0.94, ttl, transform=a.transAxes, fontsize=11,
               va="top", ha="left")
    cb = fig.colorbar(im, ax=ax[:2], location="bottom", shrink=0.9, pad=0.02)
    cb.set_label("Contribution To Peak NWL")
    # signed difference panel (RdBu_r, symmetric)
    diff = cA - cu
    dm = np.abs(diff).max()
    im2 = ax[2].imshow(diff, origin="lower", extent=ext, aspect="auto",
                       cmap="RdBu_r", vmin=-dm, vmax=dm)
    ax[2].set_xlabel("Toroidal Angle Zeta [Deg]")
    ax[2].set_ylabel("Poloidal Angle Theta [Deg]")
    ax[2].text(0.03, 0.94, "SPF Minus Unpol", transform=ax[2].transAxes,
               fontsize=11, va="top", ha="left")
    cb2 = fig.colorbar(im2, ax=ax[2], location="bottom", shrink=0.9, pad=0.02)
    cb2.set_label("Difference")
    for ext_ in ("pdf", "png"):
        fig.savefig(OUTFIG / f"culprit_{ID}.{ext_}")
    print(f"  saved {OUTFIG / f'culprit_{ID}.png'}")
    plt.close(fig)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 59509)
