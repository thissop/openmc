"""Free-streaming cross-check ON THE CONFORMAL WALL (per Ethan: the square torus is for
axisymmetric Schwartz validation only; the source must also be validated free-streaming on
the real conformal geometry). We compute the analytic point-patch free-streaming NWL
directionality A/iso and B/iso on the SAME conformal wall grid the OpenMC StellaratorSource
was ray-traced onto (quasr<ID>_conformalmap.npz), and compare per patch and by wall region.

Analytic per wall patch w: sum over source voxels s (weight sqrt(g)) of the emission kernel
about the LOCAL field Bhat_s times the geometric factor cos_inc / d^2, front-facing only:
  iso  ~ sum g_s ;  A ~ sum g_s (3/2) sin^2 th_B ;  B ~ sum g_s 2(1/4 + 3/4 cos^2 th_B)
with cos th_B = uhat . Bhat_s, uhat = (r_w - r_s)/d, g_s = w_s (n_w . uhat)/d^2.
Directionality A/iso, B/iso cancels the source normalization (as in the square-torus check).

Caveat: front-facing (n_w.uhat>0) visibility; no wall self-occlusion. Exact for convex walls;
on strongly non-convex (bean) walls the residual vs the ray-traced OpenMC map is itself the
point-patch/occlusion limit (the paper's thesis). Run on a GENTLE device for the clean check.
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stellarator_geometry as sg  # noqa: E402

CM = 100.0
GAP_FRAC = 0.30
DATA = HERE.parent / "data"


def bhat_cart(BR, Bphi, BZ, phi):
    bx = BR * np.cos(phi) - Bphi * np.sin(phi)
    by = BR * np.sin(phi) + Bphi * np.cos(phi)
    bz = BZ
    n = np.sqrt(bx * bx + by * by + bz * bz) + 1e-30
    return np.stack([bx / n, by / n, bz / n], axis=-1)


def run(ID):
    f = dict(np.load(DATA / f"quasr{ID}_vmec_fluxmap.npz"))
    R, Z, phi, sqrtg = f["R"], f["Z"], f["phi"], f["sqrtg"]
    BR, Bphi, BZ = f["BR"], f["Bphi"], f["BZ"]
    nrho, nt, nz = R.shape
    # source voxels
    Rc = R * CM; Zc = Z * CM
    xv = (Rc * np.cos(phi)).ravel(); yv = (Rc * np.sin(phi)).ravel(); zv = Zc.ravel()
    Sv = np.stack([xv, yv, zv], axis=1)                       # (Nv,3)
    wv = np.abs(sqrtg).ravel(); wv = wv / wv.sum()
    Bv = bhat_cart(BR, Bphi, BZ, phi).reshape(-1, 3)          # (Nv,3)
    # conformal wall grid (LCFS + 0.30a offset), matching conformal_wallmap
    Rl, Zl = R[-1] * CM, Z[-1] * CM
    nR, nZ = sg.poloidal_outward_normals(Rl, Zl)
    a = 0.5 * (Rl.max() - Rl.min())
    Rw, Zw = sg.offset_surface(Rl, Zl, nR, nZ, GAP_FRAC * a)
    phiw = np.linspace(0.0, 2 * np.pi, nz, endpoint=False)
    Pg = sg._verts(Rw, Zw, phiw)                              # (nt,nz,3)
    # outward wall normal = poloidal outward normal (guaranteed outward, unlike the tangent
    # cross-product which flips sign around the torus) lifted to 3D at each toroidal angle
    nRw, nZw = sg.poloidal_outward_normals(Rw, Zw)            # (nt,nz), outward in (R,Z)
    cph = np.cos(phiw)[None, :]; sph = np.sin(phiw)[None, :]
    Nw = np.stack([nRw * cph, nRw * sph, np.broadcast_to(nZw, nRw.shape)], axis=-1)
    Nw = Nw / (np.linalg.norm(Nw, axis=-1, keepdims=True) + 1e-30)
    Pw = Pg.reshape(-1, 3); Nwf = Nw.reshape(-1, 3)
    Rw_flat = np.hypot(Pw[:, 0], Pw[:, 1])                    # for inboard/outboard split
    # analytic point-patch free-streaming, chunked over wall patches
    A_iso = np.empty(Pw.shape[0]); B_iso = np.empty(Pw.shape[0])
    for i0 in range(0, Pw.shape[0], 128):
        W = Pw[i0:i0 + 128]; NW = Nwf[i0:i0 + 128]
        d = W[:, None, :] - Sv[None, :, :]                    # (k,Nv,3)  r_w - r_s
        d2 = (d * d).sum(-1)
        dist = np.sqrt(d2)
        u = d / dist[..., None]                               # emission dir
        cos_inc = np.einsum('kj,kvj->kv', NW, u)              # n_w . u
        cos_inc = np.clip(cos_inc, 0, None)                   # front-facing
        cosB = np.einsum('vj,kvj->kv', Bv, u)                 # u . Bhat_s
        g = wv[None, :] * cos_inc / d2                        # geometric weight
        iso = g.sum(1)
        A = (g * 1.5 * (1.0 - cosB ** 2)).sum(1)
        B = (g * (0.5 + 1.5 * cosB ** 2)).sum(1)
        A_iso[i0:i0 + 128] = A / (iso + 1e-30)
        B_iso[i0:i0 + 128] = B / (iso + 1e-30)
    A_iso = A_iso.reshape(nt, nz); B_iso = B_iso.reshape(nt, nz)
    # OpenMC ray-traced maps -> directionality
    c = np.load(DATA / f"quasr{ID}_conformalmap.npz")
    unpol, perp, par = c["unpol"], c["perp"], c["par"]
    mc_A = perp / (unpol + 1e-30); mc_B = par / (unpol + 1e-30)
    ok = unpol > 0.05 * unpol.max()                          # bins with statistics
    def corr(x, y):
        x = x[ok]; y = y[ok]
        return float(np.corrcoef(x, y)[0, 1])
    # region split (inboard = inner-R third, outboard = outer-R third)
    Rf = Rw_flat.reshape(nt, nz)
    inb = Rf < np.quantile(Rf, 1 / 3); outb = Rf > np.quantile(Rf, 2 / 3)
    def region(mA_an, mA_mc, mask):
        return mA_an[mask].mean(), mA_mc[mask].mean()
    print(f"=== conformal free-streaming cross-check: device {ID} (nt={nt},nz={nz}) ===")
    print(f"per-bin correlation analytic<->OpenMC:  A/iso {corr(A_iso, mc_A):.3f}   B/iso {corr(B_iso, mc_B):.3f}")
    print(f"global mean directionality  A/iso: analytic {A_iso[ok].mean():.3f}  OpenMC {mc_A[ok].mean():.3f}"
          f"   |  B/iso: analytic {B_iso[ok].mean():.3f}  OpenMC {mc_B[ok].mean():.3f}")
    for name, mask in [("inboard", inb), ("outboard", outb)]:
        aA, mA = region(A_iso, mc_A, mask); aB, mB = region(B_iso, mc_B, mask)
        print(f"  {name:9s} A/iso  an {aA:.3f}  mc {mA:.3f}  ({100*abs(aA-mA)/mA:.1f}%)"
              f"   B/iso  an {aB:.3f}  mc {mB:.3f}  ({100*abs(aB-mB)/mB:.1f}%)")
    np.savez(DATA / f"quasr{ID}_conformal_xcheck.npz",
             A_iso_an=A_iso, B_iso_an=B_iso, A_iso_mc=mc_A, B_iso_mc=mc_B, ok=ok, Rf=Rf)
    print(f"saved data/quasr{ID}_conformal_xcheck.npz")


if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 59509)
