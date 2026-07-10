"""Does SPF global steering benefit correlate with a geometry-only STANDOFF metric
and NOT with a pure-plasma SHAPE metric? (docs/notes/RIPPLE_AND_PREDICTOR.md sec 4-5)

Y  = SPF global benefit = (PF_unpol - min_a2 PF(a2)) / PF_unpol  [from transport maps]
Standoff metrics (geometry only, no transport): isotropic-flux proxy peaking, and the
  normalized standoff at the peak sightline.
Shape metrics (pure plasma, no wall): LCFS non-axisymmetric shaping amplitude, nematic Q.
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stellarator_geometry as sg  # noqa: E402

CM = 100.0
GAP_FRAC = 0.30
DEV = [1052272, 1090019, 11357, 1642553, 180790, 1854549, 1960314, 24285, 262171, 59509, 66633, 803097, 883496, 886079, 932746]
DATA = HERE.parent / "data"


def spearman(x, y):
    xr = np.argsort(np.argsort(x)).astype(float)
    yr = np.argsort(np.argsort(y)).astype(float)
    xr -= xr.mean(); yr -= yr.mean()
    return float((xr @ yr) / (np.linalg.norm(xr) * np.linalg.norm(yr) + 1e-30))


def benefit(ID):
    """SPF global peaking-reduction benefit from the conformal-wall NWL maps."""
    f = np.load(DATA / f"quasr{ID}_conformalmap.npz")
    unpol, par, dA = f["unpol"], f["par"], f["dA"]
    d = dA / dA.sum()
    def PF(m):
        return m.max() / (m * d).sum()
    a2 = np.linspace(-1, 1, 201)
    pfs = np.array([PF(unpol + t * (par - unpol)) for t in a2])
    pf0 = PF(unpol)
    imin = pfs.argmin()
    return dict(pf_unpol=pf0, pf_opt=pfs[imin], a2_opt=a2[imin],
                Y=(pf0 - pfs[imin]) / pf0)


def geometry(ID):
    """Wall points + source voxels from the VMEC fluxmap; standoff & shape metrics."""
    f = dict(np.load(DATA / f"quasr{ID}_vmec_fluxmap.npz"))
    R, Z, phi, sqrtg = f["R"], f["Z"], f["phi"], f["sqrtg"]        # (nrho,nt,nz)
    nrho, nt, nz = R.shape
    # --- source voxels (cm), weight = |sqrt(g)| (volume element ~ reactivity support) ---
    Rc, Zc, Pc = R * CM, Z * CM, phi
    xv = (Rc * np.cos(Pc)).ravel(); yv = (Rc * np.sin(Pc)).ravel(); zv = Zc.ravel()
    wv = np.abs(sqrtg).ravel(); wv = wv / wv.sum()
    Sv = np.stack([xv, yv, zv], axis=1)
    # --- conformal wall: LCFS offset outward by 0.30*a (exactly as conformal_wallmap) ---
    Rl, Zl = R[-1] * CM, Z[-1] * CM
    nR, nZ = sg.poloidal_outward_normals(Rl, Zl)
    a = 0.5 * (Rl.max() - Rl.min())
    Rw, Zw = sg.offset_surface(Rl, Zl, nR, nZ, GAP_FRAC * a)
    phiw = np.linspace(0.0, 2 * np.pi, nz, endpoint=False)
    Pw = sg._verts(Rw, Zw, phiw).reshape(-1, 3)                     # (nt*nz,3)
    # wall cell areas
    Pgrid = Pw.reshape(nt, nz, 3)
    dth, dph = 2 * np.pi / nt, 2 * np.pi / nz
    Pt = (np.roll(Pgrid, -1, 0) - np.roll(Pgrid, 1, 0)) / (2 * dth)
    Pp = (np.roll(Pgrid, -1, 1) - np.roll(Pgrid, 1, 1)) / (2 * dph)
    dA = (np.linalg.norm(np.cross(Pt, Pp), axis=-1) * dth * dph).ravel()
    dAn = dA / dA.sum()
    # --- geometry-only isotropic flux proxy G(w) = sum_v w_v / |w-r_v|^2 (no transport) ---
    # chunk over wall points to bound memory
    G = np.empty(Pw.shape[0]); dmin = np.empty(Pw.shape[0])
    for i0 in range(0, Pw.shape[0], 256):
        W = Pw[i0:i0 + 256]                                        # (k,3)
        d2 = ((W[:, None, :] - Sv[None, :, :]) ** 2).sum(-1)       # (k,nv)
        G[i0:i0 + 256] = (wv[None, :] / d2).sum(1)
        dmin[i0:i0 + 256] = np.sqrt(d2.min(1))
    # standoff metrics (geometry only)
    PF_geo = G.max() / (G * dAn).sum()
    standoff_peak_ratio = dmin[G.argmax()] / (dmin * dAn).sum()    # low => near-field spike
    spike_area = dAn[G > 0.9 * G.max()].sum()                      # area fraction of the hot spot
    # shape metrics (pure plasma, no wall)
    Rax = Rl.mean(1, keepdims=True); Zax = Zl.mean(1, keepdims=True)  # axisymmetric ref
    shaping_amp = np.sqrt(((Rl - Rax) ** 2 + (Zl - Zax) ** 2).mean()) / a
    # nematic Q of LCFS poloidal-normal orientation (director alignment)
    psi = np.arctan2(nZ, nR)
    Q = np.abs(np.exp(2j * psi).mean())
    return dict(PF_geo=PF_geo, standoff_peak_ratio=standoff_peak_ratio,
                spike_area=spike_area, shaping_amp=shaping_amp, nematic_Q=Q)


if __name__ == "__main__":
    rows = []
    for ID in DEV:
        b, g = benefit(ID), geometry(ID)
        rows.append({"ID": ID, **b, **g})
    keys = ["ID", "pf_unpol", "a2_opt", "Y", "PF_geo", "standoff_peak_ratio",
            "spike_area", "shaping_amp", "nematic_Q"]
    hdr = "  ".join(f"{k:>10}" for k in keys)
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print("  ".join(f"{r[k]:>10.4f}" if k != "ID" else f"{r[k]:>10d}" for k in keys))
    Y = np.array([r["Y"] for r in rows])
    print("\nSpearman rho of predictors vs SPF benefit Y (n=15):")
    print(f"  STANDOFF  PF_geo               : {spearman([r['PF_geo'] for r in rows], Y):+.3f}   (expect - : peakier geom => less steerable)")
    print(f"  STANDOFF  standoff_peak_ratio  : {spearman([r['standoff_peak_ratio'] for r in rows], Y):+.3f}   (expect + : more standoff at peak => more steerable)")
    print(f"  STANDOFF  spike_area           : {spearman([r['spike_area'] for r in rows], Y):+.3f}   (expect + : broader hot region => more steerable)")
    print(f"  SHAPE     shaping_amp          : {spearman([r['shaping_amp'] for r in rows], Y):+.3f}   (expect ~0: plasma shape decoupled)")
    print(f"  SHAPE     nematic_Q            : {spearman([r['nematic_Q'] for r in rows], Y):+.3f}   (expect ~0: plasma shape decoupled)")
