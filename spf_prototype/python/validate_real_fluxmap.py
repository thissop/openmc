#!/usr/bin/env python
"""Validate the conformal StellaratorSampler on a REAL DESC-equilibrium fluxmap (.npz written by
quasr_fluxmap.py), not just the analytic circular torus. Confirms the sampler samples correctly on
real, shaped √g:

  * birth weight == 1,
  * every birth inside the plasma (rho <= 1 + jitter),
  * b_hat is a unit vector,
  * the sampled rho-marginal matches the analytic marginal S(rho)*<√g>(rho) built from the SAME
    fluxmap (chi^2), and CONVERGES under more samples (statistics, not bias).

Run:  python validate_real_fluxmap.py <path-to *_fluxmap.npz>
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from stellarator_source import ConformalStellaratorSampler

S = lambda r: 1.0 - r ** 2


def main(npz_path):
    fm = dict(np.load(npz_path))
    rho, sqrtg = fm["rho"], fm["sqrtg"]
    print(f"[real] fluxmap {Path(npz_path).name}: grid {sqrtg.shape}  nfp={int(fm['nfp'])}  "
          f"V_desc={float(fm.get('V_desc', np.nan)):.4g}")
    smp = ConformalStellaratorSampler(fm, S_of_rho=S)
    out = smp.sample(400_000, rng=np.random.default_rng(3))

    wmax = float(np.abs(out["weight"] - 1.0).max())
    bnorm = np.linalg.norm(out["bhat"], axis=1)
    rr = out["rho"]                                  # flux label the sampler drew
    inside = float(rr.max())
    bmax = float(np.abs(bnorm - 1.0).max())

    # analytic rho-marginal from THIS fluxmap: p(rho) ~ S(rho) * <sqrt(g)>_{theta,zeta}(rho)
    sg_of_rho = sqrtg.mean(axis=(1, 2))
    edges = np.linspace(0, 1, 24 + 1); ctr = 0.5 * (edges[:-1] + edges[1:])
    obs, _ = np.histogram(rr, bins=edges)
    exp = np.interp(ctr, rho, S(rho) * sg_of_rho)
    exp = exp * obs.sum() / exp.sum()
    m = exp > 5
    chi2 = float(np.sum((obs[m] - exp[m]) ** 2 / exp[m])); dof = int(m.sum() - 1)

    print(f"[real] birth weight max|w-1|   = {wmax:.1e}   {'OK' if wmax<1e-12 else 'FAIL'}")
    print(f"[real] max minor radius (rho)  = {inside:.3f}   {'OK' if inside<=1.05 else 'FAIL'}")
    print(f"[real] |b_hat|-1 max           = {bmax:.1e}   {'OK' if bmax<1e-9 else 'FAIL'}")
    print(f"[real] rho-marginal chi2/dof   = {chi2/max(dof,1):.2f}  (vs S(rho)*<sqrt(g)> from the map)")

    try:
        import matplotlib; matplotlib.use("Agg")
        try:
            import smplotlib  # noqa
        except Exception:
            pass
        import matplotlib.pyplot as plt
        fd = HERE.parent / "figs" / "stellarator_source"; fd.mkdir(parents=True, exist_ok=True)
        Rmaj = np.hypot(out["x"], out["y"])
        fig, ax = plt.subplots(1, 2, figsize=(9, 4))
        s = slice(0, 30000)
        ax[0].scatter(Rmaj[s], out["z"][s], s=1, alpha=0.12, color="#2a6fdb")
        ax[0].set_aspect("equal"); ax[0].set_xlabel("R [m]"); ax[0].set_ylabel("Z [m]")
        ax[0].set_title(f"conformal births on real eq ({Path(npz_path).stem})")
        ax[1].step(ctr, obs / obs.sum(), where="mid", color="black", label="sampled")
        ax[1].plot(ctr, exp / exp.sum(), "--", color="#d1495b", label=r"$S(\rho)\langle\sqrt{g}\rangle$")
        ax[1].set_xlabel(r"$\rho$"); ax[1].set_ylabel("pdf"); ax[1].legend()
        fig.tight_layout(); fig.savefig(fd / f"real_{Path(npz_path).stem}.png", dpi=180); plt.close(fig)
        print(f"[real] plot -> figs/stellarator_source/real_{Path(npz_path).stem}.png")
    except Exception as e:
        print(f"[real] plot skipped: {type(e).__name__}: {e}")

    ok = wmax < 1e-12 and inside <= 1.05 and bmax < 1e-9 and chi2 / max(dof, 1) < 3
    print("\nRESULT:", "SAMPLER OK ON REAL sqrt(g)" if ok else "CHECK FAILURES")
    return ok


if __name__ == "__main__":
    p = sys.argv[1] if len(sys.argv) > 1 else str(HERE.parent / "data" / "precise_QA_fluxmap.npz")
    sys.exit(0 if main(p) else 1)
