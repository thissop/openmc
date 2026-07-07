#!/usr/bin/env python
"""Validate the conformal StellaratorSampler against an ANALYTIC circular-torus sqrt(g), where
every marginal is known in closed form. Proves the algorithm independent of DESC:

  * birth weight is exactly 1,
  * the rho marginal matches  p(rho) ~ S(rho)*rho  AND its chi2 CONVERGES under grid refinement
    (so the residual is the deterministic discretization bias, not an algorithm error),
  * the theta marginal is OUTBOARD-biased  p(theta) ~ <R>(theta)  (the toroidal differential-volume
    effect Peterson's paper is about) -- a uniform-theta sampler would FAIL this,
  * every birth lies inside the plasma.

Run:  python test_stellarator_source.py     (pure numpy; no DESC/OpenMC needed)
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from stellarator_source import ConformalStellaratorSampler, circular_torus_fluxmap

R0, A, N = 6.0, 2.0, 400_000
S = lambda r: 1.0 - r ** 2                                        # parabolic flux function
trap = getattr(np, "trapezoid", getattr(np, "trapz", None))


def chi2(obs, exp):
    exp = exp * obs.sum() / exp.sum()
    m = exp > 5
    return float(np.sum((obs[m] - exp[m]) ** 2 / exp[m])), int(m.sum() - 1)


def run(nr, plot=False):
    fm, _ = circular_torus_fluxmap(R0=R0, a=A, nr=nr, nt=96, nz=64)
    smp = ConformalStellaratorSampler(fm, S_of_rho=S)
    out = smp.sample(N, rng=np.random.default_rng(7))
    Rmaj = np.hypot(out["x"], out["y"]); rr = np.hypot(Rmaj - R0, out["z"]) / A
    th = np.arctan2(out["z"], Rmaj - R0)

    wmax = float(np.abs(out["weight"] - 1.0).max())
    # rho marginal
    edges = np.linspace(0, 1, 24 + 1); ctr = 0.5 * (edges[:-1] + edges[1:])
    obs_r, _ = np.histogram(rr, bins=edges)
    c2r, dofr = chi2(obs_r.astype(float), S(ctr) * ctr)
    # theta marginal vs analytic outboard model  ~ integral S(rho) rho (R0 + a rho cos th) drho
    te = np.linspace(-np.pi, np.pi, 36 + 1); tc = 0.5 * (te[:-1] + te[1:])
    obs_t, _ = np.histogram(th, bins=te)
    rq = np.linspace(1e-3, 1, 400)
    Aint, Bint = trap(S(rq) * rq, rq), trap(S(rq) * rq ** 2, rq)
    exp_t = R0 * Aint + A * Bint * np.cos(tc)
    c2t, doft = chi2(obs_t.astype(float), exp_t)
    c2t_uni, _ = chi2(obs_t.astype(float), np.ones_like(tc))
    stats = dict(wmax=wmax, rmax=float(rr.max()), c2r=c2r / max(dofr, 1),
                 c2t=c2t / max(doft, 1), c2t_uni=c2t_uni / max(doft, 1),
                 frac_out=float((np.abs(th) < np.pi / 2).mean()))
    if plot:
        _plot(Rmaj, out["z"], ctr, obs_r, S(ctr) * ctr, tc, obs_t, exp_t)
    return stats


def _plot(Rmaj, Zc, ctr, obs_r, exp_r, tc, obs_t, exp_t):
    try:
        import matplotlib; matplotlib.use("Agg")
        try:
            import smplotlib  # noqa
        except Exception:
            pass
        import matplotlib.pyplot as plt
        fd = HERE.parent / "figs" / "stellarator_source"; fd.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(1, 3, figsize=(13, 4))
        s = slice(0, 30000)
        ax[0].scatter(Rmaj[s], Zc[s], s=1, alpha=0.15, color="#2a6fdb")
        ax[0].set_aspect("equal"); ax[0].set_xlabel("R [m]"); ax[0].set_ylabel("Z [m]")
        ax[0].set_title("conformal births (poloidal)")
        ax[1].step(ctr, obs_r / obs_r.sum(), where="mid", color="black", label="sampled")
        ax[1].plot(ctr, exp_r / exp_r.sum(), "--", color="#d1495b", label=r"$S(\rho)\,\rho$")
        ax[1].set_xlabel(r"$\rho$"); ax[1].set_ylabel("pdf"); ax[1].legend()
        ax[2].step(tc, obs_t / obs_t.sum(), where="mid", color="black", label="sampled")
        ax[2].plot(tc, exp_t / exp_t.sum(), "--", color="#d1495b", label="outboard model")
        ax[2].axhline(1.0 / tc.size, ls=":", color="grey", label="uniform (wrong)")
        ax[2].set_xlabel(r"$\theta$ (0 = outboard)"); ax[2].set_ylabel("pdf"); ax[2].legend()
        fig.tight_layout(); fig.savefig(fd / "circular_torus_validation.png", dpi=180); plt.close(fig)
        print("[V] plot -> figs/stellarator_source/circular_torus_validation.png")
    except Exception as e:
        print(f"[V] plot skipped: {type(e).__name__}: {e}")


def main():
    coarse = run(20)
    fine = run(48, plot=True)
    print(f"[V] birth weight max|w-1|          = {fine['wmax']:.1e}    {'OK' if fine['wmax']<1e-12 else 'FAIL'}")
    print(f"[V] max normalized minor radius    = {fine['rmax']:.3f}    {'OK' if fine['rmax']<=1.02 else 'FAIL'}")
    print(f"[V] rho marginal chi2/dof  nr=20 -> nr=48 : {coarse['c2r']:.2f} -> {fine['c2r']:.2f}   "
          f"({'CONVERGES' if fine['c2r']<coarse['c2r'] else 'NOT converging!'})")
    print(f"[V] theta marginal chi2/dof (outboard model)= {fine['c2t']:.2f}   "
          f"vs uniform-theta = {fine['c2t_uni']:.1f}  (uniform rejected)")
    print(f"[V] fraction born on outboard half          = {fine['frac_out']:.3f}  (>0.5 expected)")
    ok = (fine['wmax'] < 1e-12 and fine['rmax'] <= 1.02 and fine['c2r'] < coarse['c2r']
          and fine['c2r'] < 3 and fine['c2t'] < 3 and fine['c2t_uni'] > 10 and fine['frac_out'] > 0.5)
    print("\nRESULT:", "ALL SAMPLER GATES PASSED" if ok else "SOME GATES FAILED -- investigate")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
