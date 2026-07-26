"""Phase F scaffolding: Bayesian co-optimization of coil standoff + non-uniform neutron-shield
allocation, at FIXED outer envelope, subject to a TBR floor.

Why Bayesian and not gradient descent: per the PI, Monte-Carlo coil-dose gradients are too noisy
to trust cell-by-cell. So we do NOT descend the noisy MC surface. Instead:

  * BOOTSTRAP the search on the CHEAP, noise-free exponential free-streaming surrogate
    (attenuation_surrogate.coil_dose) -- this is where almost all objective evaluations happen;
  * only OCCASIONALLY spend an expensive MC call to verify/correct a promising candidate,
    via the `mc_evaluate` hook (a stub a cluster run fills in; defaults to the surrogate so the
    whole thing runs end-to-end offline);
  * fit a Gaussian process to whatever we have observed and pick the next design by Expected
    Improvement. The GP smooths over the little residual noise the way a raw gradient cannot.

Objective (MINIMIZE): peak coil dose (KS soft-max of the surrogate dose field) + a TBR-floor
penalty + small smoothness / material / standoff regularizers. Peak dose is turned into a
"minimum magnet lifetime" by a simple inverse (dose ~ 1/lifetime), so minimizing peak dose
== maximizing the minimum coil lifetime, which is the real engineering objective.

Design vector = [coil standoff scalar(s) in cm] + [shield-field Fourier coefficients c_k].
The shield field feeds ThicknessField.delta(coeffs, basis) -> the fixed-envelope breeder<->shield
trade. The standoff scalar moves the coils outward: it multiplies the baseline dose by an extra
exp(-k_env * standoff) attenuation (each cm of standoff is modeled as ~one cm of shield-grade
free-streaming attenuation; k_env defaults to 1/lambda_shield). Moving coils out always helps
dose, so a small standoff engineering cost (bigger machine) is charged to keep the trade honest.

Backends, in preference order: scikit-optimize (skopt.gp_minimize) -> BoTorch -> a self-contained
numpy/scipy GP + Expected-Improvement loop implemented here. It is NEVER random search: even the
fallback fits a real GP each iteration and maximizes EI. `optimize()` reports which backend ran.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

import attenuation_surrogate as att
from thickness_field import ThicknessField

# ----------------------------------------------------------------------------------------------
# Backend detection (honest: report what actually ran).
# ----------------------------------------------------------------------------------------------
try:  # pragma: no cover - environment dependent
    import skopt  # noqa: F401
    from skopt import gp_minimize as _skopt_gp_minimize
    from skopt.space import Real as _SkReal
    _HAVE_SKOPT = True
except Exception:  # pragma: no cover
    _HAVE_SKOPT = False

try:  # pragma: no cover - environment dependent
    import botorch  # noqa: F401
    _HAVE_BOTORCH = True
except Exception:  # pragma: no cover
    _HAVE_BOTORCH = False


# ==============================================================================================
# Problem configuration
# ==============================================================================================
@dataclass
class OptProblem:
    """Everything the objective needs beyond the design vector. Bundled so the surrogate and the
    (future) MC evaluator see an identical problem definition."""
    baseline_coil_dose_field: np.ndarray      # (nTor, nPol) baseline coil dose per sightline
    tf: ThicknessField                        # the fixed-envelope trade field
    basis: np.ndarray                         # (K, nTor, nPol) Fourier basis from tf.fourier_basis
    tbr0: float                               # baseline global TBR
    tbr_floor: float = 1.05                   # TBR must not drop below this
    lam_shield: float = att.LAMBDA_DEFAULT["shield"]
    lam_breeder: float = att.LAMBDA_DEFAULT["breeder"]

    # standoff model
    n_standoff: int = 1                       # number of coil-standoff scalars in the design vector
    k_env: Optional[float] = None             # 1/cm; default 1/lam_shield (shield-grade attenuation)
    standoff_max_cm: float = 15.0             # design bound for each standoff scalar
    coeff_bound: float = 5.0                  # design bound for each shield coefficient (+/-)

    # objective weights
    rho: float = 8.0                          # KS soft-max sharpness
    penalty: float = 1.0e3                    # TBR-floor penalty weight (dominates near the floor)
    lam_smooth: float = 5.0                   # manufacturability / MHD smoothness regularizer
    lam_material: float = 0.5                 # added-shield-material (cost) regularizer
    lam_standoff: float = 2.0                 # per-cm engineering cost of moving coils out

    # TBR linearization inputs (calibrate to OpenMC; scalar broadcasts over the grid)
    breeder_sensitivity: float = 3.0e-3       # dTBR per cm of local breeder thinning
    cell_area_frac: Optional[np.ndarray] = None  # weights summing to 1; default uniform

    # lifetime conversion: lifetime ~ life_ref / peak_dose
    life_ref: float = 1.0

    def __post_init__(self):
        self.baseline_coil_dose_field = np.asarray(self.baseline_coil_dose_field, float)
        if self.k_env is None:
            self.k_env = 1.0 / self.lam_shield
        if self.cell_area_frac is None:
            n = self.baseline_coil_dose_field.size
            self.cell_area_frac = np.full(self.baseline_coil_dose_field.shape, 1.0 / n)
        self.cell_area_frac = np.asarray(self.cell_area_frac, float)
        self.k = att.trade_k(self.lam_shield, self.lam_breeder)
        self.n_coeff = int(self.basis.shape[0])
        self.n_dim = self.n_standoff + self.n_coeff
        # index (within the coeff block) of the constant / DC basis mode -- the uniform-field knob.
        # fourier_basis does NOT put it first (it emits cos(nfp*phi) before the m=n=0 term), so we
        # detect it as the basis row with (near) zero angular variance.
        var = np.array([np.var(row) for row in self.basis])
        self.dc_index = int(np.argmin(var))

    # ---- design-vector plumbing -------------------------------------------------------------
    def split(self, design):
        d = np.asarray(design, float)
        standoff = d[: self.n_standoff]
        coeffs = d[self.n_standoff:]
        return standoff, coeffs

    def uniform_field_design(self, bias=-6.0, standoff=0.0):
        """A design whose shield field is spatially uniform: only the DC coefficient is set (to
        `bias`, which passes through the sigmoid), all ripple coeffs zero. bias<<0 -> ~zero shield
        (the baseline build); bias>0 -> uniform shield added everywhere."""
        x = np.zeros(self.n_dim)
        x[: self.n_standoff] = standoff
        x[self.n_standoff + self.dc_index] = bias
        return x

    def bounds(self):
        """List of (lo, hi) per design dimension. Standoff in cm; coeffs dimensionless (feed a
        sigmoid inside ThicknessField.delta, so +/-5 already saturates the field)."""
        b = [(0.0, self.standoff_max_cm)] * self.n_standoff
        b += [(-self.coeff_bound, self.coeff_bound)] * self.n_coeff
        return b

    # ---- physics pieces ---------------------------------------------------------------------
    def effective_baseline(self, standoff):
        """Baseline coil dose after moving coils out by `standoff` cm. Multiple standoff scalars
        are treated as independent toroidal-band standoffs tiled across the toroidal axis; a single
        scalar (the common case) attenuates the whole field uniformly."""
        std_total = float(np.sum(standoff))  # simple additive standoff (cm of shield-grade path)
        return self.baseline_coil_dose_field * np.exp(-self.k_env * std_total)

    def dose_field(self, design):
        standoff, coeffs = self.split(design)
        delta = self.tf.delta(coeffs, self.basis)
        base = self.effective_baseline(standoff)
        return att.coil_dose(base, delta, self.k), delta

    def tbr_of(self, design):
        _, coeffs = self.split(design)
        delta = self.tf.delta(coeffs, self.basis)
        sens = np.full(delta.shape, self.breeder_sensitivity)
        return self.tbr0 + att.tbr_change(delta, sens, self.cell_area_frac)


# ==============================================================================================
# Objective + helpers (public API)
# ==============================================================================================
def objective(design, baseline_coil_dose_field, tf, basis, tbr0, tbr_floor=1.05,
              lam_shield=att.LAMBDA_DEFAULT["shield"], lam_breeder=att.LAMBDA_DEFAULT["breeder"],
              *, problem: Optional[OptProblem] = None, **kw):
    """Scalar to MINIMIZE. Either pass a prebuilt `problem=OptProblem(...)` (fast, avoids rebuilding
    every call) or the loose arguments (a fresh OptProblem is constructed from them + **kw).

    J = KS_softmax_peak_dose
        + penalty * max(0, tbr_floor - tbr_new)          # TBR floor (one-sided)
        + lam_smooth   * smoothness(delta)               # manufacturability
        + lam_material * added_shield_volume(delta)      # material spent
        + lam_standoff * sum(standoff)                   # bigger-machine engineering cost
    """
    p = problem or OptProblem(baseline_coil_dose_field, tf, basis, tbr0, tbr_floor,
                              lam_shield, lam_breeder, **kw)
    standoff, coeffs = p.split(design)
    dose, delta = p.dose_field(design)

    J_peak, _ = att.robust_coil_objective(dose, rho=p.rho)

    tbr_new = p.tbr0 + att.tbr_change(delta, np.full(delta.shape, p.breeder_sensitivity),
                                      p.cell_area_frac)
    tbr_pen = p.penalty * max(0.0, p.tbr_floor - tbr_new)

    reg = (p.lam_smooth * p.tf.smoothness(delta)
           + p.lam_material * p.tf.added_shield_volume(delta)
           + p.lam_standoff * float(np.sum(standoff)))

    return float(J_peak + tbr_pen + reg)


def peak_dose(design, problem: OptProblem):
    """True (hard) peak of the coil-dose field for reporting -- the soft-max is only for the
    smooth objective."""
    dose, _ = problem.dose_field(design)
    return float(np.max(dose))


def magnet_lifetime(design, problem: OptProblem):
    """Minimum coil lifetime = life_ref / peak_dose. Maximizing this == minimizing peak dose."""
    return float(problem.life_ref / (peak_dose(design, problem) + 1e-300))


def breeder_sacrificed_fraction(design, problem: OptProblem):
    """Fraction of the baseline breeder thickness traded away (area-averaged delta / t_breeder0)."""
    _, coeffs = problem.split(design)
    delta = problem.tf.delta(coeffs, problem.basis)
    return float(np.mean(delta) / problem.tf.t_b0)


# ==============================================================================================
# Self-contained Gaussian process + Expected Improvement (the numpy/scipy fallback backend)
# ==============================================================================================
class _GP:
    """Compact GP regressor: zero-mean on standardized y, anisotropic-free (single lengthscale on
    [0,1]-scaled inputs) RBF kernel + white noise. Hyperparameters (lengthscale, signal var, noise)
    are fit by maximizing the log marginal likelihood with a short multi-start L-BFGS. Deliberately
    minimal but a REAL GP -- not a stand-in for random search."""

    def __init__(self, jitter=1e-8):
        self.jitter = jitter

    def _kernel(self, A, B, ls, sv):
        d2 = np.sum(A**2, 1)[:, None] + np.sum(B**2, 1)[None, :] - 2.0 * A @ B.T
        d2 = np.maximum(d2, 0.0)
        return sv * np.exp(-0.5 * d2 / (ls**2))

    def _nll(self, theta, X, y):
        ls, sv, nz = np.exp(theta)
        K = self._kernel(X, X, ls, sv) + (nz + self.jitter) * np.eye(len(X))
        try:
            L = np.linalg.cholesky(K)
        except np.linalg.LinAlgError:
            return 1e25
        alpha = np.linalg.solve(L.T, np.linalg.solve(L, y))
        nll = 0.5 * y @ alpha + np.sum(np.log(np.diag(L))) + 0.5 * len(X) * np.log(2 * np.pi)
        return float(nll)

    def fit(self, X, y):
        from scipy.optimize import minimize
        self.X = np.asarray(X, float)
        y = np.asarray(y, float)
        self.ymu, self.ysd = y.mean(), y.std() + 1e-12
        self.y = (y - self.ymu) / self.ysd

        best = None
        # multi-start over log(lengthscale, signal_var, noise)
        starts = [np.log([0.2, 1.0, 1e-3]), np.log([0.5, 1.0, 1e-2]),
                  np.log([1.0, 1.0, 1e-4]), np.log([0.1, 0.5, 1e-3])]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for s in starts:
                try:
                    r = minimize(self._nll, s, args=(self.X, self.y), method="L-BFGS-B",
                                 bounds=[(np.log(1e-2), np.log(5.0)),
                                         (np.log(1e-3), np.log(10.0)),
                                         (np.log(1e-6), np.log(1.0))])
                    if best is None or r.fun < best.fun:
                        best = r
                except Exception:
                    continue
        theta = best.x if best is not None else starts[0]
        self.ls, self.sv, self.nz = np.exp(theta)
        K = self._kernel(self.X, self.X, self.ls, self.sv) + (self.nz + self.jitter) * np.eye(len(self.X))
        self.L = np.linalg.cholesky(K)
        self.alpha = np.linalg.solve(self.L.T, np.linalg.solve(self.L, self.y))
        return self

    def predict(self, Xs):
        Xs = np.asarray(Xs, float)
        Ks = self._kernel(Xs, self.X, self.ls, self.sv)
        mu = Ks @ self.alpha
        v = np.linalg.solve(self.L, Ks.T)
        var = self.sv - np.sum(v**2, 0)
        var = np.maximum(var, 1e-12)
        mu = mu * self.ysd + self.ymu
        sd = np.sqrt(var) * self.ysd
        return mu, sd


def _expected_improvement(mu, sd, y_best, xi=0.01):
    """EI for MINIMIZATION. y_best is the incumbent (best-so-far) objective."""
    from scipy.stats import norm
    imp = y_best - mu - xi
    z = imp / sd
    ei = imp * norm.cdf(z) + sd * norm.pdf(z)
    return np.where(sd > 1e-12, ei, 0.0)


def _scale(X, lo, hi):
    return (np.asarray(X, float) - lo) / (hi - lo)


def _unscale(U, lo, hi):
    return lo + np.asarray(U, float) * (hi - lo)


# ==============================================================================================
# The optimizer
# ==============================================================================================
def optimize(problem: OptProblem, n_bootstrap=30, n_verify=10, n_init=8,
             mc_evaluate: Optional[Callable[[np.ndarray], float]] = None,
             x0=None, seed=0, n_candidates=2000, backend="auto", verbose=True):
    """Bayesian optimization of the shield/standoff design.

    Phases:
      * n_init      : initial space-filling designs, evaluated on the cheap surrogate objective.
      * n_bootstrap : GP+EI iterations, still on the cheap surrogate.
      * n_verify    : GP+EI iterations whose objective goes through `mc_evaluate` (the expensive
                      hook). `mc_evaluate` defaults to the surrogate objective so this runs offline.

    backend: "auto" (skopt -> botorch -> gp_ei), or force one of "skopt"/"gp_ei".

    Returns a dict: best_design, best_value, peak_dose, min_magnet_lifetime, tbr, backend, history.
    `history` is a list of dicts {iter, phase, design, value}.
    """
    rng = np.random.default_rng(seed)
    lo = np.array([b[0] for b in problem.bounds()], float)
    hi = np.array([b[1] for b in problem.bounds()], float)
    dim = problem.n_dim

    def surrogate_obj(design):
        return objective(design, None, None, None, None, problem=problem)

    if mc_evaluate is None:
        mc_evaluate = surrogate_obj  # offline default: MC "verify" == surrogate

    # -- backend selection ---------------------------------------------------------------------
    chosen = backend
    if backend == "auto":
        chosen = "skopt" if _HAVE_SKOPT else ("botorch" if _HAVE_BOTORCH else "gp_ei")
    if chosen == "botorch" and not _HAVE_BOTORCH:
        chosen = "gp_ei"
    if chosen == "skopt" and not _HAVE_SKOPT:
        chosen = "gp_ei"

    # skopt path: it manages its own GP+EI; we still stage bootstrap vs verify by swapping the
    # objective callable partway through via a mutable phase counter.
    if chosen == "skopt":  # pragma: no cover - skopt not installed in this env
        return _optimize_skopt(problem, n_bootstrap, n_verify, n_init, mc_evaluate,
                               surrogate_obj, lo, hi, seed, verbose)

    # -- self-contained GP + EI ----------------------------------------------------------------
    history = []
    X, Y = [], []

    def record(design, value, phase, it):
        X.append(np.asarray(design, float))
        Y.append(float(value))
        history.append(dict(iter=it, phase=phase, design=np.asarray(design, float), value=float(value)))

    it = 0
    # initial designs (Latin-hypercube-ish: stratified uniform per dim)
    init_pts = []
    if x0 is not None:
        init_pts.append(np.clip(np.asarray(x0, float), lo, hi))
    u = (rng.permutation(np.tile(np.arange(n_init), (dim, 1)).T).astype(float)
         + rng.random((n_init, dim))) / n_init
    for row in u:
        init_pts.append(_unscale(row, lo, hi))
    for d in init_pts:
        record(d, surrogate_obj(d), "init", it); it += 1

    total = n_bootstrap + n_verify
    for step in range(total):
        phase = "bootstrap" if step < n_bootstrap else "verify"
        Us = _scale(np.array(X), lo, hi)
        gp = _GP().fit(Us, np.array(Y))
        y_best = min(Y)
        best_idx = int(np.argmin(Y))
        incU = _scale(X[best_idx], lo, hi)
        # Candidate set for the inner EI maximization (standard practice to give EI a rich pool in
        # >10-D, where a single uniform draw rarely proposes a coherent multi-lever move):
        #   (a) global uniform exploration,
        #   (b) multi-scale Gaussian perturbations of the incumbent,
        #   (c) coordinate sweeps -- vary one dim across its full range at the incumbent. These
        #       expose the dominant 1-D levers (coil standoff, DC shield bias) that a full-dim
        #       random draw drowns out.
        pool = [rng.random((n_candidates, dim))]
        for scale in (0.02, 0.08, 0.25):
            pool.append(np.clip(incU + scale * rng.standard_normal((n_candidates // 4, dim)), 0, 1))
        grid = np.linspace(0.0, 1.0, 13)
        for j in range(dim):
            sweep = np.tile(incU, (grid.size, 1))
            sweep[:, j] = grid
            pool.append(sweep)
        cand = np.vstack(pool)
        mu, sd = gp.predict(cand)
        ei = _expected_improvement(mu, sd, y_best)
        nxtU = cand[int(np.argmax(ei))]
        nxt = _unscale(nxtU, lo, hi)
        val = surrogate_obj(nxt) if phase == "bootstrap" else mc_evaluate(nxt)
        record(nxt, val, phase, it); it += 1
        if verbose and (step % 10 == 0 or step == total - 1):
            print(f"  [{phase:9s}] iter {it:3d}  best J = {min(Y):.4f}")

    best_i = int(np.argmin(Y))
    best = X[best_i]
    return dict(
        best_design=best,
        best_value=float(Y[best_i]),
        peak_dose=peak_dose(best, problem),
        min_magnet_lifetime=magnet_lifetime(best, problem),
        tbr=problem.tbr_of(best),
        backend="gp_ei",
        history=history,
    )


def _optimize_skopt(problem, n_bootstrap, n_verify, n_init, mc_evaluate, surrogate_obj,
                    lo, hi, seed, verbose):  # pragma: no cover - skopt not installed here
    space = [_SkReal(float(l), float(h)) for l, h in zip(lo, hi)]
    phase = {"n": 0}
    history = []

    def f(x):
        ph = "bootstrap" if phase["n"] < (n_init + n_bootstrap) else "verify"
        v = surrogate_obj(x) if ph != "verify" else mc_evaluate(x)
        history.append(dict(iter=phase["n"], phase=ph, design=np.asarray(x, float), value=float(v)))
        phase["n"] += 1
        return v

    res = _skopt_gp_minimize(f, space, n_calls=n_init + n_bootstrap + n_verify,
                             n_initial_points=n_init, random_state=seed, acq_func="EI")
    best = np.asarray(res.x, float)
    return dict(best_design=best, best_value=float(res.fun),
                peak_dose=peak_dose(best, problem),
                min_magnet_lifetime=magnet_lifetime(best, problem),
                tbr=problem.tbr_of(best), backend="skopt", history=history)


# ==============================================================================================
# Plotting
# ==============================================================================================
def plot_convergence(history, path="figs/bayes_convergence.png"):
    """Running-min of the objective vs evaluation, colored by phase. matplotlib Agg (headless)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import os

    vals = np.array([h["value"] for h in history], float)
    running = np.minimum.accumulate(vals)
    phases = [h["phase"] for h in history]

    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(np.arange(len(vals)), running, "-", color="0.2", lw=1.8, label="Best So Far")
    colors = {"init": "#888888", "bootstrap": "#1f77b4", "verify": "#d62728"}
    for ph in ("init", "bootstrap", "verify"):
        idx = [i for i, p in enumerate(phases) if p == ph]
        if idx:
            ax.scatter(idx, vals[idx], s=18, color=colors[ph], label=ph.capitalize(), zorder=3)
    ax.set_xlabel("Evaluation")
    ax.set_ylabel("Objective (Peak Coil Dose + Penalties)")
    ax.set_title("Bayesian Shield/Standoff Co-Optimization Convergence")
    ax.legend(frameon=False)
    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True) if os.path.dirname(path) else None
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


# ==============================================================================================
# Synthetic problem builder (used by the demo and the tests)
# ==============================================================================================
def make_synthetic_problem(nfp=4, n_coils=20, hot=1.5, ripple=0.15, tbr0=1.15,
                           n_tor=48, n_pol=36, M=2, N=0):
    """An nfp-periodic baseline coil-dose field with a HOT INBOARD BAND (poloidal angle ~pi) and a
    20-coil toroidal ripple. Returns (problem, basis, tf).

    Basis note: the dominant, shieldable feature is the POLOIDAL inboard band; the 20-coil toroidal
    ripple is too high-frequency for a smooth shield field to chase and the nfp modulation is mild.
    So the demo defaults to poloidal-only shield modes (N=0, M=2 -> 5 coeffs), which keeps the BO in
    a tractable ~6-D space AND is the physically right control. Raise N to add toroidal shield
    freedom (more expressive, but harder for BO with a fixed evaluation budget)."""
    tor = np.linspace(0.0, 360.0, n_tor, endpoint=False)
    pol = np.linspace(0.0, 360.0, n_pol, endpoint=False)
    tf = ThicknessField(nfp=nfp, toroidal_angles_deg=tor, poloidal_angles_deg=pol,
                        t_breeder0=80.0, t_shield0=20.0, t_breeder_min=15.0)
    basis = tf.fourier_basis(M=M, N=N)

    TH, PH = tf.POL, tf.TOR  # (nTor, nPol) radians
    # inboard band: peak at theta=pi (inboard midplane), von-Mises-like in poloidal angle
    inboard = np.exp(hot * (np.cos(TH - np.pi) - 1.0))          # in (0,1], peaks at theta=pi
    coil_ripple = 1.0 + ripple * np.cos(n_coils * PH)           # 20-coil toroidal ripple
    field_period = 1.0 + 0.10 * np.cos(nfp * PH)                # nfp field-period modulation
    base = 100.0 * (0.3 + inboard) * coil_ripple * field_period
    base = np.maximum(base, 1.0)

    problem = OptProblem(baseline_coil_dose_field=base, tf=tf, basis=basis, tbr0=tbr0,
                         tbr_floor=1.05, n_standoff=1)
    return problem, basis, tf


# ==============================================================================================
# Offline demo
# ==============================================================================================
def _demo():
    print("=" * 78)
    print("Bayesian shield/standoff co-optimization -- offline surrogate demo")
    print("=" * 78)
    print(f"backends: skopt={_HAVE_SKOPT}  botorch={_HAVE_BOTORCH}  "
          f"-> using {'skopt' if _HAVE_SKOPT else ('botorch' if _HAVE_BOTORCH else 'gp_ei fallback')}")

    problem, basis, tf = make_synthetic_problem()
    dim = problem.n_dim
    print(f"design dim = {dim}  ({problem.n_standoff} standoff + {problem.n_coeff} shield coeffs)")

    # starting point: zero standoff, strongly-negative DC bias -> ~zero shield field.
    x_start = problem.uniform_field_design(bias=-6.0, standoff=0.0)
    start_peak = peak_dose(x_start, problem)
    start_tbr = problem.tbr_of(x_start)
    start_life = magnet_lifetime(x_start, problem)
    print("\n-- starting design (baseline: no standoff, ~no shield trade) --")
    print(f"   peak coil dose = {start_peak:9.3f}")
    print(f"   TBR            = {start_tbr:9.4f}   (floor {problem.tbr_floor})")
    print(f"   min lifetime   = {start_life:9.4g}  (arb. units, ~1/peak dose)")

    res = optimize(problem, n_bootstrap=30, n_verify=10, n_init=8, x0=x_start, seed=1,
                   verbose=True)

    opt = res["best_design"]
    opt_peak = res["peak_dose"]
    opt_tbr = res["tbr"]
    opt_life = res["min_magnet_lifetime"]
    frac = breeder_sacrificed_fraction(opt, problem)
    standoff, _ = problem.split(opt)

    print("\n-- optimized design --")
    print(f"   backend used   = {res['backend']}")
    print(f"   coil standoff  = {float(np.sum(standoff)):9.3f} cm")
    print(f"   peak coil dose = {opt_peak:9.3f}   ({100*(1-opt_peak/start_peak):+.1f}% vs start)")
    print(f"   TBR            = {opt_tbr:9.4f}   (floor {problem.tbr_floor}; "
          f"{'OK' if opt_tbr >= problem.tbr_floor - 1e-6 else 'VIOLATED'})")
    print(f"   min lifetime   = {opt_life:9.4g}   (x{opt_life/start_life:.2f} gain)")
    print(f"   breeder sacrificed fraction = {frac:6.3f}")

    path = plot_convergence(res["history"], path="figs/bayes_convergence.png")
    print(f"\nconvergence plot saved -> {path}")

    dose_drop = 100 * (1 - opt_peak / start_peak)
    print("\nSUMMARY: peak dose {} while TBR {} floor.".format(
        f"reduced {dose_drop:.1f}%" if dose_drop > 0 else "NOT reduced",
        "held above" if opt_tbr >= problem.tbr_floor - 1e-6 else "VIOLATED"))
    return res


if __name__ == "__main__":
    _demo()
