#!/usr/bin/env python
"""Re-solve the constrained minimax candidate at the COMMON design-family budget M*=design_families.Mstar,
so all four transported families share one mass (valid equal-mass Pareto). Same LP as minimax_candidate.py
(min peak-coil t s.t. matched mass + 0<=delta<=35 over the 12 built sectors); only the mass target changes.
Writes minimax_candidate_matched.npz. Still a linear-surrogate SCREENING guess -> transport-verify."""
import os, json, numpy as np
from scipy.optimize import linprog
HERE = os.path.dirname(os.path.abspath(__file__)); SH = "/Users/tkiker/lima_share"
RHO_S, RHO_B, DMAX = 15.6, 1.94, 35.0
patch_order = ["t0p2","t0p3","t1p0","t1p1","t1p2","t1p3","t2p0","t2p1","t2p2","t2p3","t3p2","t3p3"]
frac = np.genfromtxt(f"{HERE}/b_fd_fractional.csv", delimiter=",", skip_header=1)
Rbase = np.ravel(np.load(f"{HERE}/../cluster_pull/coil_step1b_uniform.npz")["coil_fast_flux"]).astype(float)
E = frac * Rbase[None, :]
mc = json.load(open(f"{SH}/masses_local.json")); z = mc["zerodelta"]
mass_coef = np.array([(mc[p]["V_shield"] - z["V_shield"]) / 20.0 for p in patch_order]) * (RHO_S - RHO_B) / 1000.0
Mstar = float(np.load(f"{HERE}/design_families.npz")["Mstar"])           # COMMON budget (16-sector uniform 10cm)
Mmax = float(np.sum(mass_coef * DMAX))
print(f"common M*={Mstar:.0f} kg; max reachable over 12 sectors={Mmax:.0f} kg -> {'FEASIBLE' if Mstar<=Mmax else 'INFEASIBLE'}")
assert Mstar <= Mmax, "12 built sectors capped at 35cm cannot reach the common budget"
n = 12
c = np.zeros(n + 1); c[-1] = 1.0
A_ub = np.zeros((20, n + 1)); b_ub = np.zeros(20)
for j in range(20):
    A_ub[j, :n] = E[:, j] / 20.0; A_ub[j, -1] = -1.0; b_ub[j] = -Rbase[j]
A_eq = np.zeros((1, n + 1)); A_eq[0, :n] = mass_coef; b_eq = [Mstar]
r = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(0, DMAX)] * n + [(None, None)], method="highs")
assert r.success, r.message
delta = r.x[:n]; Rpred = Rbase + (E.T @ (delta / 20.0))
jb, jm = int(Rbase.argmax()), int(Rpred.argmax())
print(f"baseline peak c{11+jb} {Rbase.max():.4e}")
print(f"minimax(matched) predicted peak c{11+jm} {Rpred.max():.4e} ({100*(Rpred.max()-Rbase.max())/Rbase.max():+.1f}% vs base); "
      f"limiting-coil c{11+jb}->c{11+jm}")
print("delta cm: " + ", ".join(f"{p}:{d:.1f}" for p, d in zip(patch_order, delta)))
np.savez(f"{HERE}/minimax_candidate_matched.npz", delta=delta, Rpred=Rpred, Rbase=Rbase, Mstar=Mstar,
         patch_order=np.array(patch_order))
print("WROTE minimax_candidate_matched.npz  (mass-matched to the common family budget)")
