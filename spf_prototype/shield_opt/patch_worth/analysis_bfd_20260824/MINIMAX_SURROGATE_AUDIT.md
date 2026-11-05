# Minimax linear-surrogate audit

Common family mass: **2.163×10⁶ kg**. Each response column was measured for one isolated +20 cm patch.

| cap | enforce predicted R≥0 | feasible | predicted peak change | min predicted R | negative coils |
|---:|:---:|:---:|---:|---:|---:|
| 20 cm | no | yes | +1.2% | -1.548e-01 | 13 |
| 25 cm | no | yes | -101.2% | -1.445e-01 | 20 |
| 30 cm | no | yes | -100.1% | -1.404e-01 | 20 |
| 35 cm | no | yes | -100.0% | -1.379e-01 | 20 |
| 20 cm | yes | **no** | — | — | — |
| 35 cm | yes | **no** | — | — | — |

## Verdict

- The mass target is reachable with every patch capped at 20 cm, so mass feasibility is not the failure.
- Even at the measured 20 cm amplitude, the unconstrained linear combination predicts negative responses; enforcing non-negativity makes the equal-mass problem infeasible.
- Therefore the matched `famM` transport is an extrapolation stress test, not a trustworthy optimizer ranking.
- A production minimax design needs sequential re-linearization/direct-forward updates, a lower intervention budget, or an explicitly nonlinear surrogate trained on combined designs.
