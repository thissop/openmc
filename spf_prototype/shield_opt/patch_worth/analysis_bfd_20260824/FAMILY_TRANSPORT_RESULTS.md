# Equal-mass shield-family transport

Common design mass: **2.163×10⁶ kg**. Peak-value headline gate: **READY**. Hotspot-label gate: **READY**.

| family | peak coil | peak fast flux | rel. err. | Δ vs uniform | z | verdict |
|---|---:|---:|---:|---:|---:|---|
| uniform | c23 | 2.5449e-02 | 3.5% | +0.0% ± 4.9% | +0.00 | reference |
| clearance-driven | c23 | 2.5509e-02 | 3.3% | +0.2% ± 4.8% | +0.05 | statistically unresolved |
| response-aware | c23 | 2.5597e-02 | 3.5% | +0.6% ± 4.9% | +0.12 | statistically unresolved |
| linear-minimax stress test | c23 | 2.5222e-01 | 1.0% | +891.1% ± 35.9% | +87.96 | worsens peak |

## Integrity gates

- A peak value is headline-ready only when its relative error is ≤ 10%.
- Improvement or degradation requires |z| ≥ 2.0 after independent-error propagation.
- A hotspot label is treated as resolved only when the top coil exceeds the runner-up by ≥ 2.0 combined σ.
- Coil identity is compared by stable physical labels; each source JSON also retains the actual DAGMC tally cell IDs when split solids shift geometry IDs.
- The comparison concerns transported peak coil fast flux at matched design mass; it does not rehabilitate the linear surrogate used to propose the families.
