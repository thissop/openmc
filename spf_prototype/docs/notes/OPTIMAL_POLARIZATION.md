# Optimal SPF polarization for minimum first-wall peaking

Real OpenMC transport of the native `StellaratorSource` on each device's **VMEC** equilibrium (real √g births + real b̂), free-streaming to a torus first wall; every neutron's wall crossing captured via `surface_source_write` and histogrammed in (poloidal θ, toroidal φ). NWL is exactly linear in the emission quadrupole `a₂` (kernel `1+a₂P₂(cosθ_B)`), so `NWL(a₂)=NWL_unpol+a₂(NWL_par−NWL_unpol)` is built from the three transport maps; `a₂∈[−1,+1]` (−1 = pure perpendicular / A mode, +1 = pure parallel / B–C mode, 0 = unpolarized). Optimum = min peaking factor `PF = max(NWL)/⟨NWL⟩`.

| device | nfp | optimal a₂ | achievable (a,b,c) | PF unpol | PF optimal | reduction | MC linearity |
|---|---|---|---|---|---|---|---|
| 803097 | 3 | **-1.00** (perpendicular) | (1.00, 0.00, 0.00) | 1.59 | 1.48 | **+6.8%** | 0.041 |
| 886079 | 2 | **+0.59** (parallel) | (0.14, 0.72, 0.14) | 1.59 | 1.49 | **+6.1%** | 0.047 |
| 932746 | 3 | **-0.27** (perpendicular) | (0.46, 0.27, 0.27) | 1.25 | 1.23 | **+1.6%** | 0.044 |
| 59509 | 3 | **-0.50** (perpendicular) | (0.60, 0.20, 0.20) | 1.44 | 1.35 | **+6.2%** | 0.042 |
| 1960314 | 5 | **-1.00** (perpendicular) | (1.00, 0.00, 0.00) | 1.69 | 1.44 | **+14.6%** | 0.049 |

**MC linearity** = mean|NWL_perp+NWL_par−2·NWL_unpol| / ⟨NWL_unpol⟩ (→0 confirms the load is linear in a₂; residual is Monte-Carlo noise).

Grayscale (θ,φ) NWL/⟨NWL⟩ maps (unpolarized | optimal) per device: `figs/wallmaps/wallmap_<ID>.png`.
