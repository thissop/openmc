# Optimal SPF polarization for minimum first-wall peaking

Real OpenMC (+DAGMC) transport of the native `StellaratorSource` on each device's **VMEC** equilibrium (real √g births + real b̂), free-streaming through near-void to a **DAGMC conformal first wall** (the VMEC LCFS offset outward by 0.30·a along its poloidal normal); every neutron's wall crossing captured via `surf_source_write` and mapped (KD-tree on the wall grid) to (poloidal θ, toroidal φ). NWL is exactly linear in the emission quadrupole `a₂` (kernel `1+a₂P₂(cosθ_B)`), so `NWL(a₂)=NWL_unpol+a₂(NWL_par−NWL_unpol)` is built from the three transport maps; `a₂∈[−1,+1]` (−1 = pure perpendicular / A mode, +1 = pure parallel / B–C mode, 0 = unpolarized). Optimum = min peaking factor `PF = max(NWL)/⟨NWL⟩`.

| device | nfp | optimal a₂ | achievable (a,b,c) | PF unpol | PF optimal | reduction | MC linearity |
|---|---|---|---|---|---|---|---|
| 803097 | 3 | **-0.83** (perpendicular) | (0.84, 0.08, 0.08) | 1.40 | 1.38 | **+1.4%** | 0.054 |
| 886079 | 2 | **+0.83** (parallel) | (0.06, 0.89, 0.06) | 1.78 | 1.49 | **+16.1%** | 0.057 |
| 932746 | 3 | **+0.23** (parallel) | (0.26, 0.48, 0.26) | 1.38 | 1.32 | **+4.4%** | 0.053 |
| 59509 | 3 | **+0.22** (parallel) | (0.26, 0.48, 0.26) | 1.37 | 1.31 | **+4.4%** | 0.052 |
| 1960314 | 5 | **+0.60** (parallel) | (0.13, 0.73, 0.13) | 1.60 | 1.51 | **+6.0%** | 0.058 |

**MC linearity** = mean|NWL_perp+NWL_par−2·NWL_unpol| / ⟨NWL_unpol⟩ (→0 confirms the load is linear in a₂; residual is Monte-Carlo noise).

Grayscale (θ,φ) NWL/⟨NWL⟩ maps (unpolarized | optimal) per device: `figs/wallmaps/wallmap_<ID>.png`.
