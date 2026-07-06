# Blanket & reactor parameters (paper-grounded, cited)

For Tiers 5–6 (scattering + FLiBe blanket/TBR). Every value below is traced to a
paper. The configuration is **ARC/LIBRA-class** (your group's designs), with a
**FLiBe** breeder — the closest published analog is Peterson's **LIBRA**.

## Primary references
- **Sorbom et al. 2015, "ARC"** — Fusion Eng. Des. 100, 378; arXiv:1409.3540.
- **Segantin et al. 2020, "Optimization of TBR in ARC"** — Fusion Eng. Des. 154, 111531; PSFC/JA-20-56.
- **Peterson et al. 2022, "LIBRA"** (FLiBe liquid-immersion blanket) — Fusion Sci. Technol., doi:10.1080/15361055.2022.2078136.
- **Dunn, Woller, Peterson 2025, "ClLiF"** — Fusion Sci. Technol., doi:10.1080/15361055.2025.2504844.
- FLiBe properties: ACS J. Chem. Eng. Data (2023); Romatoski salt review (MIT).
- TBR target/benchmarks: EU DEMO HCPB≈1.18 / WCLL≈1.15 (goal 1.15); FNG HCPB OpenMC benchmark (doi:10.1080/15361055.2025.2567167).

## Cited parameter table
| parameter | value | source |
|---|---|---|
| Major radius R₀ | 3.3 m | Sorbom 2015 |
| Minor radius a | 1.1 m | Sorbom 2015 |
| Toroidal field B₀ | 9.2 T | Sorbom 2015 |
| Fusion power | ~500 MW | Sorbom 2015 |
| ARC layer stack (cavity→out) | 0.1 cm W / 1 cm Inconel 718 / 2 cm FLiBe channel / 1 cm Be, VV immersed in bulk FLiBe tank | Sorbom 2015 |
| FLiBe operating temperature | 900 K (inlet 800 K) | Sorbom 2015 |
| ARC TBR | ≥1.1 | Sorbom 2015 |
| ARC TBR (baseline, recomputed) | ~1.07; 1.08±0.004 with 1 cm Be | Segantin 2020 |
| **Optimal Li-6 enrichment** | **20–50%** (not the previously assumed 90%) | Segantin 2020 |
| Structural materials ranked | Eurofer97 ≈ V-15Cr-5Ti > Inconel 718 (for TBR) | Segantin 2020 |
| Be multiplier reaction | ⁹Be(n,2n), threshold ~1.9 MeV | Segantin 2020 |
| FLiBe composition | 2LiF·BeF₂ = Li₂BeF₄; atoms Li:Be:F = 2:1:4 (28.6/14.3/57.1 at%) | ACS JCED 2023; LIBRA |
| FLiBe density (operating) | ρ(600 °C)=1.98, ρ(700 °C)=1.96 g/cm³ → use **1.94 g/cm³** | ACS JCED 2023 |
| FLiBe melting point | 459 °C | ACS JCED 2023 |
| Natural Li-6 abundance | ~7.5 at% | ACS JCED 2023 |
| **LIBRA blanket thickness** | **~40 cm** FLiBe (min radial), 1000 kg | Peterson 2022 |
| **LIBRA TBR** | **~1.4 at 50% Li-6**; gains plateau beyond ~30% | Peterson 2022 |
| LIBRA concept | liquid-immersion FLiBe tank surrounding the VV | Peterson 2022 |

## Chosen model configuration (scoping; square cross-section retained)
Geometry: keep the **validated square-cross-section box-torus** (so the
free-streaming "sandwich" anchor maps 1:1 to Tier-2b), scaled to **ARC-class
size** (minor radius ≈ 1 m). Layers lining the cavity wall, outward:

| layer | material | thickness | basis |
|---|---|---|---|
| first wall | **tungsten** (ρ=19.3) | 0.1 cm | Sorbom (your W choice) |
| structure | reduced-activation steel (Fe-9Cr, Eurofer-like; ρ≈7.8) | 1 cm | Segantin (Eurofer > Inconel); avoids Nb/Mo data |
| multiplier | beryllium (ρ=1.85) | 1 cm | Sorbom/Segantin |
| breeder | **FLiBe** Li₂BeF₄ (ρ=1.94), Li-6 enriched | **40 cm** | LIBRA |

- **Li-6 enrichment:** baseline **30%** (LIBRA plateau / Segantin 20–50%), with a
  **scan** {7.5% (natural), 30, 50, 90} — the scan doubles as verification (TBR
  rises then plateaus, matching LIBRA's reported behavior).
- **TBR we expect:** order **1.1–1.4** for a ~40 cm FLiBe blanket — i.e., a
  number we can sanity-check against LIBRA/ARC rather than an analytic truth.

## Caveats (state in results)
- 294 K cross sections used with 900 K FLiBe density (Doppler/thermal-scatter
  effects small at 14 MeV; note the inconsistency).
- Single-zone scoping blanket; no coolant channels, no depletion, no tritium
  extraction; Inconel approximated by Fe-9Cr steel.
- Square cross-section (aspect 2.0) kept for the anchor, not ARC's D-shape / aspect 3.
