# Radial-build validation vs. literature (2026-07-14)

Validation of the QH baseline radial build against published reactor neutronics. [F] = from
fetched primary PDF; [I] = engineering inference.

## Current build (107 cm plasma->coil) — per-layer verdict
| Layer | Current | Verdict | Fix |
|---|---|---|---|
| First wall | 4 cm SOLID W (19.3) | **WRONG** — ~2 mfp of W down-scatters/reflects source neutrons before the breeder, depressing TBR | 2 mm W armor + ~3 cm RAFM steel, homogenized ~7-8 g/cc effective |
| Breeder | 50 cm FLiBe, 90% Li6, NO multiplier | TBR<1 is **REALISTIC** (0.94 credible) | add **Be multiplier 1-2 cm** (+0.03 TBR/cm); re-optimize Li6 to ~30-50% (saturates, 90% not optimal) |
| Back wall | 3 cm RAFM | ok (ARIES-CS 5 cm) | 3-5 cm |
| Shield | 40 cm WC (15.6) | **VALIDATED** — ARIES-CS used 34 cm WC; density correct | keep WC; upgrade path W2B5 (10x better) or WC+TiH2 composite |
| Vac vessel | 10 cm steel, NO gap | **WRONG** — real VV 20-82 cm; MUST have cryostat gap + thermal shield before cold coil | 20-28 cm steel+borated water + 2 cm gap + 2-4 cm thermal shield |
| Total | 107 cm | **TOO THIN** — below ARIES-CS floor 130.7 cm | ~130 cm (QH has 162 cm coil-plasma available, so it fits easily) |

## Recommended literature-anchored build (inboard/min; thicker breeding elsewhere)
| Layer | Material | Thickness | Density | Citation |
|---|---|---|---|---|
| First wall | 2 mm W + ~3 cm RAFM, FLiBe-cooled | ~3.2 cm | ~7-8 g/cc | Kuang 2018 (ARC 3mm W/Inconel); Energies 14:1541 (DEMO-WCLL 2mm W+25mm EUROFER); UWFDM-1336 (ARIES-CS 3.8 cm) |
| **Be multiplier (ADD)** | Beryllium | 1-2 cm | 1.85 | Segantin 2020 (FED, PSFC/JA-20-56): +0.03 TBR/cm, lifts 1.01->1.07 |
| Breeder | FLiBe, Li6 ~30-90% | 50 cm | 1.94 | Sorbom 2015 (ARC); Segantin 2020 |
| Back wall | RAFM (EUROFER/F82H) | 3-5 cm | 7.9 | UWFDM-1336 |
| Shield | WC (or W2B5) | 34-50 cm | 15.6 (WC) | UWFDM-1336 (ARIES-CS 34 cm WC); Windsor 2021 NF 61:086018 (W2B5 10x) |
| Gap | cryostat/assembly | 2 cm | void | UWFDM-1336 |
| Vac vessel | double-wall SS316 + borated water | 20-28 cm | ~7.9 | Ioki 1998 (ITER 45-82 cm); UWFDM-1336 (28 cm) |
| Thermal shield/coil case | cooled steel + insulator | 2-4 cm | -- | EU-DEMO VVTS; UWFDM-1336 |
| **Total to winding pack** | | **~115-135 cm** | | matches ARIES-CS Delta_min 130.7 cm |

## HTS magnet radiation limits (these SET the shield thickness)
| Quantity | Value | Source |
|---|---|---|
| Fast fluence (E>0.1 MeV) REBCO/Nb3Sn | 3e18 (ARC) - 1e19 n/cm2 (ARIES-CS) | Sorbom 2015; El-Guebaly FST 54 (2008) 747 |
| dpa in Cu stabilizer | ~6e-3 (ARIES-CS) | El-Guebaly; UWFDM-872 |
| Insulator dose | 1e7 Gy (epoxy/cyanate) - 1e9 Gy (polyimide) | ITER; ARIES-CS |
| Peak winding-pack heating | ~1-2 mW/cm3 | ITER; ARIES-CS |

CAVEAT [F]: recent in-situ cryogenic 14-MeV REBCO irradiation shows Ic-degradation onset
ORDERS below the legacy 3e18 n/cm2 threshold (cold degrades ~1.6x faster than warm) —
Fischer SuST 31 (2018) 044006; MIT SuST 38 (2025) 015012. Shields sized to the classic
number may be optimistic for HTS lifetime.

## Bottom line
Shield material/thickness is fine (WC, ARIES-CS-validated). TBR=0.94 is physically real, not
a bug — the fix is a Be multiplier + Li6 re-optimization, NOT more FLiBe. The FW (solid W) and
VV (thin, no gap) are wrong and must be corrected. Total standoff should be ~130 cm, which the
QH coil set (162 cm min coil-plasma) accommodates with room to spare — our 107 cm left space
unused. Re-baseline quantitative results on the corrected build; the closed-loop demo itself is
robust to the baseline (relative comparison on a fixed grid).
