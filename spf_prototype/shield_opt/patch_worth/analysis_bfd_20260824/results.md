# B_FD reconstruction — results (2026-08-24)

**Baseline:** `coil_step1b_uniform.npz` (shield 42.55 cm). **Response:** cell-integrated coil_fast_flux (NOT density). **Patches:** 12 single +20 cm sectors.

## Caveats (charter §0.5, §Statistical)
- patch shield-outside=40cm vs baseline 42.55cm -> ~2.55cm uniform confound in absolute dR
- patches and baseline both used no weight windows with matched settings/seed; there is no WW mismatch
- single disjoint-sector patches -> NOT an additive basis; combined designs need a new run

## Significance
- 98% of dR entries exceed 2.0sigma (independent-error bound).

## Peak migration & control (per patch)

baseline peak: coil 20 = 3.4521e-02, TBR 1.073

| patch | (i,j) | added·cm | dR_peak | new_peak | peakX | ->cell | dTBR% | class |
|---|---|---|---|---|---|---|---|---|
| t3p3 | (3,3) | 500 | -2.80e-02 | 6.722e-03 | 5.1x | 25 | +2.4 | reduce |
| t0p3 | (0,3) | 400 | -2.78e-02 | 6.843e-03 | 5.0x | 25 | -9.2 | reduce |
| t0p2 | (0,2) | 400 | -2.80e-02 | 6.898e-03 | 5.0x | 30 | -9.1 | reduce |
| t3p2 | (3,2) | 500 | -2.80e-02 | 6.937e-03 | 5.0x | 30 | -2.2 | reduce |
| t1p3 | (1,3) | 300 | -2.50e-02 | 9.506e-03 | 3.6x | 20 | +2.5 | reduce |
| t2p0 | (2,0) | 300 | -2.43e-02 | 1.032e-02 | 3.3x | 25 | +8.5 | tradeoff |
| t1p0 | (1,0) | 300 | -2.47e-02 | 3.875e-02 | 0.9x | 30 | -1.9 | tradeoff |
| t1p1 | (1,1) | 240 | -2.16e-02 | 4.153e-02 | 0.8x | 15 | -11.4 | tradeoff |
| t1p2 | (1,2) | 300 | -4.98e-03 | 4.158e-02 | 0.8x | 30 | -5.3 | tradeoff |
| t2p1 | (2,1) | 240 | -2.23e-02 | 4.896e-02 | 0.7x | 30 | -4.0 | tradeoff |
| t2p3 | (2,3) | 300 | +1.42e-02 | 4.973e-02 | 0.7x | 25 | -17.1 | backfire |
| t2p2 | (2,2) | 300 | -1.72e-02 | 4.999e-02 | 0.7x | 15 | -15.5 | tradeoff |

## SVD of uncertainty-whitened E (12x20)
- singular values: [739.4, 581.6, 245.5, 119.3, 78.5, 37.6, 34.4, 13.4, 8.8, 3.1, 1.8, 0.5]
- effective rank (entropy): 4.66
- cumulative energy: 1SV=0.56, 2=0.91, 3=0.98
- top-mode common-sign fraction: 0.75 (>0.8 => common reduction mode; ~0.5 => tradeoff)

## Single-intervention control authority (superposition NOT assumed)
- best single patch: **t3p3**  peak 3.4521e-02 -> 6.7222e-03 (**5.1x**, migrates to cell 25), dTBR +2.4%
- patches that reduce the peak (no significant off-peak rise): 5
- tradeoff patches (lower peak, raise another): 6

## OUTCOME
**1. HIDDEN STRUCTURED CONTROL: single sector patches significantly reduce the peak (placement-dependent)**

Note: the absolute peakX is baseline-confounded by the ~2.55 cm shield-outside offset; the placement-dependent STRUCTURE (which sectors reduce vs backfire) is the robust result. There is no WW-vs-no-WW mismatch. Combined designs and clean derivatives still require a consistently baselined OpenMC run.
