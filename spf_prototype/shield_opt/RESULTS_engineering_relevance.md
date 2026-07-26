# Engineering-relevance (blanket-fit) filter + geometric selection effect

**Question (user):** for the 200-300-device magnet-side study, keep only devices where a real
~1 m+ radial build (blanket+shield) physically fits between plasma and coils — exclude the tiny
(~16 cm) plasma-coil-standoff devices — and **quantify the geometric selection effect** the filter
imposes, since it confounds any law derived on the filtered set.

## Method
Reactor scaling = ARIES-CS (minor radius a=1.704 m, Kappel convention). Radial build = the validated
8-layer corrected build ≈ **1.29 m** (FW 3.2 + mult 2 + breeder 50 + backwall 4 + shield 40 + gap 2 +
VV 25 + TS 3 cm); also a relaxed **1.0 m** cut. Two criteria, because "does a blanket fit" is not
scaling-neutral — and **that ambiguity IS the selection effect**:
- **(A) As-shipped:** the QUASR coil solution's actual gap, `gap = (d_min/a)·a_reactor`. "Room with THESE coils."
- **(B) Min-achievable:** Kappel `L_REGCOIL` (closest coils CAN sit) / the `L_gradB` proxy. The floor.

## Results
**(A) As-shipped QUASR coils (366 QA/QH):**
- gap ≥ 1.29 m: **327/366 keep (89%)** — QA 79%, QH 97%.
- gap ≥ 1.0 m: 344/366 (94%).
- **Selection effect — drops COMPACT devices:** aspect median *kept* = **12.0** vs *dropped* = **3.3–4.0**.
  High-aspect thin plasmas sit far from QUASR's generously-placed coils → pass; low-aspect **compact
  (reactor-attractive)** devices get squeezed out. The as-shipped filter **biases toward high aspect**,
  i.e. toward reactor-*un*attractive geometry. The ~11% dropped are the user's "16 cm" devices
  (d_min/a < 0.76).

**(B) Min-achievable separation (Kappel-45, incl QI/tokamak/stellarator):**
- L_REGCOIL ≥ 1.29 m: **28/45 (62%)** — QA 3/4, QH 7/10, **QI 4/4**, tokamak 2/5.
- **Selection effect — favors low nfp / QI:** corr(nfp, L_REGCOIL) = **−0.472**, corr(nfp, L_gradB) = −0.464.
  Small nfp ⇒ more blanket room (Kappel's reactor-attractive trend). The min-achievable filter
  **preferentially keeps low-nfp / QI / QA**.

## The finding
The two criteria bias in **different directions** (as-shipped → high aspect; min-achievable → low nfp/QI),
so the blanket-fit filter is **not population-neutral**, and the *criterion choice sets the bias*. Any
"geometry → magnet-shieldability" law fit on the filtered set inherits this. For the magnet-side study we
adopt criterion (A) as operative (we study devices with THEIR coils) but **must report results
stratified/controlled for aspect and nfp**, and cross-check against (B) (min-achievable / L_gradB) so the
conclusion is not an artifact of who survived the cut. This is the confound the free-streaming zoo law
(RESULTS_zoo_coil_law.md) already partly addressed (standoff is NOT a concentration driver; elongation is).

## Artifacts
- `shield_opt/data/engineering_relevant_devices.csv` — 366 QA/QH with gap_m + blanket_fit flags (327 pass @1.29 m).
- `figs/engineering_relevance.png` — gap distribution (as-shipped) + nfp→separation (min-achievable).
- QUASR-366 is QA/QH only; **QI blanket-fit devices** (median L_REGCOIL 3.58 m → fit easily) come from
  Kappel-45 + ConStellaration. Building the full 200-300 QI/QA/QH engineering-relevant set needs d_min
  computed for MORE devices (incl QI coils) — the scaled magnet-side pipeline's selection step.
