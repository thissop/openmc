# Polarization port onto Ethan's branch — plan & status

Ethan's `stellarator_source` branch already has `TokamakSource` and `StellaratorSource` (Miller /
Fourier geometry, **no polarization**). Your `TokamakSource` is the *same PR#3999 base*, so the
polarization layer grafts on rather than replacing anything. We add the SPF layer + a B-hat backend;
we do **not** add a second source class, and your `spf_core_integration.patch` does NOT `git apply`
here (different base) — this is a manual, source-mapped graft.

Source of truth for every hunk: `openmc/spf_prototype/native_spf/` (Python `source.py` +181,
C++ `source.cpp` +208, `include/openmc/source.h` +43, `openmc/__init__.py` +2).

## Status
- [x] `spin_fractions_to_abc(d_plus,d_zero,d_minus,t_plus,t_minus)` -> `openmc/source.py` (verified: non-pol=1/3,1/3,1/3; sum=1). Auto-exported (`source.py` has no `__all__`; `__init__` does `import *`).
- [x] **TokamakSource Python** — DONE: `polarization`/`field_model`/`safety_factor`/`field_sign` kwargs
      + assignments + property setters (normalize / convert-from-dict / renormalize-with-warn) +
      `_validate` cross-field check (pitched requires safety_factor) + `populate_xml_element`
      serialization + `from_xml_element` parse + `cls(...)` kwargs. File parses; imports present.
      Functional roundtrip test pending a build (needs compiled openmc to import).
- [ ] **TokamakSource C++** (`src/source.cpp`, `include/openmc/source.h`) — Eq.2 coefficient derivation
      (`w_perp=3/4 a`, `w_par=2/3 b+1/3 c`, `p_perp=a/eta`), `sample_polarized_direction(bhat, seed)`
      (stratified P2 shapes: sin^2 and 1/4+3/4cos^2, inverse-CDF), `field_direction(r,alpha,phi)`
      (toroidal or q-pitched), gated on `polarized_`; unpolarized path bit-identical to Ethan's.
- [ ] **`openmc/__init__.py`** — export `spin_fractions_to_abc` (StellaratorSource already exported).
- [ ] **StellaratorSource polarization** — same kwargs; B-hat from the flux-surface field. Two backends:
      (a) fluxmap-B (your current path), (b) **on-the-fly B from iota+lambda** (Ethan's ask) — the wout
      already carries `iotaf`+`lmns`; the validation harness's `read_wout` surfaces them.
- [ ] **Tests** in Ethan's style (`tests/unit_tests/test_source_tokamak.py`): polarized-mode corner
      values (eta = 2/3, 1, 2/3, 1/3), B/C-share-shape, isotropic==unpolarized, roundtrip.

## Quick win (Ethan flagged, ~immediate)
Tokamak spin-pol with a **real** b-hat: his TokamakSource ran capstone with b-hat=toroidal; add the
q-profile + derivatives so b-hat carries the poloidal tilt (`field_model='pitched'`, `safety_factor`).
Your Python layer already has these kwargs — the work is the C++ `field_direction` pitched branch.

## Commit discipline
Accumulate all of the above into ONE clean commit off `eepeterson/stellarator_source`, then PR. The
`spf_validation/` harness may ride along as supporting evidence or be split out — decide at PR time.
