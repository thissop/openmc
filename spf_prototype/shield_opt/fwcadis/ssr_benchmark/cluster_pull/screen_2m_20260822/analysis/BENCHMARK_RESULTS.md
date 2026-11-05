# SSR variance-reduction benchmark

Surface-source responses are multiplied by the bank's total starting weight per master source history. Bias therefore tests the full deep-freeze replay, including the deliberate black-inside approximation.

- Surface-bank sites: 8,291,969
- Surface-bank total weight: 8226805.9
- Surface weight per master history: 2.0567015

| mode | coil bins > 0 | median rel. err. | median FOM | transport (s) | tracks/source | RMS absolute bias | RMS pattern error | master-compatible bins |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| analog | 1/5 | 1 | 0.00341 | 293 | 1 | 1 | 1.43 | 0/5 |
| magic | 3/5 | 1 | 0.00313 | 319 | 1.16 | 0.999 | 0.744 | 0/5 |
| reciprocity | 1/5 | 1 | 0.00348 | 287 | 1.06 | 1 | 1.43 | 0/5 |

The pilot is a split-safety and wiring check. A production verdict requires all five coil bins to be positive and stable across independent seeds.
