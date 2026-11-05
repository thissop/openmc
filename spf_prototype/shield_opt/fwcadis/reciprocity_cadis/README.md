# Reciprocity-CADIS fallback

This is Path B from the 2026-08-21 handoff. It runs continuous-energy forward
transport from an isotropic point cloud sampled along the first five QH coil
filaments. The reciprocal flux map approximates the adjoint importance for a
fast-flux response in those coils. The script then applies the same inverse,
volume-normalized, one-half normalization used by OpenMC's CADIS generator and
writes `ww_reciprocity_cadis.h5`.

Bins with no samples or relative error above 0.5 are disabled with OpenMC's
negative-bound convention. Rejected bins are excluded before normalization so
a single noisy low-flux bin cannot collapse all useful windows. This path is a
spatial reciprocity approximation, not a replacement for a fully energy-adjoint
solve; its value is that it retains the real continuous-energy DAGMC transport
while avoiding MGXS and random-ray fragility.

The batch script is staged and submitted only through the AdvancedAgents
Ginsburg protocol. Its default pilot is 400,000 histories on a 16 x 16 x 8
mesh. The `result.json` reliable-bin fraction determines whether a finer or
larger pilot is justified.
