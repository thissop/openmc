# One-period native FW-CADIS

This is the real-DAGMC continuation of Path A after the canonical CSG smoke
test passed. It fixes two defects in the archived prototype: the target is
provided as an `openmc.Tallies` collection, and the requested ray count is set
after `convert_to_random_ray()` so it is not overwritten by OpenMC's geometric
default.

The generated CASMO-2 infinite-medium library is stored once as
`mgxs_casmo2.h5` beside the staged script and reused by later Slurm runs. The
preflight record proves the final energy mode, target list, ray count, and
random-ray configuration immediately before transport.

Random-ray fixed-source mode accepts only `IndependentSource`. The native
`StellaratorSource` therefore cannot be passed through directly. This runner
reads the same `spf_fluxmap_v1` payload, weights grid locations by the native
default emission law `(1-rho^2) sqrt(g)`, deterministically resamples at most
2,000 positions, folds them into the first field period, and supplies them as
weighted isotropic 14.06-MeV `Point` sources. Random ray locates each point's
source region, so the fluxmap distribution is retained without relying on a
spatial distribution type the solver would reject. Polarization is outside the
scalar random-ray solve.

Only the eight physical-density materials use the infinite-medium generator.
The `gap` and `Vacuum` hydrogen placeholders have negligible optical depth in
the device; running them inside the generator's 1000-m reflective box causes
maximum-event histories. They therefore receive explicit two-group
`1e-12 cm^-1` absorption cross sections and zero scatter, which behaves as void
on the device length scale without poisoning the MGXS generation step.

The current 10 x 10 x 10 mesh is only an integration pilot. It cannot resolve
the 40-cm shield through the full device bounding box. Production use still
requires roughly eight radial shield layers (about 5 cm each) in the ParaStell
geometry, followed by a mesh-resolution study.
