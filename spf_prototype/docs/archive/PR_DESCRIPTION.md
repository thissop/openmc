This adds spin-polarized D-T emission on top of Dr. Ethan Peterson's `stellarator_source`
branch. Fusion neutrons are born with their direction sampled about the local magnetic field
`B̂` from the `(a, b, c)` collision-mode fractions (Schwartz 2025), through a single shared
sampler that also drives Dr. Peterson's `TokamakSource`. With no `polarization` argument
nothing changes, so it is fully backward compatible (same isotropic path, same RNG draw
order).

What's in it:

- `StellaratorSource`: `polarization` with `B̂` reconstructed on the fly from `iota` and
  `lambda`. `from_vmec` reads them from the wout (half and full mesh); `from_desc` pulls them
  from the DESC Fourier-Zernike basis and iota profile.
- `TokamakSource`: `polarization` with a toroidal or q-pitched field.
- `spin_fractions_to_abc` helper, user-guide and API docs, and 56 unit tests.

Validation: the compiled sampler recovers the analytic angular moments (`⟨cos²θ⟩` = 1/5, 1/3,
7/15) for both sources, and the reconstructed field matches VMEC's exported `bsupumnc`
(~1e-4°) and DESC's `compute("B")` (sub-degree, converging with resolution). The full DESC to
sampler chain is checked end-to-end on a QA and a QH device.

I also compared evaluating the stellarator geometry and magnetic field analytically (on the
fly) against interpolating them from a precomputed grid, looking at speed, memory, and
sampling accuracy for a uniform grid and an error-driven adaptive grid.
