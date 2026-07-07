#!/usr/bin/env python
"""Analytic Schwartz NWL oracle (via anarrima) for a square-cross-section torus with a parabolic
(1-r^2) plasma and a toroidal field. Emits the per-mode inboard/outboard MIDPLANE NWL and the
+43% / -22% reference the native SPF TokamakSource / StellaratorSource must reproduce in transport
(the capstone comparison). Mirrors anarrima/examples/square_torus.py exactly.

Local; needs anarrima importable (/Users/tkiker/Documents/GitHub/anarrima) + jax.
"""
import sys
import numpy as np

sys.path.insert(0, "/Users/tkiker/Documents/GitHub/anarrima/src")
from anarrima import reduced_integrals  # noqa: E402
import jax.numpy as jnp  # noqa: E402
from jax import jit  # noqa: E402

g_HI = jit(jnp.vectorize(reduced_integrals.g_HI))
g_HA = jit(jnp.vectorize(reduced_integrals.g_HA))
g_HAa = jit(jnp.vectorize(reduced_integrals.g_HAa))
g_HBa = jit(jnp.vectorize(reduced_integrals.g_HBa))

# --- geometry + plasma, identical to anarrima's square_torus.py ---
u_inboard, w_outboard, z_topbot = 0.4, 1.6, 0.6
minor_radius, major_radius = 0.5, 1.0
nr = nz = 80
ar = np.linspace(0.5, 1.5, nr); az = np.linspace(-0.5, 0.5, nz)
dr = (ar.max() - ar.min()) / nr; dz = (az.max() - az.min()) / nz
P, Zg = np.meshgrid(ar, az)
rho2 = (P - major_radius) ** 2 + Zg ** 2
inl = rho2 < minor_radius ** 2
strength_grid = np.where(inl, 1.0 - rho2 / minor_radius ** 2, 0.0) * dr * dz
p = P[inl]; zz = Zg[inl]; strength = strength_grid[inl]

phi_i = np.arccos(u_inboard / p)
phi_o = np.arccos(u_inboard / p) + np.arccos(u_inboard / w_outboard)


def irr_in(z, kind):
    zr = z - zz
    g = {"I": g_HI(p, zr, u_inboard, phi_i),
         "A": g_HA(p, zr, u_inboard, phi_i),
         "B": g_HBa(p, zr, u_inboard, 0, 0, phi_i)}[kind]
    return float(np.dot(np.asarray(g), strength))


def irr_out(z, kind):
    zr = zz - z
    g = {"I": -g_HI(p, zr, w_outboard, phi_o),
         "A": -g_HAa(p, zr, w_outboard, 0, 0, phi_o),
         "B": -g_HBa(p, zr, w_outboard, 0, 0, phi_o)}[kind]
    return float(np.dot(np.asarray(g), strength))


def main():
    z0 = 0.0
    # mode combination from square_torus.py: fI=(1/3)I, fa=(3/4)A, fb=(2/3)B
    fI_in = (1 / 3) * irr_in(z0, "I"); fI_out = (1 / 3) * irr_out(z0, "I")
    perp_in = (2 / 3) * (3 / 4) * irr_in(z0, "A")     # A rescaled x2/3 for fixed total power
    perp_out = (2 / 3) * (3 / 4) * irr_out(z0, "A")
    par_in = (2 / 3) * irr_in(z0, "B"); par_out = (2 / 3) * irr_out(z0, "B")
    print("=== anarrima oracle: square torus, parabolic (1-r^2) plasma, toroidal field, midplane z=0 ===")
    print(f"inboard  : unpol {fI_in:.4g} | perp/A {perp_in:.4g} ({100*(perp_in/fI_in-1):+.1f}%) "
          f"| par/BC {par_in:.4g} ({100*(par_in/fI_in-1):+.1f}%)")
    print(f"outboard : unpol {fI_out:.4g} | perp/A {perp_out:.4g} ({100*(perp_out/fI_out-1):+.1f}%) "
          f"| par/BC {par_out:.4g} ({100*(par_out/fI_out-1):+.1f}%)")
    print("\nSchwartz reference: perp/A ~ +43% inboard, -22% outboard; par/BC the mirror.")
    return dict(inboard=dict(unpol=fI_in, perp=perp_in, par=par_in),
                outboard=dict(unpol=fI_out, perp=perp_out, par=par_out))


if __name__ == "__main__":
    main()
