"""Phase D scaffolding: the smooth t_shield(theta,phi) control field and the
fixed-envelope breeder<->shield trade that feeds ParaStell's per-layer thickness_matrix.

Design variable is a low-dimensional smooth field (Fourier on the angular torus), NOT a
per-cell thickness (Monte-Carlo noise + meshing pathology make that unstable). The trade is
strictly volume-neutral in the radial direction: whatever thickness we ADD to the shield we
REMOVE from the breeder at the same (theta,phi), so the outer envelope (and the coils) never
move. Constraints: 0 <= delta <= delta_max, and breeder never thinner than t_breeder_min.
"""
from __future__ import annotations

import numpy as np


class ThicknessField:
    def __init__(self, nfp, toroidal_angles_deg, poloidal_angles_deg,
                 t_breeder0, t_shield0, delta_max=None, t_breeder_min=10.0):
        self.nfp = int(nfp)
        self.tor = np.radians(np.asarray(toroidal_angles_deg, float))
        self.pol = np.radians(np.asarray(poloidal_angles_deg, float))
        # ParaStell convention: thickness_matrix is (n_toroidal, n_poloidal)
        self.TOR, self.POL = np.meshgrid(self.tor, self.pol, indexing="ij")
        self.shape = self.TOR.shape
        self.t_b0 = float(t_breeder0)
        self.t_s0 = float(t_shield0)
        self.t_b_min = float(t_breeder_min)
        # max shield addition is bounded by how much breeder we can spare
        self.delta_max = float(delta_max) if delta_max is not None else (self.t_b0 - self.t_b_min)

    # ---- basis -----------------------------------------------------------
    def fourier_basis(self, M=2, N=1):
        """Stellarator-symmetric cos + sin modes cos/sin(m*theta - n*nfp*phi),
        m in 0..M, n in -N..N. Returns (K, nTor, nPol)."""
        basis = []
        for m in range(M + 1):
            for n in range(-N, N + 1):
                arg = m * self.POL - n * self.nfp * self.TOR
                basis.append(np.cos(arg))
                if not (m == 0 and n == 0):
                    basis.append(np.sin(arg))
        return np.asarray(basis)

    # ---- delta field (bounded, non-negative) -----------------------------
    def delta(self, coeffs, basis):
        """delta(theta,phi) = delta_max * sigmoid(sum_k c_k basis_k), then clipped so the
        breeder stays >= t_breeder_min. Smooth, in [0, delta_max]. Zero field <-> coeffs that
        make sigmoid~0 (large negative bias c_0)."""
        g = np.tensordot(np.asarray(coeffs, float), basis, axes=(0, 0))
        d = self.delta_max / (1.0 + np.exp(-g))
        return np.minimum(d, self.t_b0 - self.t_b_min)

    def gaussian_bump(self, theta0_deg, phi0_deg, amp, width_deg=30.0):
        """A localized smooth shield-thickening bump toward a target (theta0,phi0) --
        the natural warm-start once Phase C hands us a coil hotspot location. Periodic."""
        def dperiodic(a, b):
            d = (a - b + np.pi) % (2 * np.pi) - np.pi
            return d
        w = np.radians(width_deg)
        dth = dperiodic(self.POL, np.radians(theta0_deg))
        dph = dperiodic(self.nfp * self.TOR, self.nfp * np.radians(phi0_deg)) / self.nfp
        d = amp * np.exp(-(dth ** 2 + dph ** 2) / (2 * w ** 2))
        return np.minimum(d, self.t_b0 - self.t_b_min)

    # ---- the fixed-envelope trade ---------------------------------------
    def trade(self, delta):
        """shield += delta, breeder -= delta. Returns ParaStell-ready
        (shield_thickness_matrix, breeder_thickness_matrix). Sum is invariant."""
        t_s = self.t_s0 + delta
        t_b = self.t_b0 - delta
        assert np.all(t_b >= self.t_b_min - 1e-9), "breeder below minimum"
        assert np.allclose(t_s + t_b, self.t_s0 + self.t_b0), "envelope not conserved"
        return t_s, t_b

    def radial_build_dict(self, delta, base_dict):
        """Insert the traded breeder/shield thickness matrices into a baseline radial-build
        dict (leaves first_wall/back_wall/vacuum_vessel untouched)."""
        t_s, t_b = self.trade(delta)
        rb = {k: dict(v) for k, v in base_dict.items()}
        rb["breeder"]["thickness_matrix"] = t_b
        rb["shield"]["thickness_matrix"] = t_s
        return rb

    # ---- regularizers ----------------------------------------------------
    def smoothness(self, delta):
        """Mean squared angular gradient -- the manufacturability / MHD-gentleness penalty."""
        dph = np.gradient(delta, axis=0)
        dth = np.gradient(delta, axis=1)
        return float(np.mean(dth ** 2 + dph ** 2))

    def added_shield_volume(self, delta):
        """Proxy for how much shield material we spent (area-integral of delta)."""
        return float(np.mean(delta))  # per-unit-area; scale by wall area downstream
