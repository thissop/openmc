"""Rigorous tests for the adjoint -> shield-placement pipeline (adjoint_placement.py).

Covers: plasma-source rasterization, the contributon (kills the central-vacuum
artifact), voxel->angle mapping, the (theta,phi) placement priority (shape, peak
location, periodicity, normalization), and the INTEGRATION with ThicknessField (the
optimizer's control field) -- the peak feeds gaussian_bump and the trade conserves the
envelope. A final integration test on the real cell-15 map asserts the placement peak
points at the coil (skipped if the committed map is absent).
"""
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import adjoint_placement as ap  # noqa: E402
from thickness_field import ThicknessField  # noqa: E402


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
MESH = dict(ll=np.array([-1000.0, -1000.0, -400.0]),
            ur=np.array([1000.0, 1000.0, 400.0]),
            dim=np.array([20, 20, 16]))


def _empty_map():
    d = dict(lower_left=MESH["ll"], upper_right=MESH["ur"], dimension=MESH["dim"],
             importance=np.zeros(tuple(MESH["dim"])), response="flat",
             coil_cells=np.array([15]))
    return d


def _set_voxel(d, x, y, z, val):
    xc, yc, zc = ap._voxel_centers(d["lower_left"], d["upper_right"], d["dimension"])
    i = int(np.argmin(np.abs(xc - x))); j = int(np.argmin(np.abs(yc - y)))
    k = int(np.argmin(np.abs(zc - z)))
    d["importance"][i, j, k] += val
    return xc[i], yc[j], zc[k]


# --------------------------------------------------------------------------- #
# plasma source rasterization
# --------------------------------------------------------------------------- #
def test_plasma_source_conserves_weight_and_lands_in_annulus():
    # a ring of plasma points at R=600 cm (scale=1), varying phi, z=0
    n = 200
    phi = np.linspace(0, 2 * np.pi, n, endpoint=False)
    R = np.full(n, 600.0); Z = np.zeros(n); sg = np.full(n, 2.0)
    fm = dict(R=R, Z=Z, phi=phi, sqrtg=sg)
    S = ap.plasma_source_on_mesh(fm, MESH["ll"], MESH["ur"], MESH["dim"], scale=1.0)
    assert np.isclose(S.sum(), sg.sum()), "histogram must conserve total source weight"
    # all mass at R~600: nonzero voxels lie in an annulus, none in the far corners
    xc, yc, zc = ap._voxel_centers(MESH["ll"], MESH["ur"], MESH["dim"])
    X, Y, _ = np.meshgrid(xc, yc, zc, indexing="ij")
    Rv = np.hypot(X, Y)
    assert S[Rv < 400].sum() == 0.0, "no source near the axis"
    assert S[S > 0].size > 0 and (Rv[S > 0].min() > 400) and (Rv[S > 0].max() < 800)


def test_contributon_kills_central_vacuum():
    # Bare importance nonzero EVERYWHERE (incl. the optically-thin central vacuum,
    # where the real bare adjoint flux artefactually peaks). The plasma source S(r)
    # lives only in an annulus, so the contributon C=S*imp must vanish in the centre
    # and survive in the annulus -- independent of any voxel-edge alignment.
    d = _empty_map()
    d["importance"][:] = 1.0e6
    phi = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    fm = dict(R=np.full(400, 650.0), Z=np.zeros(400), phi=phi, sqrtg=np.ones(400))
    C = ap.contributon(d, fm, scale=1.0)
    xc, yc, zc = ap._voxel_centers(MESH["ll"], MESH["ur"], MESH["dim"])
    X, Y, _ = np.meshgrid(xc, yc, zc, indexing="ij")
    Rv = np.hypot(X, Y)
    Cimp = C["importance"]
    assert Cimp[Rv < 300].sum() == 0.0, "contributon must be zero in the central vacuum"
    assert Cimp[Rv > 400].sum() > 0.0, "contributon must survive in the plasma annulus"
    # and the plasma-source mask really is what carved it
    assert np.all((Cimp > 0) == (C["plasma_source"] > 0))


# --------------------------------------------------------------------------- #
# voxel -> angle mapping
# --------------------------------------------------------------------------- #
def test_importance_angles_known_voxel():
    d = _empty_map()
    _set_voxel(d, 500.0, 500.0, 0.0, 1.0)     # phi = 45 deg
    phi, theta, w, R0 = ap.importance_angles(d, R0=0.0)
    assert len(phi) == 1
    assert np.isclose(np.degrees(phi[0]), 45.0, atol=3.0)
    assert np.isclose(np.degrees(theta[0]), 0.0, atol=3.0)   # z=0 -> midplane


def test_importance_angles_auto_R0_is_weighted_mean_radius():
    d = _empty_map()
    _set_voxel(d, 800.0, 0.0, 0.0, 3.0)
    _set_voxel(d, 0.0, 400.0, 0.0, 1.0)
    _, _, _, R0 = ap.importance_angles(d, R0=None)
    assert np.isclose(R0, (3 * 800 + 1 * 400) / 4.0, atol=60.0)   # ~ within a voxel


# --------------------------------------------------------------------------- #
# placement priority field
# --------------------------------------------------------------------------- #
def test_priority_shape_matches_thicknessfield():
    d = _empty_map(); _set_voxel(d, 600.0, 0.0, 0.0, 1.0)
    tor = np.linspace(0, 360, 24, endpoint=False)
    pol = np.linspace(0, 360, 32, endpoint=False)
    tf = ThicknessField(nfp=4, toroidal_angles_deg=tor, poloidal_angles_deg=pol,
                        t_breeder0=100.0, t_shield0=20.0)
    res = ap.placement_priority(d, tor, pol, R0=400.0)
    assert res["priority"].shape == tf.shape == (24, 32)


def test_priority_peaks_at_known_toroidal_angle():
    d = _empty_map()
    # concentrate importance in the phi ~ 135 deg sector (x<0, y>0)
    for _ in range(1):
        _set_voxel(d, -500.0, 500.0, 0.0, 10.0)
    tor = np.linspace(0, 360, 72, endpoint=False)   # 5 deg resolution
    pol = np.linspace(0, 360, 72, endpoint=False)
    res = ap.placement_priority(d, tor, pol, R0=0.0, width_deg=10.0)
    dphi = ((res["phi0_deg"] - 135.0 + 180) % 360) - 180
    assert abs(dphi) <= 10.0, f"priority peak {res['phi0_deg']} not at 135 deg"


def test_priority_periodic_across_zero():
    d = _empty_map()
    _set_voxel(d, 700.0, 0.0, 0.0, 1.0)    # phi = 0 deg (wrap boundary)
    tor = np.linspace(0, 360, 72, endpoint=False)
    pol = np.linspace(0, 360, 72, endpoint=False)
    res = ap.placement_priority(d, tor, pol, R0=0.0, width_deg=10.0)
    dphi = ((res["phi0_deg"] - 0.0 + 180) % 360) - 180
    assert abs(dphi) <= 10.0, "peak should sit at phi~0 without an edge artifact"


def test_priority_normalized_to_unity():
    d = _empty_map(); _set_voxel(d, 600.0, 100.0, 50.0, 1.0)
    res = ap.placement_priority(d, np.linspace(0, 360, 36, endpoint=False),
                                np.linspace(0, 360, 36, endpoint=False), R0=300.0)
    assert np.isclose(res["priority"].max(), 1.0)
    assert res["priority"].min() >= 0.0


def test_coil_angles_known():
    phi, th = ap.coil_angles((-500.0, 0.0, 0.0), R0=0.0)
    assert np.isclose(phi, 180.0, atol=1e-6)


# --------------------------------------------------------------------------- #
# INTEGRATION with the optimizer control field
# --------------------------------------------------------------------------- #
def test_gaussian_bump_at_priority_peak_and_envelope_conserved():
    d = _empty_map(); _set_voxel(d, -500.0, 500.0, 120.0, 5.0)
    tor = np.linspace(0, 360, 48, endpoint=False)
    pol = np.linspace(0, 360, 48, endpoint=False)
    tf = ThicknessField(nfp=4, toroidal_angles_deg=tor, poloidal_angles_deg=pol,
                        t_breeder0=100.0, t_shield0=20.0, t_breeder_min=10.0)
    res = ap.placement_priority(d, tor, pol, R0=300.0)
    # feed the recovered target into the optimizer's localized-thickening primitive
    delta = tf.gaussian_bump(res["theta0_deg"], res["phi0_deg"], amp=1.0, width_deg=30.0)
    assert delta.shape == tf.shape
    # the bump's peak must sit at the priority target. Toroidally the bump is
    # replicated every field period (nfp-fold stellarator symmetry), so it matches
    # phi0 MODULO 360/nfp -- which still covers the coil (all periods are equivalent).
    k = np.unravel_index(np.argmax(delta), delta.shape)
    period = 360.0 / tf.nfp
    dphi = (np.degrees(tf.tor[k[0]]) - res["phi0_deg"]) % period
    dphi = min(dphi, period - dphi)
    assert dphi <= 15.0, f"bump toroidal peak off by {dphi} deg (mod {period})"
    # poloidal is 2pi-periodic (no nfp folding) -> direct match
    dth = (np.degrees(tf.pol[k[1]]) - res["theta0_deg"]) % 360.0
    dth = min(dth, 360.0 - dth)
    assert dth <= 15.0, f"bump poloidal peak off by {dth} deg"
    # the breeder<->shield trade must conserve the fixed envelope
    t_s, t_b = tf.trade(delta)
    assert np.allclose(t_s + t_b, tf.t_s0 + tf.t_b0)
    assert np.all(t_b >= tf.t_b_min - 1e-9)


# --------------------------------------------------------------------------- #
# END-TO-END on the real committed map (the physics check)
# --------------------------------------------------------------------------- #
REAL_MAP = os.path.join(os.path.dirname(HERE), "..", "data", "adjoint", "coil15_fil_flat.npz")
REAL_FLUX = os.path.join(os.path.dirname(HERE), "..", "data", "qh_freestream_fluxmap.npz")
COIL15 = (-891.95, -743.44, -74.47)


@pytest.mark.skipif(not (os.path.exists(REAL_MAP) and os.path.exists(REAL_FLUX)),
                    reason="committed cell-15 adjoint map / fluxmap not present")
def test_real_contributon_placement_points_at_coil15():
    field = ap.contributon(REAL_MAP, REAL_FLUX, scale=100.0)
    tor = np.linspace(0, 360, 48, endpoint=False)
    pol = np.linspace(0, 360, 48, endpoint=False)
    res = ap.placement_priority(field, tor, pol)
    cphi, _ = ap.coil_angles(COIL15, res["R0"])
    dphi = ((res["phi0_deg"] - cphi + 180) % 360) - 180
    assert abs(dphi) <= 45.0, (
        f"placement peak phi0={res['phi0_deg']:.0f} not aligned with coil-15 "
        f"phi={cphi:.0f} (|dphi|={abs(dphi):.0f} deg)")


@pytest.mark.skipif(not (os.path.exists(REAL_MAP) and os.path.exists(REAL_FLUX)),
                    reason="committed cell-15 adjoint map / fluxmap not present")
def test_contributon_beats_bare_flux_on_real_map():
    """The contributon must point at the coil markedly better than the bare flux."""
    coil = np.array(COIL15); u = coil / np.linalg.norm(coil)
    d = np.load(REAL_MAP)
    ll, ur, dim = d["lower_left"], d["upper_right"], d["dimension"]
    xc, yc, zc = ap._voxel_centers(ll, ur, dim)
    X, Y, Z = np.meshgrid(xc, yc, zc, indexing="ij")

    def cos_to_coil(W):
        t = W.sum(); cen = np.array([(X * W).sum(), (Y * W).sum(), (Z * W).sum()]) / t
        return float(cen @ coil / (np.linalg.norm(cen) * np.linalg.norm(coil)))

    bare = cos_to_coil(d["importance"].astype(float))
    C = ap.contributon(REAL_MAP, REAL_FLUX, scale=100.0)["importance"]
    con = cos_to_coil(C)
    assert con > 0.7, f"contributon centroid cos(coil)={con:.2f} should point at coil"
    assert con > bare + 0.3, f"contributon ({con:.2f}) must beat bare flux ({bare:.2f})"
