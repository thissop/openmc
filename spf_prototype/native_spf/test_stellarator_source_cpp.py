#!/usr/bin/env python
"""Unit tests for the native C++ ``StellaratorSource`` (native_spf/stellarator_source.cpp).

A self-contained NumPy MIRROR of the exact C++ cascade math -- the marginal ->
conditional CDF construction, the searchsorted(side="right") cell selection, the
in-cell jitter, and the trilinear interpolation (theta->zeta->rho term order,
2*pi-periodic theta/zeta, clamped rho) -- is reproduced here and checked, on a
small synthetic fluxmap, against the VALIDATED reference algorithm
``ConformalStellaratorSampler`` (python/stellarator_source.py):

  * the C++ CDF cascade is built the same way as the reference (bit-consistent
    to ~1e-12 despite the sequential-loop vs numpy-vectorized summation);
  * given the SAME uniform draw sequence, the mirror and the reference select
    the IDENTICAL discrete cells (i, j, k)  -- "same seeds -> same cells";
  * the trilinear-interpolated positions and grid b-hat match the reference to
    ~1e-12 (proves the interp transcription and the theta/zeta period wrap);
  * every birth weight is exactly 1 (rejection-free cascade) and every birth is
    inside the plasma (rho <= 1 + jitter tolerance);
  * the portable spf_fluxmap_v1 .meta/.bin round-trips bit-exactly.

The mirror reflects the C++ STRUCTURE (not the reference), so a transcription bug
in stellarator_source.cpp that this test also mirrors would still be caught by
the comparison against the independent reference. The bit-for-bit C++ prn-stream
parity is a separate build-time gate (see README_stellarator.md), exactly as for
the SPF direction sampler.

Run:  pytest -q native_spf/test_stellarator_source_cpp.py
  or: python native_spf/test_stellarator_source_cpp.py     (pure numpy; no build)
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PYDIR = HERE.parent / "python"
sys.path.insert(0, str(PYDIR))

from stellarator_source import ConformalStellaratorSampler, circular_torus_fluxmap  # noqa: E402
from quasr_fluxmap import write_fluxmap_bin  # noqa: E402


# Small synthetic fluxmap with a KNOWN analytic sqrt(g) (circular torus).
NR, NT, NZ = 8, 16, 12
R0, A = 6.0, 2.0
S_of_rho = lambda r: 1.0 - r ** 2  # noqa: E731  (reference default profile)


def make_map():
    fm, _ = circular_torus_fluxmap(R0=R0, a=A, nr=NR, nt=NT, nz=NZ, nfp=1)
    return fm


# ===========================================================================
# NumPy mirror of native_spf/stellarator_source.cpp
# ===========================================================================

def mirror_build_cdfs(fm, S_fn):
    """Mirror StellaratorSource::build_cdf_cascade using the C++ SEQUENTIAL
    accumulation (triple loops), NOT numpy's vectorized sum. Comparing this
    against the reference's numpy CDFs is the real test that the C++ summation /
    axis order / per-slice normalization reproduce the reference."""
    rho = np.asarray(fm["rho"], float)
    sqrtg = np.asarray(fm["sqrtg"], float)
    nr, nt, nz = sqrtg.shape
    Svec = np.maximum(np.array([S_fn(r) for r in rho], float), 0.0)
    w = Svec[:, None, None] * sqrtg  # per-cell weight S(rho)*sqrt(g)

    # marginal p(rho): cdf_rho[nr+1]
    cdf_rho = np.zeros(nr + 1)
    for i in range(nr):
        m = 0.0
        for j in range(nt):
            for k in range(nz):
                m += w[i, j, k]
        cdf_rho[i + 1] = cdf_rho[i] + m
    cdf_rho /= cdf_rho[-1]

    # conditional p(zeta|rho): cdf_zeta[nr, nz+1]
    cdf_zeta = np.zeros((nr, nz + 1))
    for i in range(nr):
        for k in range(nz):
            s = 0.0
            for j in range(nt):
                s += w[i, j, k]
            cdf_zeta[i, k + 1] = cdf_zeta[i, k] + s
        dz = cdf_zeta[i, nz]
        cdf_zeta[i] /= dz if dz > 0 else 1.0

    # conditional p(theta|rho,zeta): cdf_theta[nr, nt+1, nz]
    cdf_theta = np.zeros((nr, nt + 1, nz))
    for i in range(nr):
        for k in range(nz):
            for j in range(nt):
                cdf_theta[i, j + 1, k] = cdf_theta[i, j, k] + w[i, j, k]
            dt = cdf_theta[i, nt, k]
            cdf_theta[i, :, k] /= dt if dt > 0 else 1.0

    return cdf_rho, cdf_zeta, cdf_theta


def cells_from_cdfs(cdf_rho, cdf_zeta, cdf_theta, u_rho, u_zeta, u_theta):
    """Mirror the cascade cell selection: searchsorted(side="right") - 1, clamped.
    (side="right" == std::upper_bound == the C++ upper_bound_strided.)"""
    nr = cdf_rho.size - 1
    nz = cdf_zeta.shape[1] - 1
    nt = cdf_theta.shape[1] - 1
    i = np.clip(np.searchsorted(cdf_rho, u_rho, side="right") - 1, 0, nr - 1)
    k = np.clip(np.array([np.searchsorted(cdf_zeta[ii], uu, side="right") - 1
                          for ii, uu in zip(i, u_zeta)]), 0, nz - 1)
    j = np.clip(np.array([np.searchsorted(cdf_theta[ii, :, kk], uu, side="right") - 1
                          for ii, kk, uu in zip(i, k, u_theta)]), 0, nt - 1)
    return i, j, k


def _wrap_index(x, grid, period):
    n = grid.size
    d = period / n
    t = (x - grid[0]) / d
    ft = np.floor(t)
    i0 = (ft.astype(int) % n + n) % n
    frac = t - ft
    i1 = (i0 + 1) % n
    return i0, i1, frac


def mirror_position_and_bhat(fm, i, j, k, u_jr, u_jt, u_jz):
    """Mirror the jitter + trilinear interp of StellaratorSource::sample."""
    rho = np.asarray(fm["rho"], float)
    theta = np.asarray(fm["theta"], float)
    zeta = np.asarray(fm["zeta"], float)
    nr = rho.size
    drho = rho[1] - rho[0]
    dth = theta[1] - theta[0]
    dze = zeta[1] - zeta[0]
    # jitter (order rho, theta, zeta -- matches source.cpp)
    rr = rho[i] + (u_jr - 0.5) * drho
    tt = theta[j] + (u_jt - 0.5) * dth
    zz = zeta[k] + (u_jz - 0.5) * dze
    # brackets: rho clamped (side="left"/lower_bound), theta/zeta 2*pi-periodic
    ir = np.clip(np.searchsorted(rho, rr) - 1, 0, nr - 2)
    fr = (rr - rho[ir]) / drho
    jt0, jt1, ft = _wrap_index(tt, theta, 2 * np.pi)
    kz0, kz1, fz = _wrap_index(zz, zeta, 2 * np.pi)

    def interp(A):  # term order theta -> zeta -> rho
        c00 = A[ir, jt0, kz0] * (1 - ft) + A[ir, jt1, kz0] * ft
        c01 = A[ir, jt0, kz1] * (1 - ft) + A[ir, jt1, kz1] * ft
        c0 = c00 * (1 - fz) + c01 * fz
        d00 = A[ir + 1, jt0, kz0] * (1 - ft) + A[ir + 1, jt1, kz0] * ft
        d01 = A[ir + 1, jt0, kz1] * (1 - ft) + A[ir + 1, jt1, kz1] * ft
        d0 = d00 * (1 - fz) + d01 * fz
        return c0 * (1 - fr) + d0 * fr

    R = interp(np.asarray(fm["R"], float))
    Z = interp(np.asarray(fm["Z"], float))
    PH = interp(np.asarray(fm["phi"], float))
    x = R * np.cos(PH)
    y = R * np.sin(PH)
    br = interp(np.asarray(fm["BR"], float))
    bp = interp(np.asarray(fm["Bphi"], float))
    bz = interp(np.asarray(fm["BZ"], float))
    bx = br * np.cos(PH) - bp * np.sin(PH)
    by = br * np.sin(PH) + bp * np.cos(PH)
    bn = np.sqrt(bx ** 2 + by ** 2 + bz ** 2)
    bn = np.where(bn > 0, bn, 1.0)
    bhat = np.stack([bx / bn, by / bn, bz / bn], axis=-1)
    return np.column_stack([x, y, Z]), rr, bhat


def read_fluxmap_bin(stem):
    """Mirror StellaratorSource::read_fluxmap: parse .meta + unpack .bin."""
    meta = {}
    for line in Path(stem + ".meta").read_text().splitlines():
        if line.strip():
            key, val = line.split()
            meta[key] = val
    assert meta["magic"] == "spf_fluxmap_v1"
    nr, nt, nz = int(meta["nR"]), int(meta["ntheta"]), int(meta["nzeta"])
    raw = np.fromfile(stem + ".bin", dtype="<f8")
    off = 0

    def take(n):
        nonlocal off
        a = raw[off:off + n]
        off += n
        return a

    rho, theta, zeta = take(nr), take(nt), take(nz)
    n3 = nr * nt * nz
    take3 = lambda: take(n3).reshape(nr, nt, nz)  # noqa: E731
    sqrtg, R, Z, phi = take3(), take3(), take3(), take3()
    BR, Bphi, BZ = take3(), take3(), take3()
    assert off == raw.size, f"trailing bytes in {stem}.bin"
    return dict(rho=rho, theta=theta, zeta=zeta, sqrtg=sqrtg, R=R, Z=Z, phi=phi,
                BR=BR, Bphi=Bphi, BZ=BZ, nfp=int(meta["nfp"]))


# ===========================================================================
# Tests
# ===========================================================================

def test_binary_roundtrip():
    """write_fluxmap_bin -> read_fluxmap_bin recovers every array bit-exactly and
    proves the .bin layout the C++ reader expects."""
    fm = make_map()
    with tempfile.TemporaryDirectory() as d:
        stem = str(Path(d) / "syn_fluxmap")
        write_fluxmap_bin(stem, fm["rho"], fm["theta"], fm["zeta"], fm["sqrtg"],
                          fm["R"], fm["Z"], fm["phi"], fm["BR"], fm["Bphi"],
                          fm["BZ"], fm["nfp"], sign_sqrtg=1.0)
        rd = read_fluxmap_bin(stem)
    for key in ("rho", "theta", "zeta", "sqrtg", "R", "Z", "phi", "BR", "Bphi", "BZ"):
        assert np.array_equal(np.asarray(fm[key], "<f8"), rd[key]), f"{key} mismatch"
    assert rd["nfp"] == int(fm["nfp"])
    assert rd["sqrtg"].shape == (NR, NT, NZ)


def test_cdf_cascade_matches_reference():
    """The C++ (sequential-loop) CDF construction reproduces the reference's
    numpy-vectorized CDFs to ~1e-12: same axis order, same per-slice norm."""
    fm = make_map()
    smp = ConformalStellaratorSampler(fm, S_of_rho=S_of_rho)
    cdf_rho, cdf_zeta, cdf_theta = mirror_build_cdfs(fm, S_of_rho)
    assert np.allclose(cdf_rho, smp._cdf_rho, atol=1e-12, rtol=0)
    assert np.allclose(cdf_zeta, smp._cdf_zeta, atol=1e-12, rtol=0)
    assert np.allclose(cdf_theta, smp._cdf_theta, atol=1e-12, rtol=0)
    # Monotone, ends at 1 (the invariants the C++ hard-asserts).
    assert np.all(np.diff(cdf_rho) >= -1e-15)
    assert abs(cdf_rho[-1] - 1.0) < 1e-12


def test_same_seeds_same_cells():
    """Given the SAME uniform draw sequence, the mirror-of-C++ cascade selects
    the IDENTICAL cells as the reference _sample_cells."""
    fm = make_map()
    smp = ConformalStellaratorSampler(fm, S_of_rho=S_of_rho)
    n = 20_000
    # Reference draws rho(n), zeta(n), theta(n) in that order from the rng.
    ref_i, ref_j, ref_k = smp._sample_cells(n, np.random.default_rng(0))
    # Replicate the SAME draws, then run them through the mirror CDFs.
    rng = np.random.default_rng(0)
    u_rho, u_zeta, u_theta = rng.random(n), rng.random(n), rng.random(n)
    cdf_rho, cdf_zeta, cdf_theta = mirror_build_cdfs(fm, S_of_rho)
    mir_i, mir_j, mir_k = cells_from_cdfs(cdf_rho, cdf_zeta, cdf_theta,
                                          u_rho, u_zeta, u_theta)
    assert np.array_equal(ref_i, mir_i)
    assert np.array_equal(ref_j, mir_j)
    assert np.array_equal(ref_k, mir_k)


def test_positions_and_bhat_match_reference():
    """Full jitter + trilinear interp: mirror positions and grid b-hat match the
    reference sample() to ~1e-12 (validates the interp + period wrap)."""
    fm = make_map()
    smp = ConformalStellaratorSampler(fm, S_of_rho=S_of_rho)
    n = 20_000
    seed = 12345
    out = smp.sample(n, rng=np.random.default_rng(seed), jitter=True)
    # Reference draw order: rho-cell, zeta-cell, theta-cell, then jitter r,theta,z.
    rng = np.random.default_rng(seed)
    u_rho, u_zeta, u_theta = rng.random(n), rng.random(n), rng.random(n)
    u_jr, u_jt, u_jz = rng.random(n), rng.random(n), rng.random(n)
    cdf_rho, cdf_zeta, cdf_theta = mirror_build_cdfs(fm, S_of_rho)
    i, j, k = cells_from_cdfs(cdf_rho, cdf_zeta, cdf_theta, u_rho, u_zeta, u_theta)
    xyz, rr, bhat = mirror_position_and_bhat(fm, i, j, k, u_jr, u_jt, u_jz)
    ref_xyz = np.column_stack([out["x"], out["y"], out["z"]])
    assert np.allclose(xyz, ref_xyz, atol=1e-10, rtol=0)
    assert np.allclose(bhat, out["bhat"], atol=1e-10, rtol=0)


def test_unit_weight_and_inside_plasma():
    """Rejection-free cascade => birth weight exactly 1; every birth inside the
    plasma (rho <= 1 + jitter tolerance)."""
    fm = make_map()
    smp = ConformalStellaratorSampler(fm, S_of_rho=S_of_rho)
    out = smp.sample(50_000, rng=np.random.default_rng(7), jitter=True)
    assert float(np.abs(out["weight"] - 1.0).max()) < 1e-12
    Rmaj = np.hypot(out["x"], out["y"])
    rho_born = np.hypot(Rmaj - R0, out["z"]) / A
    assert rho_born.max() <= 1.0 + 1.0 / NR  # one jitter half-cell past the edge
    # b-hat is a unit vector everywhere.
    assert np.allclose(np.linalg.norm(out["bhat"], axis=1), 1.0, atol=1e-12)


def test_grid_bhat_is_toroidal_for_circular_torus():
    """On the analytic circular torus BR=BZ=0, Bphi>0, so the grid b-hat equals
    the pure toroidal phi-hat = (-sin PH, cos PH, 0) -- the field_model='toroidal'
    override must agree with field_model='fluxmap' here."""
    fm = make_map()
    smp = ConformalStellaratorSampler(fm, S_of_rho=S_of_rho)
    out = smp.sample(10_000, rng=np.random.default_rng(3), jitter=True)
    PH = out["phi"]
    phihat = np.column_stack([-np.sin(PH), np.cos(PH), np.zeros_like(PH)])
    assert np.allclose(out["bhat"], phihat, atol=1e-10, rtol=0)


# ===========================================================================
# Standalone runner (mirrors the style of test_stellarator_source.py)
# ===========================================================================

def main():
    checks = [
        ("spf_fluxmap_v1 .meta/.bin round-trip", test_binary_roundtrip),
        ("CDF cascade == reference (~1e-12)", test_cdf_cascade_matches_reference),
        ("same seeds -> same cells", test_same_seeds_same_cells),
        ("positions + b-hat == reference", test_positions_and_bhat_match_reference),
        ("unit weight + inside plasma", test_unit_weight_and_inside_plasma),
        ("grid b-hat == toroidal phi-hat", test_grid_bhat_is_toroidal_for_circular_torus),
    ]
    ok = True
    for name, fn in checks:
        try:
            fn()
            print(f"[V] {name:42s}  OK")
        except AssertionError as e:
            ok = False
            print(f"[V] {name:42s}  FAIL  {e}")
    print("\nRESULT:", "ALL C++-MIRROR GATES PASSED" if ok
          else "SOME GATES FAILED -- investigate")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
