import numpy as np
import pytest

import openmc
import openmc.stats

from tests.unit_tests import assert_sample_mean


def circular_tokamak_modes(r0=600.0, a=150.0, n_rho=17):
    """Axisymmetric circular torus expressed in stellarator Fourier form."""
    rho = np.linspace(0.0, 1.0, n_rho)
    mode_m = [0, 1]
    mode_n = [0, 0]
    rmnc = np.column_stack([np.full(n_rho, r0), a * rho])
    zmns = np.column_stack([np.zeros(n_rho), a * rho])
    return rho, mode_m, mode_n, rmnc, zmns


def rotating_ellipse_modes(r0=600.0, a=150.0, e=50.0, n_rho=17, nfp=5):
    """Classic rotating-ellipse stellarator:

    R = R0 + a*rho*cos(theta) + e*rho*cos(theta - nfp*zeta)
    Z =      a*rho*sin(theta) - e*rho*sin(theta - nfp*zeta)
    """
    rho = np.linspace(0.0, 1.0, n_rho)
    mode_m = [0, 1, 1]
    mode_n = [0, 0, 1]
    rmnc = np.column_stack([np.full(n_rho, r0), a * rho, e * rho])
    zmns = np.column_stack([np.zeros(n_rho), a * rho, -e * rho])
    return rho, mode_m, mode_n, rmnc, zmns


def make_source(**kwargs):
    """Build a valid StellaratorSource, overriding defaults via kwargs."""
    rho, mode_m, mode_n, rmnc, zmns = rotating_ellipse_modes()
    params = dict(
        rho=rho,
        emission_density=(1.0 - rho**2),
        mode_m=mode_m,
        mode_n=mode_n,
        rmnc=rmnc,
        zmns=zmns,
        num_field_periods=5,
        energy=openmc.stats.muir(e0=14.08e6, m_rat=5.0, kt=2.0e4),
    )
    params.update(kwargs)
    return openmc.StellaratorSource(**params)


def make_model(src):
    sphere = openmc.Sphere(r=2000.0, boundary_type='vacuum')
    return openmc.Model(
        geometry=openmc.Geometry([openmc.Cell(region=-sphere)]),
        settings=openmc.Settings(
            particles=100, batches=1, run_mode='fixed source', source=src),
    )


def reference_moments(rho, mode_m, mode_n, rmnc, zmns, nfp, emission_density,
                      n_theta=128, n_zeta=128):
    """Compute moments of the exact density S(rho)*R*|tau| by quadrature.

    The Fourier coefficients here are linear in rho, so the piecewise-linear
    coefficient interpolation used by the implementation is exact and the
    continuous density can be used as a reference.
    """
    rq = np.linspace(0.0, 1.0, 201)[1:]  # exclude axis (zero measure)
    theta = np.linspace(0.0, 2 * np.pi, n_theta, endpoint=False)
    zeta = np.linspace(0.0, 2 * np.pi / nfp, n_zeta, endpoint=False)
    rg, tg, zg = np.meshgrid(rq, theta, zeta, indexing='ij')

    # Coefficients are c_k(rho) = slope_k * rho + const_k
    slopes_r = (rmnc[-1] - rmnc[0])  # over rho in [0, 1]
    consts_r = rmnc[0]
    slopes_z = (zmns[-1] - zmns[0])
    consts_z = zmns[0]

    R = np.zeros_like(rg)
    Rr = np.zeros_like(rg)
    Rt = np.zeros_like(rg)
    Z = np.zeros_like(rg)
    Zr = np.zeros_like(rg)
    Zt = np.zeros_like(rg)
    for k, (m, n) in enumerate(zip(mode_m, mode_n)):
        ang = m * tg - n * nfp * zg
        c, s = np.cos(ang), np.sin(ang)
        rc = consts_r[k] + slopes_r[k] * rg
        zs = consts_z[k] + slopes_z[k] * rg
        R += rc * c
        Rr += slopes_r[k] * c
        Rt += -m * rc * s
        Z += zs * s
        Zr += slopes_z[k] * s
        Zt += m * zs * c

    S = np.interp(rg, rho, emission_density)
    f = S * np.abs(R * (Rr * Zt - Rt * Zr))
    w = f / f.sum()
    return {
        'rho': (w * rg).sum(),
        'R': (w * R).sum(),
        'Z': (w * Z).sum(),
        'Z2': (w * Z**2).sum(),
    }


def test_stellarator_source_roundtrip():
    src = make_source(strength=2.0, time=openmc.stats.Uniform(0.0, 1e-6))

    elem = src.to_xml_element()
    assert elem.get('type') == 'stellarator'

    new = openmc.SourceBase.from_xml_element(elem)
    assert isinstance(new, openmc.StellaratorSource)
    assert new.num_field_periods == src.num_field_periods
    assert new.strength == src.strength
    np.testing.assert_allclose(new.rho, src.rho)
    np.testing.assert_allclose(new.emission_density, src.emission_density)
    np.testing.assert_array_equal(new.mode_m, src.mode_m)
    np.testing.assert_array_equal(new.mode_n, src.mode_n)
    np.testing.assert_allclose(new.rmnc, src.rmnc)
    np.testing.assert_allclose(new.zmns, src.zmns)
    assert new.rmns is None
    assert new.zmnc is None
    assert len(new.energy) == 1
    assert isinstance(new.time, openmc.stats.Uniform)
    assert new.time.a == src.time.a
    assert new.time.b == src.time.b


def test_stellarator_source_asymmetric_roundtrip():
    rho, mode_m, mode_n, rmnc, zmns = rotating_ellipse_modes()
    src = make_source(rmns=0.01 * rmnc, zmnc=0.01 * zmns)
    new = openmc.SourceBase.from_xml_element(src.to_xml_element())
    np.testing.assert_allclose(new.rmns, src.rmns)
    np.testing.assert_allclose(new.zmnc, src.zmnc)


def test_stellarator_source_multiple_energies():
    rho, mode_m, mode_n, rmnc, zmns = rotating_ellipse_modes(n_rho=5)
    energies = [openmc.stats.muir(e0=14.08e6, m_rat=5.0, kt=kt)
                for kt in (1.0e4, 1.5e4, 2.0e4, 2.5e4, 3.0e4)]
    src = make_source(rho=rho, emission_density=np.ones_like(rho),
                      rmnc=rmnc, zmns=zmns, energy=energies)
    assert len(src.energy) == len(rho)
    assert src.time is None

    new = openmc.SourceBase.from_xml_element(src.to_xml_element())
    assert len(new.energy) == len(rho)
    assert new.time is None


@pytest.mark.parametrize("kwargs, match", [
    (dict(rho=np.linspace(0.1, 1.0, 17)), "must start at 0"),
    (dict(rho=np.linspace(0.0, 0.9, 17)), "must end at 1"),
    (dict(emission_density=np.ones(5)), "same length as rho"),
    (dict(emission_density=np.zeros(17)), "must contain a positive value"),
    (dict(mode_m=[0, 1]), "same length as mode_m"),
    (dict(mode_m=[0, -1, 1]), "must be >= 0"),
    (dict(rmnc=np.ones((17, 2))), "must have shape"),
    (dict(rmns=np.ones((17, 3))), "must both be given"),
    (dict(energy=[openmc.stats.muir(14.08e6, 5.0, 2.0e4)] * 2),
     "Number of energy distributions"),
    (dict(num_field_periods=0), "num_field_periods"),
])
def test_stellarator_source_invalid(kwargs, match):
    with pytest.raises(ValueError, match=match):
        make_source(**kwargs)


def test_stellarator_source_axisymmetric_sampling(run_in_tmpdir):
    """Circular axisymmetric plasma has analytic moments.

    With S = 1 the joint density ~ (R0 + a*rho*cos(theta)) * a^2 * rho gives
    E[rho] = 2/3 and E[R] = R0 + a^2/(4*R0).
    """
    r0, a = 600.0, 150.0
    rho, mode_m, mode_n, rmnc, zmns = circular_tokamak_modes(r0, a)
    src = make_source(
        rho=rho, emission_density=np.ones_like(rho), mode_m=mode_m,
        mode_n=mode_n, rmnc=rmnc, zmns=zmns, num_field_periods=1,
        energy=openmc.stats.delta_function(14.07e6),
    )
    model = make_model(src)

    sites = model.sample_external_source(30_000)
    xyz = np.array([site.r for site in sites])
    major_r = np.hypot(xyz[:, 0], xyz[:, 1])
    minor_rho = np.hypot(major_r - r0, xyz[:, 2]) / a

    assert minor_rho.max() < 1.0 + 1e-9
    assert_sample_mean(minor_rho, 2.0 / 3.0)
    assert_sample_mean(major_r, r0 + a**2 / (4.0 * r0))


@pytest.mark.parametrize('orientation', [1, -1])
def test_stellarator_source_rotating_ellipse_sampling(run_in_tmpdir, orientation):
    """Compare sampled moments of a 3-D configuration against quadrature of
    the exact density S(rho) * R * |tau|."""
    r0, a, e, nfp = 600.0, 150.0, 50.0, 5
    rho, mode_m, mode_n, rmnc, zmns = rotating_ellipse_modes(
        r0, a, e, nfp=nfp)
    zmns *= orientation
    emission_density = 1.0 - rho**2
    src = make_source(
        rho=rho, emission_density=emission_density, mode_m=mode_m,
        mode_n=mode_n, rmnc=rmnc, zmns=zmns, num_field_periods=nfp,
        energy=openmc.stats.delta_function(14.07e6),
    )
    model = make_model(src)

    sites = model.sample_external_source(30_000)
    xyz = np.array([site.r for site in sites])
    major_r = np.hypot(xyz[:, 0], xyz[:, 1])
    z = xyz[:, 2]

    ref = reference_moments(
        rho, mode_m, mode_n, rmnc, zmns, nfp, emission_density)
    assert_sample_mean(major_r, ref['R'])
    assert_sample_mean(z, ref['Z'])
    assert_sample_mean(z**2, ref['Z2'])


@pytest.mark.parametrize('file_format', ['netcdf3', 'hdf5'])
def test_stellarator_source_from_vmec(run_in_tmpdir, file_format):
    """Round-trip a synthetic rotating-ellipse wout file through from_vmec."""
    import h5py
    from scipy.io import netcdf_file

    ns, nfp = 16, 5
    r0_m, a_m, e_m = 6.0, 1.5, 0.5  # VMEC files are in meters
    s = np.linspace(0.0, 1.0, ns)
    rho = np.sqrt(s)
    xm = np.array([0, 1, 1])
    xn = np.array([0, 0, nfp])  # xn includes the nfp factor
    rmnc = np.column_stack([np.full(ns, r0_m), a_m * rho, e_m * rho])
    zmns = np.column_stack([np.zeros(ns), a_m * rho, -e_m * rho])

    if file_format == 'netcdf3':
        with netcdf_file('wout_test.nc', 'w') as ds:
            ds.createDimension('radius', ns)
            ds.createDimension('mn_mode', len(xm))
            ds.createVariable('rmnc', 'd', ('radius', 'mn_mode'))[:] = rmnc
            ds.createVariable('zmns', 'd', ('radius', 'mn_mode'))[:] = zmns
            ds.createVariable('xm', 'd', ('mn_mode',))[:] = xm
            ds.createVariable('xn', 'd', ('mn_mode',))[:] = xn
            ds.createVariable('nfp', 'i', ()).data.fill(nfp)
            ds.createVariable('lasym__logical__', 'i', ()).data.fill(0)
    else:
        # Exercise the HDF5 fallback without constructing a full NetCDF4 file.
        with h5py.File('wout_test.nc', 'w') as ds:
            for name, value in dict(rmnc=rmnc, zmns=zmns, xm=xm, xn=xn,
                                    nfp=nfp, lasym__logical__=0).items():
                ds[name] = value

    src = openmc.StellaratorSource.from_vmec(
        'wout_test.nc',
        emission_density=lambda r: 1.0 - r**2,
        energy=openmc.stats.delta_function(14.07e6),
    )
    assert src.num_field_periods == nfp
    np.testing.assert_allclose(src.rho, rho)
    np.testing.assert_array_equal(src.mode_m, xm)
    np.testing.assert_array_equal(src.mode_n, [0, 0, 1])
    np.testing.assert_allclose(src.rmnc, 100.0 * rmnc)  # m to cm
    np.testing.assert_allclose(src.zmns, 100.0 * zmns)
    np.testing.assert_allclose(src.emission_density, 1.0 - rho**2)
    assert src.rmns is None


@pytest.fixture(scope='module', params=[0, 1], ids=['axisymmetric', 'stellarator'])
def desc_equilibrium(request):
    pytest.importorskip('desc')
    from desc.equilibrium import Equilibrium

    # Circular torus plus a nontrivial l=3 radial contribution.
    # Replace the initial geometry below, so skip DESC's initial-guess checks.
    eq = Equilibrium(L=3, M=2, N=request.param, NFP=3, sym=False,
                     check_orientation=False, ensure_nested=False)
    r = np.zeros(eq.R_basis.num_modes)
    z = np.zeros(eq.Z_basis.num_modes)
    r[eq.R_basis.get_idx(0, 0, 0)] = 6.0
    r[eq.R_basis.get_idx(1, 1, 0)] = 1.5
    z[eq.Z_basis.get_idx(1, -1, 0)] = 1.5
    r[eq.R_basis.get_idx(3, -1, 0)] = 0.02
    z[eq.Z_basis.get_idx(2, 0, 0)] = 0.03
    if eq.N:
        # Exercise all product-basis sign combinations and m=0 harmonics.
        for m, n in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            r[eq.R_basis.get_idx(3, m, n)] = 0.02
            z[eq.Z_basis.get_idx(3, m, n)] = 0.03
        r[eq.R_basis.get_idx(0, 0, -1)] = 0.01
        z[eq.Z_basis.get_idx(0, 0, 1)] = 0.01
    eq.R_lmn, eq.Z_lmn = r, z
    return eq


@pytest.mark.parametrize('sym', [False, True])
def test_stellarator_source_from_desc(desc_equilibrium, sym):
    from desc.grid import Grid

    eq = desc_equilibrium
    if sym:
        # Keep the shared fixture asymmetric for the other conversion tests.
        eq = eq.copy()
        eq.change_resolution(sym=True)
    src = openmc.StellaratorSource.from_desc(
        eq, lambda r: 1-r**2, openmc.stats.delta_function(14.07e6), n_rho=7)
    assert src.num_field_periods == eq.NFP
    np.testing.assert_allclose(src.emission_density, 1-src.rho**2)
    assert (src.rmns is None) == sym

    # Compare physical surfaces against DESC's full basis evaluation, which
    # independently exercises radial evaluation, mode signs, and units.
    rng = np.random.default_rng(640)
    theta, zeta = rng.uniform(0, 2*np.pi, (2, 30))
    phase = src.mode_m[:, None]*theta - src.mode_n[:, None]*eq.NFP*zeta
    # Evaluate all surfaces together to avoid repeated DESC/JAX dispatch.
    nodes = np.column_stack([np.repeat(src.rho, theta.size),
                             np.tile(theta, src.rho.size),
                             np.tile(zeta, src.rho.size)])
    # Supplying a Grid preserves rho=0 without the default axis shift.
    grid = Grid(nodes)
    shape = (src.rho.size, theta.size)
    expected_r = np.asarray(eq.R_basis.evaluate(grid) @ eq.R_lmn).reshape(shape)*100
    expected_z = np.asarray(eq.Z_basis.evaluate(grid) @ eq.Z_lmn).reshape(shape)*100
    r = src.rmnc @ np.cos(phase)
    z = src.zmns @ np.sin(phase)
    if not sym:
        r += src.rmns @ np.sin(phase)
        z += src.zmnc @ np.cos(phase)
    np.testing.assert_allclose(r, expected_r, atol=1e-10, rtol=1e-12)
    np.testing.assert_allclose(z, expected_z, atol=1e-10, rtol=1e-12)


@pytest.mark.parametrize('family', [False, True])
def test_stellarator_source_desc_file(desc_equilibrium, tmp_path, family):
    from desc.equilibrium import EquilibriaFamily

    eq = desc_equilibrium
    obj = EquilibriaFamily(eq.copy(), eq) if family else eq
    if family:
        obj[0].R_lmn = obj[0].R_lmn * 2
    path = tmp_path/'equilibrium.h5'
    obj.save(path)
    args = dict(emission_density=lambda r: 1-r**2,
                energy=openmc.stats.delta_function(14.07e6), n_rho=5)
    reference = openmc.StellaratorSource.from_desc(eq, **args)
    for value in (obj, path, str(path)):
        source = openmc.StellaratorSource.from_desc(value, **args)
        assert source.num_field_periods == reference.num_field_periods
        for name in ('rho', 'mode_m', 'mode_n', 'rmnc', 'rmns', 'zmnc', 'zmns'):
            np.testing.assert_allclose(getattr(source, name),
                                       getattr(reference, name))


def test_stellarator_source_desc_invalid(desc_equilibrium, tmp_path):
    from desc.equilibrium import EquilibriaFamily

    args = dict(emission_density=lambda r: 1-r**2,
                energy=openmc.stats.delta_function(14.07e6))
    with pytest.raises(ValueError, match='empty'):
        openmc.StellaratorSource.from_desc(EquilibriaFamily(), **args)
    with pytest.raises(TypeError, match='eq'):
        openmc.StellaratorSource.from_desc(object(), **args)
    with pytest.raises(TypeError, match='n_rho'):
        openmc.StellaratorSource.from_desc(desc_equilibrium, n_rho=3.5, **args)
    with pytest.raises(ValueError, match='n_rho'):
        openmc.StellaratorSource.from_desc(desc_equilibrium, n_rho=1, **args)
    with pytest.raises(FileNotFoundError):
        openmc.StellaratorSource.from_desc(tmp_path/'missing.h5', **args)
    path = tmp_path/'surface.h5'
    desc_equilibrium.get_surface_at(rho=1.0).save(path)
    with pytest.raises(TypeError, match='eq'):
        openmc.StellaratorSource.from_desc(path, **args)


# ---------------------------------------------------------------------------
# Spin-polarized D-T emission (field direction from iota + lambda)
# ---------------------------------------------------------------------------

def make_polarized_source(**kwargs):
    """Return a circular-torus StellaratorSource with optional VMEC field data."""
    rho = np.linspace(0.0, 1.0, 5)
    # The third mode has zero geometry coefficient but carries a genuinely 3-D
    # lambda harmonic in polarized tests.
    mode_m = np.array([0, 1, 1])
    mode_n = np.array([0, 0, 1])
    rmnc = np.zeros((len(rho), len(mode_m)))
    zmns = np.zeros_like(rmnc)
    rmnc[:, 0] = 600.0
    rmnc[:, 1] = 150.0 * rho
    zmns[:, 1] = 150.0 * rho
    params = dict(
        rho=rho,
        emission_density=1.0 - 0.5 * rho**2,
        mode_m=mode_m,
        mode_n=mode_n,
        rmnc=rmnc,
        zmns=zmns,
        num_field_periods=2,
        energy=openmc.stats.delta_function(14.07e6),
    )
    params.update(kwargs)
    return openmc.StellaratorSource(**params)


def polarized_field_data():
    field_rho = np.array([0.1, 0.4, 0.7, 0.95])
    iota = 0.35 + 0.1 * field_rho
    lmns = np.zeros((len(field_rho), 3))
    lmns[:, 2] = 0.04 + 0.06 * field_rho
    return field_rho, iota, lmns


def test_stellarator_source_unpolarized_xml_unchanged():
    elem = make_polarized_source().to_xml_element()
    for name in ('polarization', 'field_rho', 'iota', 'lmns', 'lmnc'):
        assert elem.find(name) is None


def test_stellarator_source_polarization_roundtrip():
    field_rho, iota, lmns = polarized_field_data()
    src = make_polarized_source(
        polarization=(0.5, 0.3, 0.2), field_rho=field_rho, iota=iota,
        lmns=lmns)

    new = openmc.SourceBase.from_xml_element(src.to_xml_element())
    assert isinstance(new, openmc.StellaratorSource)
    assert new.polarization == pytest.approx(src.polarization)
    np.testing.assert_allclose(new.field_rho, field_rho)
    np.testing.assert_allclose(new.iota, iota)
    np.testing.assert_allclose(new.lmns, lmns)
    assert new.lmnc is None


def test_stellarator_source_from_vmec_uses_native_half_mesh(monkeypatch):
    """Classic wout lambda/iota row zero is padding, not an axis value."""
    ns = 5
    rho = np.sqrt(np.linspace(0.0, 1.0, ns))
    rmnc = np.zeros((ns, 2))
    zmns = np.zeros_like(rmnc)
    rmnc[:, 0] = 6.0
    rmnc[:, 1] = 1.5 * rho
    zmns[:, 1] = 1.5 * rho
    lmns = np.zeros((ns, 2))
    lmns[1:, 1] = [0.01, 0.02, 0.03, 0.04]
    data = dict(
        rmnc=rmnc, zmns=zmns, xm=np.array([0, 1]),
        xn=np.array([0, 0]), nfp=np.asarray(1),
        lasym__logical__=np.asarray(0), lmns=lmns,
        iotas=np.array([0.0, 0.31, 0.32, 0.33, 0.34]),
        iotaf=np.linspace(0.3, 0.35, ns),
    )
    monkeypatch.setattr(openmc.StellaratorSource, '_read_wout',
                        staticmethod(lambda path: data))

    src = openmc.StellaratorSource.from_vmec(
        'unused.nc', emission_density=np.ones(ns),
        energy=openmc.stats.delta_function(14.07e6),
        polarization=(1, 0, 0))
    expected_rho = np.sqrt((np.arange(1, ns) - 0.5) / (ns - 1))
    np.testing.assert_allclose(src.field_rho, expected_rho)
    np.testing.assert_allclose(src.iota, data['iotas'][1:])
    np.testing.assert_allclose(src.lmns, lmns[1:])


@pytest.mark.parametrize("kwargs, match", [
    (dict(polarization=(1, 0, 0)), "requires field_rho"),
    (dict(field_rho=[0.1, 0.9], iota=[0.4, 0.4],
          lmns=np.zeros((2, 3))), "only used with polarization"),
    (dict(polarization=(1, 0, 0), field_rho=[0.1, 0.9],
          iota=[0.4], lmns=np.zeros((2, 3))),
     "iota and field_rho"),
    (dict(polarization=(1, 0, 0), field_rho=[0.1, 0.9],
          iota=[0.4, 0.4], lmns=np.zeros((2, 2))), "lmns must have shape"),
])
def test_stellarator_source_invalid_field_data(kwargs, match):
    with pytest.raises(ValueError, match=match):
        make_polarized_source(**kwargs)


def _interpolate_clamped(x, xp, values):
    """Linear interpolation with the same endpoint clamping as the C++ path."""
    return np.interp(x, xp, values)


@pytest.mark.parametrize(("polarization", "expected_cos2"), [
    ((1.0, 0.0, 0.0), 1.0 / 5.0),
    ((1.0, 1.0, 1.0), 1.0 / 3.0),
    ((0.0, 1.0, 0.0), 7.0 / 15.0),
])
def test_stellarator_source_polarized_direction_sampling(
    run_in_tmpdir, polarization, expected_cos2
):
    """Check the compiled SPF sampler about an independently rebuilt B-hat."""
    R0, a, nfp = 600.0, 150.0, 2
    field_rho, iota_grid, lmns = polarized_field_data()
    src = make_polarized_source(
        polarization=polarization, field_rho=field_rho, iota=iota_grid,
        lmns=lmns)

    sphere = openmc.Sphere(r=2000.0, boundary_type='vacuum')
    model = openmc.Model(
        geometry=openmc.Geometry([openmc.Cell(region=-sphere)]),
        settings=openmc.Settings(
            particles=100, batches=1, run_mode='fixed source', source=src),
    )
    sites = model.sample_external_source(30_000)
    xyz = np.array([site.r for site in sites])
    directions = np.array([site.u for site in sites])

    zeta = np.arctan2(xyz[:, 1], xyz[:, 0])
    major_r = np.hypot(xyz[:, 0], xyz[:, 1])
    theta = np.arctan2(xyz[:, 2], major_r - R0)
    rho = np.hypot(major_r - R0, xyz[:, 2]) / a

    iota = _interpolate_clamped(rho, field_rho, iota_grid)
    lambda_coeff = _interpolate_clamped(rho, field_rho, lmns[:, 2])
    phase = theta - nfp * zeta
    lambda_theta = lambda_coeff * np.cos(phase)
    lambda_zeta = -nfp * lambda_coeff * np.cos(phase)
    b_theta = iota - lambda_zeta
    b_zeta = 1.0 + lambda_theta

    # Circular-torus coordinate basis, evaluated independently of the C++
    # Fourier derivative implementation.
    e_r = np.column_stack((np.cos(zeta), np.sin(zeta), np.zeros_like(zeta)))
    e_phi = np.column_stack((-np.sin(zeta), np.cos(zeta), np.zeros_like(zeta)))
    e_z = np.column_stack((np.zeros_like(zeta), np.zeros_like(zeta),
                           np.ones_like(zeta)))
    e_theta = (-a * rho * np.sin(theta))[:, None] * e_r \
        + (a * rho * np.cos(theta))[:, None] * e_z
    e_zeta = major_r[:, None] * e_phi
    field = b_theta[:, None] * e_theta + b_zeta[:, None] * e_zeta
    field /= np.linalg.norm(field, axis=1)[:, None]

    cos2 = np.einsum('ij,ij->i', directions, field)**2
    assert_sample_mean(cos2, expected_cos2)
    np.testing.assert_allclose(np.linalg.norm(directions, axis=1), 1.0,
                               rtol=0.0, atol=2e-14)
