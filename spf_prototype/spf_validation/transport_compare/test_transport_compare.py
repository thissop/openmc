from pathlib import Path

import numpy as np
import pytest

from analyze_transport_compare import independent_difference
from run_transport_compare import cells_from_npz, parse_cell_spec


def test_parse_cell_spec_accepts_ranges_and_ids():
    assert parse_cell_spec("11-13,20,18") == [11, 12, 13, 18, 20]


@pytest.mark.parametrize("spec", ["", "13-11", "11,11"])
def test_parse_cell_spec_rejects_invalid_input(spec):
    with pytest.raises(ValueError):
        parse_cell_spec(spec)


def test_cells_from_npz_recognizes_known_schema(tmp_path: Path):
    path = tmp_path / "cells.npz"
    np.savez(path, magnet_cells=np.array([14, 12, 13]))
    assert cells_from_npz(path) == [12, 13, 14]


def test_independent_difference_uses_quadrature_standard_error():
    delta, standard_error, z = independent_difference(
        np.array([10.0]), np.array([3.0]),
        np.array([14.0]), np.array([4.0]))
    np.testing.assert_allclose(delta, [4.0])
    np.testing.assert_allclose(standard_error, [5.0])
    np.testing.assert_allclose(z, [0.8])
