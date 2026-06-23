"""Shared fixtures for tier-1 tests: put spf_prototype/python on sys.path, build
the standalone C++ driver once per session, and expose a run_driver() helper."""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]          # .../openmc
PYDIR = REPO / "spf_prototype" / "python"
SRCDIR = REPO / "spf_prototype" / "src"
BUILDDIR = REPO / "spf_prototype" / "build"
DRIVER = BUILDDIR / "spf_driver"

sys.path.insert(0, str(PYDIR))


def _build_driver() -> Path:
    BUILDDIR.mkdir(parents=True, exist_ok=True)
    cmd = [
        "g++", "-O2", "-std=c++17",
        f"-I{REPO / 'include'}", f"-I{SRCDIR}",
        str(SRCDIR / "standalone_driver.cpp"),
        str(REPO / "src" / "random_lcg.cpp"),
        "-o", str(DRIVER),
    ]
    subprocess.run(cmd, check=True)
    return DRIVER


@pytest.fixture(scope="session")
def driver() -> Path:
    if not DRIVER.exists():
        _build_driver()
    return DRIVER


@pytest.fixture(scope="session")
def run_driver(driver):
    def _run(a, b, c, B, N, seed, mode="invcdf") -> np.ndarray:
        out = subprocess.run(
            [str(driver), repr(float(a)), repr(float(b)), repr(float(c)),
             repr(float(B[0])), repr(float(B[1])), repr(float(B[2])),
             str(int(N)), str(int(seed)), mode],
            check=True, capture_output=True, text=True,
        )
        # columns: cz ux uy uz
        return np.array([[float(t) for t in line.split()]
                         for line in out.stdout.strip().splitlines()])
    return _run
