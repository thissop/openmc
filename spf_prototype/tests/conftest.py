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
FIELD_DRIVER = BUILDDIR / "spf_field_driver"

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


def _build_field_driver() -> Path:
    """spf_field_driver: spf_field.hpp only (OpenMC-free, no random_lcg)."""
    BUILDDIR.mkdir(parents=True, exist_ok=True)
    cmd = [
        "g++", "-O2", "-std=c++17", f"-I{SRCDIR}",
        str(SRCDIR / "spf_field_driver.cpp"),
        "-o", str(FIELD_DRIVER),
    ]
    subprocess.run(cmd, check=True)
    return FIELD_DRIVER


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


@pytest.fixture(scope="session")
def field_driver() -> Path:
    if not FIELD_DRIVER.exists():
        _build_field_driver()
    return FIELD_DRIVER


@pytest.fixture(scope="session")
def run_field_driver(field_driver):
    def _run(bmode, positions, **params) -> np.ndarray:
        args = [str(field_driver), bmode]
        if bmode == "constant":
            args += [repr(float(params["b"][0])), repr(float(params["b"][1])),
                     repr(float(params["b"][2]))]
        elif bmode == "angled":
            args += [repr(float(params["alpha"])), repr(float(params["beta"]))]
        elif bmode == "fieldmap":
            args += [str(params["path"])]
        stdin = "\n".join(f"{x!r} {y!r} {z!r}" for x, y, z in positions)
        out = subprocess.run(args, input=stdin, check=True,
                             capture_output=True, text=True)
        return np.array([[float(t) for t in line.split()]
                         for line in out.stdout.strip().splitlines()])
    return _run
