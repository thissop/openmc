"""Python API for the native ``StellaratorSource`` (mirrors ``TokamakSource``).

This is the Python-side counterpart of ``native_spf/stellarator_source.{h,cpp}``:
a :class:`openmc.SourceBase` subclass that serializes a stellarator plasma
neutron source to/from XML. The C++ core reads that XML plus a portable
``spf_fluxmap_v1`` file (``<stem>.meta`` + ``<stem>.bin``, written by
``python/quasr_fluxmap.py:write_fluxmap_bin``) and samples birth positions from
the flux-map CDF cascade and birth directions from the Schwartz P2 spin-polarized
mixture about the local field b-hat read from the same grid.

Staging note: this file is written as a standalone importable module so it can be
unit-tested and, when merged upstream, dropped into ``openmc/source.py`` (add a
``'stellarator'`` branch to ``SourceBase.from_xml_element`` -- see
``README_stellarator.md``). ``polarization=None`` (the default) writes no
polarization XML, so the C++ takes the isotropic draw (backward compatible).
"""
from __future__ import annotations

from numbers import Real
from typing import Any
import warnings

import lxml.etree as ET
import numpy as np

import openmc
import openmc.checkvalue as cv
from openmc.checkvalue import PathLike
from openmc.source import SourceBase
from openmc.stats.univariate import Univariate
from openmc.utility_funcs import input_path
from openmc._xml import get_text


def spin_fractions_to_abc(d_plus, d_zero, d_minus, t_plus, t_minus):
    """Schwartz 2025 Eq. 1: map deuteron/triton spin-projection fractions to the
    (a, b, c) collision-mode fractions.

        a = d_plus*t_plus + d_minus*t_minus
        b = d_zero
        c = d_plus*t_minus + d_minus*t_plus
    """
    a = d_plus * t_plus + d_minus * t_minus
    b = d_zero
    c = d_plus * t_minus + d_minus * t_plus
    return (a, b, c)


class StellaratorSource(SourceBase):
    """A source representing spin-polarized neutron emission from a stellarator.

    Birth positions are sampled from a real-equilibrium flux map (the
    ``spf_fluxmap_v1`` file) via a rejection-free marginal->conditional CDF
    cascade ``p(rho) p(zeta|rho) p(theta|rho,zeta)`` (so every birth weight is
    exactly 1). Birth directions are isotropic when ``polarization is None``, or
    else drawn from the Schwartz P2 mode-mixture about the local magnetic-field
    direction b-hat read from the SAME grid.

    Parameters
    ----------
    fluxmap : path-like
        Stem of the ``spf_fluxmap_v1`` file (``<stem>.meta`` + ``<stem>.bin``).
        A trailing ``.meta``/``.bin`` is accepted and stripped.
    polarization : None, sequence of float, or dict, optional
        Spin-polarization. ``None`` (default) => isotropic emission. A 3-sequence
        is the collision-mode fractions ``(a, b, c)`` (Schwartz Eq. 1). A dict of
        ``d_plus, d_zero, d_minus, t_plus, t_minus`` is mapped to ``(a, b, c)``
        via :func:`spin_fractions_to_abc`. Fractions must be >= 0; they are
        renormalized to sum 1 (a warning is issued if the input sum deviates).
    field_model : {'fluxmap', 'toroidal'}, optional
        How the C++ obtains b-hat when polarized. ``'fluxmap'`` (default) reads
        and normalizes the gridded field ``(BR, Bphi, BZ)``. ``'toroidal'``
        overrides it with the pure toroidal direction phi-hat (the axisymmetric
        limit used to cross-check against :class:`openmc.TokamakSource`).
    energy : openmc.stats.Univariate or sequence, optional
        Birth energy distribution. Defaults to monoenergetic 14.06 MeV. (The C++
        prototype uses a single distribution for all radii.)
    time : openmc.stats.Univariate, optional
        Time distribution. If ``None``, particles are born at ``t=0``.
    emission : tuple(sequence, sequence), optional
        Optional radial emission profile ``(rho_grid, density)`` overriding the
        default parabolic ``S(rho) = 1 - rho^2``. Only the shape matters.
    strength : float, optional
        Source strength (default 1.0).
    constraints : dict, optional
        See :class:`openmc.SourceBase`.

    Attributes
    ----------
    type : str
        Indicator of source type: ``'stellarator'``.
    """

    def __init__(
        self,
        fluxmap: PathLike,
        polarization: Any | None = None,
        field_model: str = 'fluxmap',
        energy: Univariate | None = None,
        time: Univariate | None = None,
        emission: tuple | None = None,
        strength: float = 1.0,
        constraints: dict[str, Any] | None = None,
    ):
        super().__init__(strength=strength, constraints=constraints)
        self.fluxmap = fluxmap
        self.field_model = field_model
        self.polarization = polarization
        self.energy = (energy if energy is not None
                       else openmc.stats.Discrete([14.06e6], [1.0]))
        self.time = time
        self.emission = emission

    @property
    def type(self) -> str:
        return 'stellarator'

    # -- fluxmap ------------------------------------------------------------
    @property
    def fluxmap(self) -> PathLike:
        return self._fluxmap

    @fluxmap.setter
    def fluxmap(self, value: PathLike):
        cv.check_type('fluxmap', value, PathLike)
        s = str(value)
        if s.endswith('.meta') or s.endswith('.bin'):
            s = s.rsplit('.', 1)[0]
        self._fluxmap = input_path(s)

    # -- field_model --------------------------------------------------------
    @property
    def field_model(self) -> str:
        return self._field_model

    @field_model.setter
    def field_model(self, value: str):
        cv.check_value('field_model', value, ('fluxmap', 'toroidal'))
        self._field_model = value

    # -- polarization -------------------------------------------------------
    @property
    def polarization(self):
        return self._polarization

    @polarization.setter
    def polarization(self, value):
        if value is None:
            self._polarization = None
            return
        if isinstance(value, dict):
            a, b, c = spin_fractions_to_abc(**value)
        else:
            seq = tuple(float(x) for x in value)
            if len(seq) != 3:
                raise ValueError(
                    'polarization must be a 3-tuple (a, b, c) or a spin-fraction '
                    'dict')
            a, b, c = seq
        if not (a >= 0 and b >= 0 and c >= 0):
            raise ValueError('polarization fractions a, b, c must each be >= 0')
        s = a + b + c
        if not (s > 0):
            raise ValueError('polarization a + b + c must be > 0')
        if abs(s - 1.0) > 1e-6:
            warnings.warn(
                f'polarization (a, b, c) sums to {s:.6g}; renormalizing to 1.')
        self._polarization = (a / s, b / s, c / s)

    # -- energy -------------------------------------------------------------
    @property
    def energy(self) -> list[Univariate]:
        return self._energy

    @energy.setter
    def energy(self, value):
        if isinstance(value, Univariate):
            self._energy = [value]
        else:
            cv.check_iterable_type('energy distributions', value, Univariate)
            self._energy = list(value)

    # -- time ---------------------------------------------------------------
    @property
    def time(self) -> Univariate | None:
        return self._time

    @time.setter
    def time(self, value):
        if value is not None:
            cv.check_type('time distribution', value, Univariate)
        self._time = value

    # -- emission -----------------------------------------------------------
    @property
    def emission(self):
        return self._emission

    @emission.setter
    def emission(self, value):
        if value is None:
            self._emission = None
            return
        rho, dens = value
        rho = np.asarray(rho, dtype=float)
        dens = np.asarray(dens, dtype=float)
        if rho.ndim != 1 or rho.size < 2 or rho.shape != dens.shape:
            raise ValueError(
                'emission must be (rho_grid, density) 1-D arrays of equal '
                'length >= 2')
        if not np.all(np.diff(rho) > 0):
            raise ValueError('emission rho grid must be strictly increasing')
        if np.any(dens < 0):
            raise ValueError('emission density must be >= 0')
        self._emission = (rho, dens)

    # -- XML ----------------------------------------------------------------
    def populate_xml_element(self, element):
        ET.SubElement(element, 'fluxmap').text = str(self.fluxmap)

        if self.emission is not None:
            rho, dens = self.emission
            ET.SubElement(element, 'emission_rho').text = \
                ' '.join(str(r) for r in rho)
            ET.SubElement(element, 'emission_density').text = \
                ' '.join(str(d) for d in dens)

        for dist in self.energy:
            element.append(dist.to_xml_element('energy'))

        if self.time is not None:
            element.append(self.time.to_xml_element('time'))

        # Polarization XML is emitted ONLY when polarized (absent => the C++
        # sees polarized_ = false and takes the unchanged isotropic draw).
        if self.polarization is not None:
            ET.SubElement(element, 'polarization').text = \
                ' '.join(str(x) for x in self.polarization)
            ET.SubElement(element, 'field_model').text = self.field_model

    @classmethod
    def from_xml_element(cls, elem: ET.Element, meshes=None) -> 'StellaratorSource':
        fluxmap = get_text(elem, 'fluxmap')

        pol_text = get_text(elem, 'polarization')
        polarization = (tuple(float(x) for x in pol_text.split())
                        if pol_text else None)
        field_model = get_text(elem, 'field_model') or 'fluxmap'

        energy = [Univariate.from_xml_element(e) for e in elem.findall('energy')]
        if len(energy) == 1:
            energy = energy[0]
        elif len(energy) == 0:
            energy = None

        time_elem = elem.find('time')
        time = (Univariate.from_xml_element(time_elem)
                if time_elem is not None else None)

        er = get_text(elem, 'emission_rho')
        ed = get_text(elem, 'emission_density')
        emission = None
        if er and ed:
            emission = ([float(x) for x in er.split()],
                        [float(x) for x in ed.split()])

        constraints = cls._get_constraints(elem)
        strength_text = get_text(elem, 'strength')
        strength = float(strength_text) if strength_text else 1.0

        return cls(
            fluxmap=fluxmap,
            polarization=polarization,
            field_model=field_model,
            energy=energy,
            time=time,
            emission=emission,
            strength=strength,
            constraints=constraints,
        )
