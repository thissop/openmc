#!/usr/bin/env python
"""Equilibrium-field b_hat provider -- the G3 fix that admits reversal devices.

The coil (Biot-Savart) field carries each device's real B.n residual on the LCFS, so
the field_audit HARD-GATE excludes the worst-fit devices -- which are exactly the
high-iota QH configs with toroidal-sense reversal (e.g. device 1190023), the only ones
where the first-moment C decouples from the parity-correct nematic order S_phi
(docs/THEORY.md, SCOPE_AND_CAVEATS.md). Taking b_hat from the VMEC/Boozer EQUILIBRIUM
field instead -- a flux function, tangent to the LCFS by construction -- removes the
B.n residual so those devices pass the audit and can adjudicate C-vs-tensor.

EquilibriumField wraps the SAME portable field map (spf_fieldmap_v1) the compiled
source reads, produced offline by quasr_equilibrium_field.py (DESC; login node), and
exposes the CoilField interface (.B, .bhat) so field_audit / coherence_metrics /
biotsavart_field.write_fieldmap consume it unchanged. The map stores only the unit
DIRECTION (all the .so needs); a solved MHD equilibrium is divergence-free by
construction, so the finite-difference div-B gate is not meaningful here and is skipped
(div_free_by_construction) -- the boundary B.n gate DOES apply and is the real check.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

import fieldmap


class EquilibriumField:
    div_free_by_construction = True    # solved equilibrium => skip the FD div-B gate

    def __init__(self, stem, meta=None):
        stem = str(stem)
        if not (Path(stem + ".meta").exists() and Path(stem + ".bin").exists()):
            raise FileNotFoundError(
                f"equilibrium field map {stem}.(meta|bin) not found. Pre-generate it on "
                "a DESC machine (login node / desc venv):\n"
                f"    python spf_prototype/python/quasr_equilibrium_field.py <ID> "
                f"{Path(stem).name}\n"
                "field_source='equilibrium' requires this map (see docs/THEORY.md, G3).")
        self._fm = fieldmap.FieldMapField(stem)
        self.meta = meta or {}
        self.stem = stem

    def bhat(self, points):
        """(Q,3) unit vectors, interpolated IDENTICALLY to the compiled source (looping
        the parity-checked scalar FieldMapField.bhat, so C matches the .so's field)."""
        pts = np.asarray(points, float)
        out = np.empty((len(pts), 3), float)
        for i, r in enumerate(pts):
            out[i] = self._fm.bhat((float(r[0]), float(r[1]), float(r[2])))
        return out

    def B(self, points):
        """Direction only (the map stores unit b_hat). |B| magnitude is not carried:
        the audit's boundary gate normalizes by |B| (so direction suffices) and the div
        gate is skipped (div_free_by_construction)."""
        return self.bhat(points)
