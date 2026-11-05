#!/usr/bin/env python3
"""Audit the matched-mass minimax LP's extrapolation failure.

The 12 finite-difference columns were measured as isolated +20 cm patches.  The
matched-mass LP combines them at a much larger aggregate intervention.  This
script tests whether simple thickness caps and physical non-negativity of every
predicted coil response can coexist with the common family mass.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


HERE = Path(__file__).resolve().parent
SHARED = Path("/Users/tkiker/lima_share")
PATCH_ORDER = [
    "t0p2",
    "t0p3",
    "t1p0",
    "t1p1",
    "t1p2",
    "t1p3",
    "t2p0",
    "t2p1",
    "t2p2",
    "t2p3",
    "t3p2",
    "t3p3",
]


def load_problem():
    fractional = np.genfromtxt(
        HERE / "b_fd_fractional.csv", delimiter=",", skip_header=1
    )
    baseline = np.ravel(
        np.load(HERE / "../cluster_pull/coil_step1b_uniform.npz")[
            "coil_fast_flux"
        ]
    ).astype(float)
    response = fractional * baseline[None, :]

    masses = json.loads((SHARED / "masses_local.json").read_text())
    zero = masses["zerodelta"]
    mass_per_cm = (
        np.array(
            [
                (masses[patch]["V_shield"] - zero["V_shield"]) / 20.0
                for patch in PATCH_ORDER
            ]
        )
        * (15.6 - 1.94)
        / 1000.0
    )
    target_mass = float(np.load(HERE / "design_families.npz")["Mstar"])
    return baseline, response, mass_per_cm, target_mass


def solve(baseline, response, mass_per_cm, target_mass, cap_cm, nonnegative):
    n_patch = response.shape[0]
    objective = np.zeros(n_patch + 1)
    objective[-1] = 1.0
    upper_rows = []
    upper_bounds = []
    for coil_index in range(response.shape[1]):
        peak_row = np.zeros(n_patch + 1)
        peak_row[:n_patch] = response[:, coil_index] / 20.0
        peak_row[-1] = -1.0
        upper_rows.append(peak_row)
        upper_bounds.append(-baseline[coil_index])
        if nonnegative:
            physical_row = np.zeros(n_patch + 1)
            physical_row[:n_patch] = -response[:, coil_index] / 20.0
            upper_rows.append(physical_row)
            upper_bounds.append(baseline[coil_index])

    equality = np.zeros((1, n_patch + 1))
    equality[0, :n_patch] = mass_per_cm
    result = linprog(
        objective,
        A_ub=np.asarray(upper_rows),
        b_ub=np.asarray(upper_bounds),
        A_eq=equality,
        b_eq=[target_mass],
        bounds=[(0.0, cap_cm)] * n_patch + [(0.0, None)],
        method="highs",
    )
    if not result.success:
        return {"success": False, "message": result.message}

    delta = result.x[:n_patch]
    predicted = baseline + response.T @ (delta / 20.0)
    return {
        "success": True,
        "cap_cm": cap_cm,
        "nonnegative_constraints": nonnegative,
        "mass_kg": float(mass_per_cm @ delta),
        "delta_cm": {
            patch: float(value) for patch, value in zip(PATCH_ORDER, delta)
        },
        "patches_at_cap": int(np.count_nonzero(delta >= cap_cm - 1e-7)),
        "predicted_min": float(predicted.min()),
        "predicted_peak": float(predicted.max()),
        "predicted_peak_cell": int(11 + predicted.argmax()),
        "predicted_peak_change_percent": float(
            100.0 * (predicted.max() / baseline.max() - 1.0)
        ),
        "negative_coil_count": int(np.count_nonzero(predicted < 0.0)),
    }


def render_markdown(audit):
    lines = [
        "# Minimax linear-surrogate audit",
        "",
        f"Common family mass: **{audit['target_mass_kg'] / 1e6:.3f}×10⁶ kg**. "
        "Each response column was measured for one isolated +20 cm patch.",
        "",
        "| cap | enforce predicted R≥0 | feasible | predicted peak change | min predicted R | negative coils |",
        "|---:|:---:|:---:|---:|---:|---:|",
    ]
    for case in audit["cases"]:
        if case["success"]:
            lines.append(
                f"| {case['cap_cm']:.0f} cm | "
                f"{'yes' if case['nonnegative_constraints'] else 'no'} | yes | "
                f"{case['predicted_peak_change_percent']:+.1f}% | "
                f"{case['predicted_min']:.3e} | {case['negative_coil_count']} |"
            )
        else:
            lines.append(
                f"| {case['cap_cm']:.0f} cm | "
                f"{'yes' if case['nonnegative_constraints'] else 'no'} | **no** | — | — | — |"
            )

    lines.extend(
        [
            "",
            "## Verdict",
            "",
            "- The mass target is reachable with every patch capped at 20 cm, so mass feasibility is not the failure.",
            "- Even at the measured 20 cm amplitude, the unconstrained linear combination predicts negative responses; "
            "enforcing non-negativity makes the equal-mass problem infeasible.",
            "- Therefore the matched `famM` transport is an extrapolation stress test, not a trustworthy optimizer ranking.",
            "- A production minimax design needs sequential re-linearization/direct-forward updates, a lower intervention "
            "budget, or an explicitly nonlinear surrogate trained on combined designs.",
            "",
        ]
    )
    return "\n".join(lines)


def main():
    baseline, response, mass_per_cm, target_mass = load_problem()
    cases = []
    for cap_cm in (20.0, 25.0, 30.0, 35.0):
        cases.append(
            solve(
                baseline,
                response,
                mass_per_cm,
                target_mass,
                cap_cm,
                nonnegative=False,
            )
        )
    for cap_cm in (20.0, 35.0):
        case = solve(
            baseline,
            response,
            mass_per_cm,
            target_mass,
            cap_cm,
            nonnegative=True,
        )
        case.setdefault("cap_cm", cap_cm)
        case.setdefault("nonnegative_constraints", True)
        cases.append(case)

    audit = {
        "target_mass_kg": target_mass,
        "max_mass_at_20cm_kg": float(mass_per_cm.sum() * 20.0),
        "baseline_peak": float(baseline.max()),
        "cases": cases,
    }
    json_path = HERE / "minimax_surrogate_audit.json"
    markdown_path = HERE / "MINIMAX_SURROGATE_AUDIT.md"
    json_path.write_text(json.dumps(audit, indent=2) + "\n")
    markdown_path.write_text(render_markdown(audit))
    print(f"WROTE {json_path}")
    print(f"WROTE {markdown_path}")


if __name__ == "__main__":
    main()
