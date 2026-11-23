"""Recompute low-frequency two-channel samples for the two nominal six-segment profiles."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
UPSTREAM = ROOT / "upstream"
FREQUENCIES = (0.2, 0.05, 0.0125)
CHANNELS = ("R", "Tline", "TqT", "TqL")
CT = 2 * 1.3 / (3 * np.pi)
CL = 1.3 / (3 * np.pi * np.sqrt(3)**3)


def verify_sources():
    manifest = json.loads((ROOT / "provenance.json").read_text())
    for name, record in manifest["files"].items():
        data = (UPSTREAM / name).read_bytes()
        git_blob = b"blob " + str(len(data)).encode() + b"\0" + data
        if (hashlib.sha256(data).hexdigest() != record["sha256"] or
                hashlib.sha1(git_blob).hexdigest() != record["git_blob"]):
            raise ValueError(f"Input hash mismatch for {name}")


def direct_response(omega, profile, real_green, branches):
    """Cross-check with an independently assembled full line-and-force matrix."""
    green = real_green + 1j * sum(branches)
    stiffness = 1.69
    n = len(profile)
    num_line, num_corner = n + 1, 2 * n
    total = num_line + num_corner
    hessian = np.zeros((total, total))
    for j, strength in enumerate(profile):
        for sy, sz in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
            corner = np.array([j + 0.5, sy / 2, sz / 2])
            for endpoint in (j, j + 1):
                direction = corner - np.array([endpoint, 0, 0])
                direction /= np.linalg.norm(direction)
                extension = np.zeros(total)
                extension[endpoint] = -direction[0]
                extension[num_line + 2 * j] = direction[0] / 2
                transverse_direction = np.array([0, sy, sz]) / (2 * np.sqrt(2))
                extension[num_line + 2 * j + 1] = direction @ transverse_direction
                hessian += strength * np.outer(extension, extension)
    upper = hessian[:num_line, :num_line]
    v = hessian[num_line:, :num_line]
    u = hessian[num_line:, num_line:]
    wave_number = 2 * np.arcsin(omega / (2 * np.sqrt(stiffness)))
    phase = np.exp(1j * wave_number)
    gamma = 2 * stiffness * phase.imag
    h = (omega**2 - 2 * stiffness) * np.eye(num_line, dtype=complex) - upper
    h += stiffness * (np.eye(num_line, k=1) + np.eye(num_line, k=-1))
    h[0, 0] += stiffness * phase
    h[-1, -1] += stiffness * phase
    matrix = np.block([[h, v.T @ green], [-v, np.eye(num_corner) + u @ green]])
    rhs = np.zeros(total)
    rhs[0] = np.sqrt(gamma)
    state = np.linalg.solve(matrix, rhs)
    force = state[num_line:]
    emitted = [float(2 * np.vdot(force, b @ force).real) for b in branches]
    return np.array([abs(1 - 1j * np.sqrt(gamma) * state[0])**2,
                     gamma * abs(state[num_line - 1])**2,
                     emitted[0] + emitted[1], emitted[2]])


def compute():
    verify_sources()
    from upstream.low_frequency import GREF, scatter

    ramp = np.arange(1., 7.) / 7
    profiles = {"reference_taper": GREF, "ramp": ramp, "swap": ramp[[0, 1, 4, 3, 2, 5]]}
    rows = []
    for level in ("refined", "validation"):
        folder = UPSTREAM / level
        baseline = json.loads((folder / "rows.json").read_text())
        for omega in FREQUENCIES:
            with np.load(folder / f"kernel_{omega:.8f}.npz", allow_pickle=False) as archive:
                real_green = archive["PV"]
                branches = [archive[f"Q{i}"] for i in range(3)]
            for name, profile in profiles.items():
                originals = scatter(omega, 1., real_green, branches, profile)
                powers = np.array([originals[key] for key in CHANNELS])
                direct = direct_response(omega, profile, real_green, branches)
                np.testing.assert_allclose(direct, powers, rtol=1e-10, atol=2e-14)
                if abs(sum(powers) - 1) > 1e-12 or min(powers) < 0:
                    raise ValueError("Power balance or positivity failure")
                ref_error = None
                if name == "reference_taper":
                    ref = next(r for r in baseline if r["omega"] == omega and r["s"] == 1)
                    expected = np.array([ref[key] for key in CHANNELS])
                    np.testing.assert_allclose(powers, expected, rtol=1e-10, atol=2e-14)
                    ref_error = float(np.max(abs(powers - expected)))
                rows.append(dict(kernel=level, omega=omega, profile=name,
                                 **dict(zip(CHANNELS, map(float, powers))),
                                 ratio_qT=float(powers[2] / (CT * omega**2)),
                                 ratio_qL=float(powers[3] / (CL * omega**2)),
                                 closure_error=float(sum(powers) - 1),
                                 second_solve_max_abs_error=float(np.max(abs(powers - direct))),
                                 taper_baseline_max_abs_error=ref_error))
    return {"CT": CT, "CL": CL, "rows": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="Print full-precision outputs and checks.")
    args = parser.parse_args()
    result = compute()
    if args.json:
        print(json.dumps(result, indent=2, allow_nan=False))
    else:
        print("kernel        omega    profile            T_qT/(CT*w^2)   T_qL/(CL*w^2)")
        for row in result["rows"]:
            print(f'{row["kernel"]:12} {row["omega"]:7.4f}  {row["profile"]:16} '
                  f'{row["ratio_qT"]:15.8f} {row["ratio_qL"]:15.8f}')
        print("Maximum absolute power-balance residual:",
              max(abs(r["closure_error"]) for r in result["rows"]))
        print("Maximum absolute discrepancy between the two line solves (any channel):",
              max(r["second_solve_max_abs_error"] for r in result["rows"]))
        print("Maximum absolute discrepancy with archived reference taper powers:",
              max(r["taper_baseline_max_abs_error"] or 0.0 for r in result["rows"]))
