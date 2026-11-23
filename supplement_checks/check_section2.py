"""Check the Section 2 coefficients and three low-frequency scattering samples.

The saved bulk matrices are quadrature estimates. The uniform limit is
established by the proof in Supplement S2.
"""
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SECTION3 = HERE.parent / "section3"
FREQUENCIES = (0.2, 0.05, 0.0125)


def load_junction():
    source = SECTION3 / "independent/scripts/junction.py"
    spec = spec_from_file_location("section3_junction", source)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {source}")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_acoustic_tensor():
    """Compare all coefficients of the quadratic bulk symbol with Eq. (S11).

    Work in units K = m_b = a = 1. Coordinate directions and pair sums
    determine all six coefficients of each homogeneous quadratic polynomial.
    Each matrix entry is computed using exact rational arithmetic.
    """
    axes = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    diagonals = [tuple(x + sign*y for x, y in zip(u, v))
                 for u, v in combinations(axes, 2) for sign in (-1, 1)]
    samples = axes + [tuple(x + y for x, y in zip(u, v))
                      for u, v in combinations(axes, 2)]
    for k in samples:
        norm2 = sum(x*x for x in k)
        for i in range(3):
            for j in range(3):
                entry = sum(Fraction(d[i]*d[j] * sum(x*y for x, y in zip(k, d))**2,
                                     sum(x*x for x in d)) for d in axes + diagonals)
                expected = (norm2 if i == j else 0) + 2*k[i]*k[j]
                if entry != expected:
                    raise AssertionError(f"Acoustic tensor at {k}, entry {i},{j}: {entry}")
    print("All quadratic coefficients of the nine-family bulk symbol agree.")


def relative_operator_errors(matrices, omega, limits):
    return [
        np.linalg.norm(matrix / omega - limit, ord=2) / np.linalg.norm(limit, ord=2)
        for matrix, limit in zip(matrices, limits)
    ]


def main():
    check_acoustic_tensor()
    junction = load_junction()
    # Archive units: K = m_b = a = m_q = 1 and K_q = 1.69.
    chain_stiffness = 1.69
    impedance = np.sqrt(chain_stiffness)
    speeds = np.array([1.0, np.sqrt(3.0)])
    # Angular weights of the transverse pair and longitudinal mode.
    mu, weights = np.polynomial.legendre.leggauss(4)
    angles = 2 * np.pi * np.array([weights @ (1 - mu**2), weights @ mu**2])
    np.testing.assert_allclose(angles, [8*np.pi/3, 4*np.pi/3])
    b = np.pi * angles / ((2*np.pi)**3 * 2 * speeds**3)
    coefficients = 4 * impedance * b
    np.testing.assert_allclose(coefficients, [2*impedance/(3*np.pi),
                                               impedance/(3*np.pi*speeds[1]**3)])

    ramp = np.arange(1.0, 7.0) / 7
    profiles = [ramp, ramp[[0, 1, 4, 3, 2, 5]], *list(3.0 * np.eye(6))]
    # In the archive's orthonormal axial/radial basis this vector gives the
    # physical axial resultant of a corner-force vector on six segments.
    p = np.tile([2.0, 0.0], 6)
    reduced_limits = b[:, None, None] * np.outer(p, p)
    # In Cartesian coordinates the vector resultant includes all three axes.
    sites, _, _, _ = junction.assemble(ramp)
    resultant = np.tile(np.eye(3), (1, len(sites)))
    cartesian_limits = b[:, None, None] * (resultant.T @ resultant)
    largest_force_residual = 0.0
    power_errors, reduced_errors, cartesian_errors = [], [], []
    for omega in FREQUENCIES:
        file = HERE / "data" / f"green_w{omega:g}_N512.npz"
        errors = []
        with np.load(file, allow_pickle=False) as data:
            re_green, transverse, longitudinal = data["ReG"], data["AT"], data["AL"]
            matrices = [junction.site_matrix(sites, band) for band in (transverse, longitudinal)]
            cartesian_errors.append(relative_operator_errors(matrices, omega, cartesian_limits))
            for profile in profiles:
                state = junction.solve(profile, omega, re_green, transverse, longitudinal)
                k, q = state["k"], state["q"]
                gamma = 2 * chain_stiffness * np.sin(k)
                # Unit incident displacement. Internal chain forces cancel.
                end_force = (-1j * gamma
                             + chain_stiffness * np.expm1(1j*k) * (q[0] + q[-1])
                             + omega**2 * q.sum())
                residual = abs(state["Fx"] - end_force) / gamma
                if not np.isfinite(residual):
                    raise AssertionError(f"Nonfinite endpoint-force residual at frequency {omega}")
                largest_force_residual = max(largest_force_residual, residual)
                fractions = np.array([state["TqT"], state["TqL"]])
                errors.append(np.abs(fractions / (coefficients * omega**2) - 1))
        power_errors.append(np.max(errors, axis=0))
        file = SECTION3 / "paired/low_frequency/upstream/refined" / f"kernel_{omega:.8f}.npz"
        with np.load(file, allow_pickle=False) as data:
            spectral = np.stack([data["Q0"] + data["Q1"], data["Q2"]])
        reduced_errors.append(relative_operator_errors(spectral, omega, reduced_limits))

    if not largest_force_residual < 1e-11:
        raise AssertionError(f"Endpoint and spring force sums differ: {largest_force_residual:g}")
    checks = (("power fractions", power_errors),
              ("reduced bulk matrices", reduced_errors),
              ("Cartesian bulk matrices", cartesian_errors))
    for name, errors in checks:
        values = np.array(errors)
        # Tolerances below describe just these archived examples.
        if not (np.all(np.diff(values, axis=0) < 0) and np.all(values[-1] < 0.01)):
            raise AssertionError(f"Unexpected relative errors for {name}: {values}")
        print(f"{name}: largest relative error in each group at the smallest frequency "
              f"(qT, qL) = {values[-1, 0]:.2g}, {values[-1, 1]:.2g}")
    print(f"Endpoint and spring force sums agree in all {len(profiles)*len(FREQUENCIES)} cases "
          f"(largest residual relative to the drive: {largest_force_residual:.2g}).")


if __name__ == "__main__":
    main()
