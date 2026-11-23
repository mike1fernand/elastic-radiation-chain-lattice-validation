"""Check one linearized contact segment and the finite-profile controls.

Run with Python 3: python supplement_checks/check_section1.py
All calculations use exact fractions. The checks concern equilibrium geometry
and the contact energy; they do not compute a scattering solution.
"""
from fractions import Fraction
from itertools import product
import unittest


SPACING = Fraction(5, 3)
STIFFNESS = Fraction(7, 5)
SIGNS = tuple(product((-1, 1), repeat=2))
LINE = ((0, 0, 0), (SPACING, 0, 0))
CORNERS = tuple((SPACING / 2, sy * SPACING / 2, sz * SPACING / 2)
                for sy, sz in SIGNS)
# State = two scalar line displacements, then four 3-vector corner displacements.
DOF = 2 + 3 * len(CORNERS)
ZERO = (Fraction(0),) * DOF


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def minus(a, b):
    return tuple(x - y for x, y in zip(a, b))


def corner_slice(c):
    return slice(2 + 3 * c, 5 + 3 * c)


def energy(state):
    """Sum the harmonic bond energies using equilibrium positions directly."""
    total = Fraction(0)
    for j, endpoint in enumerate(LINE):
        for c, corner in enumerate(CORNERS):
            bond = minus(corner, endpoint)
            relative = minus(state[corner_slice(c)], (state[j], 0, 0))
            extension_squared = dot(bond, relative) ** 2 / dot(bond, bond)
            total += STIFFNESS * extension_squared / 2
    return total


def forces(state):
    """Evaluate the force law using the signed components of the unit bonds."""
    result = list(ZERO)
    for c, (sy, sz) in enumerate(SIGNS):
        ux, uy, uz = state[corner_slice(c)]
        for j, sx in enumerate((1, -1)):
            # The unit vector is (sx, sy, sz)/sqrt(3). One 1/sqrt(3)
            # occurs in the extension and the other in the force direction.
            projection = sx * (ux - state[j]) + sy * uy + sz * uz
            strength = STIFFNESS * projection / 3
            result[j] += strength * sx
            for axis, sign in enumerate((sx, sy, sz)):
                result[2 + 3 * c + axis] -= strength * sign
    return tuple(result)


def shifted(state, index, delta):
    result = list(state)
    result[index] += delta
    return tuple(result)


class ContactChecks(unittest.TestCase):
    def test_distinct_bonds_and_lengths(self):
        bonds = [(line, corner) for line in LINE for corner in CORNERS]
        self.assertEqual(len(set(bonds)), 8)
        for line, corner in bonds:
            bond = minus(corner, line)
            self.assertEqual(dot(bond, bond), 3 * SPACING ** 2 / 4)

    def test_forces_as_energy_derivatives(self):
        # The force map is linear. Testing every basis input and every
        # derivative checks the complete 14-by-14 segment matrix.
        h = Fraction(1, 13)
        for basis in range(DOF):
            state = shifted(ZERO, basis, 1)
            force = forces(state)
            for k in range(DOF):
                derivative = (energy(shifted(state, k, h))
                              - energy(shifted(state, k, -h))) / (2 * h)
                self.assertEqual(force[k], -derivative, (basis, k))
            self.assertEqual(sum(force[:2]) + sum(force[2::3]), 0)

    def test_common_axial_translation(self):
        t = Fraction(3, 7)
        translation = (t, t) + (t, 0, 0) * len(CORNERS)
        self.assertEqual(energy(translation), 0)
        self.assertEqual(forces(translation), ZERO)

    def test_two_link_force_on_undisplaced_corner(self):
        # The two axial basis inputs span all motions of the segment endpoints.
        for left, right in ((Fraction(1), 0), (0, Fraction(1))):
            state = (left, right) + ZERO[2:]
            actual = forces(state)
            for c, (sy, sz) in enumerate(SIGNS):
                expected = tuple(STIFFNESS * v / 3 for v in
                                 (left + right, (left - right) * sy,
                                  (left - right) * sz))
                self.assertEqual(actual[corner_slice(c)], expected)

    def test_line_driven_radial_resultant(self):
        # Here corners are at equilibrium. Extension about a fixed midpoint
        # pulls each corner inward; the four radial vectors sum to zero.
        for delta in (Fraction(1, 17), Fraction(-1, 17)):
            state = (-delta / 2, delta / 2) + ZERO[2:]
            force = forces(state)
            for c, (sy, sz) in enumerate(SIGNS):
                self.assertEqual(force[corner_slice(c)],
                                 (0, -STIFFNESS * delta * sy / 3,
                                  -STIFFNESS * delta * sz / 3))
            self.assertEqual(sum(force[3::3]), 0)
            self.assertEqual(sum(force[4::3]), 0)

    def test_same_end_exchange(self):
        ramp = tuple(Fraction(x, 7) for x in (1, 2, 3, 4, 5, 6))
        swap = tuple(Fraction(x, 7) for x in (1, 2, 5, 4, 3, 6))
        self.assertEqual(sorted(ramp), sorted(swap))
        self.assertEqual([i for i, (a, b) in enumerate(zip(ramp, swap), 1)
                          if a != b], [3, 5])
        for profile in (ramp, swap):
            self.assertEqual(sum(profile), 3)
            self.assertEqual(8 * sum(profile), 24)
            self.assertEqual((profile[0], profile[-1]),
                             (Fraction(1, 7), Fraction(6, 7)))


if __name__ == "__main__":
    unittest.main(verbosity=2)
