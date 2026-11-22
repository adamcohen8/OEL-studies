"""Stage 1 checks; numerical references integrate the volume, not the faces."""

import itertools
import unittest
from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss

from cube_bridge import cube_accel


@lru_cache(maxsize=None)
def rule(n):
    x, w = leggauss(n)
    return (x + 1.0) / 2.0, w / 2.0


def exterior_reference(q, n):
    """Direct tensor quadrature of the unit-cube volume kernel."""
    x, w = rule(n)
    s = np.stack(np.meshgrid(x - 0.5, x - 0.5, x - 0.5, indexing='ij'), axis=-1)
    r = s - np.asarray(q)
    r2 = np.sum(r * r, axis=-1)
    weights = w[:, None, None] * w[None, :, None] * w[None, None, :]
    return np.sum(r * (weights / (r2 * np.sqrt(r2)))[..., None], axis=(0, 1, 2))


def closed_cube_reference(q, n):
    """Volume quadrature at interior or exact boundary points.

    Split the volume into boxes meeting at q, dropping empty boxes. In each
    box, divide into three regions according to which normalized coordinate
    is largest. In region k the displacement is t*v, with v[k] the signed
    extent and the remaining v components given by the other signed extents
    times two independent variables in [0,1]. The Jacobian contains t**2,
    which cancels the kernel's t**-2. Integrate t out exactly, then apply
    two-dimensional Gaussian quadrature to the remaining smooth integrand.
    """
    q = np.asarray(q, dtype=float)
    x, w = rule(n)
    a, b = np.meshgrid(x, x, indexing='ij')
    weights = w[:, None] * w[None, :]
    result = np.zeros(3)
    for ends in itertools.product((-0.5, 0.5), repeat=3):
        ext = np.asarray(ends) - q
        jac = abs(np.prod(ext))
        if jac == 0.0:
            continue
        for k in range(3):
            i, j = (k + 1) % 3, (k + 2) % 3
            v = np.empty((n, n, 3))
            v[..., k] = ext[k]
            v[..., i] = ext[i] * a
            v[..., j] = ext[j] * b
            r2 = np.sum(v * v, axis=-1)
            result += np.sum(v * (weights * jac / (r2 * np.sqrt(r2)))[..., None], axis=(0, 1))
    return result


class CubeBridgeTests(unittest.TestCase):
    def check_reference(self, point, reference):
        coarse = reference(point, 32)
        fine = reference(point, 64)
        actual = cube_accel(point, 1.0, 1.0)
        print(f'\npoint={point} bridge={actual.tolist()} reference={fine.tolist()} '
              f'ref_delta={np.max(np.abs(coarse - fine)):.3g} '
              f'bridge_delta={np.max(np.abs(actual - fine)):.3g}', flush=True)
        np.testing.assert_allclose(coarse, fine, rtol=0.0, atol=2e-13)
        np.testing.assert_allclose(actual, fine, rtol=0.0, atol=4e-13)
        if any(point):
            self.assertLess(np.dot(point, actual), 0.0)

    def test_exterior_volume_reference(self):
        for p in [(0.72, 0.12, -0.22), (1.25, -0.9, 0.65), (4.0, 0.0, 0.0)]:
            with self.subTest(point=p):
                self.check_reference(p, exterior_reference)

    def test_interior_and_boundary_volume_reference(self):
        for p in [(0.0, 0.0, 0.0), (0.13, -0.21, 0.07), (0.5, 0.0, 0.0),
                  (-0.5, 0.18, -0.11), (0.5, -0.5, 0.19), (0.5, -0.5, 0.5)]:
            with self.subTest(point=p):
                self.check_reference(p, closed_cube_reference)

    def test_face_center_one_sided_limits(self):
        ref = closed_cube_reference((0.5, 0.0, 0.0), 64)
        on = cube_accel((0.5, 0.0, 0.0), 1.0, 1.0)
        np.testing.assert_allclose(on, ref, rtol=0.0, atol=4e-13)
        for direction in (-1, 1):
            errors = []
            for eps in (1e-3, 1e-6, 1e-9):
                actual = cube_accel((0.5 + direction * eps, 0.0, 0.0), 1.0, 1.0)
                np.testing.assert_array_equal(actual[1:], np.zeros(2))
                errors.append(abs(actual[0] - ref[0]))
            print(f'\nface_center direction={direction:+d} errors={errors}', flush=True)
            self.assertLess(errors[1], errors[0] / 100.0)
            self.assertLess(errors[2], errors[1] / 100.0)
            self.assertLess(errors[2], 2e-8)

    def test_symmetry_and_scaling(self):
        p = np.array([0.23, 0.5, -0.41])
        unit = cube_accel(p, 1.0, 1.0)
        perm = [2, 0, 1]
        signs = np.array([-1, 1, -1])
        np.testing.assert_allclose(cube_accel(signs * p[perm], 1.0, 1.0),
                                   signs * unit[perm], rtol=0.0, atol=5e-15)
        np.testing.assert_allclose(cube_accel(8.0 * p, 8.0, 128.0),
                                   2.0 * unit, rtol=0.0, atol=5e-15)
        np.testing.assert_array_equal(cube_accel((0, 0, 0), 2.0, 3.0), np.zeros(3))

    def test_invalid_inputs(self):
        for p, side, mu in [((0, 0), 1, 1), ((0, 0, np.nan), 1, 1),
                            ((0, 0, 0), 0, 1), ((0, 0, 0), 1, -1),
                            ((0, 0, 0), 1, np.inf), ((4.01, 0, 0), 1, 1)]:
            with self.subTest(position=p, side=side, mu=mu):
                with self.assertRaises(ValueError):
                    cube_accel(p, side, mu)


if __name__ == '__main__':
    unittest.main()
