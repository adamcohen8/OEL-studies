import unittest
from unittest.mock import patch

import numpy as np
from sim.dynamics.orbit.accelerations import accel_two_body

from dust_bridge import acceleration


class AccelerationTests(unittest.TestCase):
    # Inputs are in km, km/s, and km^3/s^2; expected accelerations are in km/s^2.
    def test_zero_beta_delegates_to_native_gravity(self):
        r = np.array([3.0, 4.0, 0.0])
        expected = accel_two_body(r, 125.0)
        with patch('dust_bridge.accel_two_body', wraps=accel_two_body) as native:
            actual = acceleration(r, [1, -2, 2], 125, 0, 100)
        native.assert_called_once()
        called_r, called_mu = native.call_args.args
        np.testing.assert_array_equal(called_r, r)
        self.assertEqual(called_mu, 125.0)
        np.testing.assert_array_equal(actual, expected)

    def test_stationary_direct_pressure(self):
        for beta, expected_x in [(0.25, -3.75), (1.0, 0.0), (1.25, 1.25)]:
            with self.subTest(beta=beta):
                actual = acceleration([5, 0, 0], [0, 0, 0], 125, beta, 100)
                np.testing.assert_allclose(actual, [expected_x, 0, 0], rtol=0, atol=1e-13)

    def test_tangential_velocity_has_no_radial_drag(self):
        actual = acceleration([5, 0, 0], [0, 2, 0], 125, 0.25, 100)
        np.testing.assert_allclose(actual, [-3.75, -0.025, 0], rtol=0, atol=1e-13)

    def test_radial_velocity_has_two_drag_contributions(self):
        # At beta=1 static forces cancel, isolating the signed drag.
        for vx, expected_x in [(2.0, -0.2), (-2.0, 0.2)]:
            with self.subTest(vx=vx):
                actual = acceleration([5, 0, 0], [vx, 0, 0], 125, 1, 100)
                np.testing.assert_allclose(actual, [expected_x, 0, 0], rtol=0, atol=1e-13)

    def test_mixed_velocity_off_axis(self):
        # Radial speed is -1 km/s; drag adds (-0.005, 0.035, -0.025).
        actual = acceleration([3, 4, 0], [1, -2, 2], 125, 0.25, 100)
        np.testing.assert_allclose(actual, [-2.255, -2.965, -0.025], rtol=0, atol=1e-13)

    def test_invalid_inputs(self):
        valid = ([5, 0, 0], [0, 2, 0], 125, 0.25, 100)
        for index, invalid in [(0, [0, 0, 0]), (0, [5, 0]), (0, [np.inf, 0, 0]),
                               (1, [0, np.nan, 0]), (1, [[0, 2, 0]]),
                               (2, 0), (2, np.inf), (3, -0.1), (3, np.nan),
                               (4, 0), (4, -100)]:
            with self.subTest(index=index, invalid=invalid):
                args = list(valid)
                args[index] = invalid
                with self.assertRaises(ValueError):
                    acceleration(*args)


if __name__ == '__main__':
    unittest.main()
