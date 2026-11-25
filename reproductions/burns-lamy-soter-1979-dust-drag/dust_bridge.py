'''Heliocentric gravity and first-order Poynting–Robertson acceleration.'''

import numpy as np
from sim.dynamics.orbit.accelerations import accel_two_body


def _vector3(value, name):
    raw = np.asarray(value)
    if raw.shape != (3,) or np.iscomplexobj(raw):
        raise ValueError(f'{name} must be a finite real 3-vector')
    vector = raw.astype(float, copy=False)
    if not np.all(np.isfinite(vector)):
        raise ValueError(f'{name} must be a finite real 3-vector')
    return vector


def _scalar(value, name, *, allow_zero=False):
    raw = np.asarray(value)
    if raw.ndim != 0 or np.iscomplexobj(raw):
        raise ValueError(f'{name} must be a finite real scalar')
    number = float(raw)
    if not np.isfinite(number) or number < 0 or (number == 0 and not allow_zero):
        domain = 'nonnegative' if allow_zero else 'positive'
        raise ValueError(f'{name} must be finite and {domain}')
    return number


def acceleration(r_km, v_km_s, mu_km3_s2, beta, c_km_s):
    '''Return total heliocentric inertial acceleration as a 3-vector in km/s^2.

    Position is in km, velocity and light speed are in km/s, mu is in
    km^3/s^2, and beta is dimensionless. The origin is excluded. The
    radiation force retains terms through first order in v/c and is
    physically applicable when the speed is small compared to c.
    '''
    r = _vector3(r_km, 'r_km')
    v = _vector3(v_km_s, 'v_km_s')
    mu = _scalar(mu_km3_s2, 'mu_km3_s2')
    beta = _scalar(beta, 'beta', allow_zero=True)
    c = _scalar(c_km_s, 'c_km_s')

    radius_km = np.hypot.reduce(r)
    if radius_km == 0:
        raise ValueError('r_km must be nonzero')

    gravity_km_s2 = np.asarray(accel_two_body(r, mu), dtype=float)
    if beta == 0:
        return gravity_km_s2

    r_hat = r / radius_km
    v_over_c = v / c
    radial_v_over_c = np.dot(v_over_c, r_hat)
    radiation_scale_km_s2 = beta * ((mu / radius_km) / radius_km)
    # Direct pressure is -beta times central gravity. Both velocity terms
    # enter with a minus sign; the radial projection is signed.
    return ((1.0 - beta) * gravity_km_s2
            - radiation_scale_km_s2 * (v_over_c + radial_v_over_c * r_hat))
