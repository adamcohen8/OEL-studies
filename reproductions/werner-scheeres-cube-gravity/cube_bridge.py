"""Newtonian acceleration for a homogeneous, axis-aligned cube.

The cube is centered at the origin. Inputs are in km and km**3/s**2;
outputs are in km/s**2. The numerical domain is the closed box
abs(position[i] / side) <= 4 for each coordinate, which includes the full
study region within four side lengths of the center.

We normalize to side 1. For field point q and source point s,
(s_i - q_i)/|s-q|**3 = -d(1/|s-q|)/ds_i. The divergence theorem thus gives
component i as the lower-face integral of 1/|s-q| minus the upper-face
integral, multiplied by mu/side**2. If q is on or in the cube, excising a
ball of radius eps adds a boundary term bounded by 4*pi*eps, which tends
to zero. The remaining face integrals are finite even when q lies on a
face, edge, or corner. The primitive below uses its continuous limits
there, so no displacement of the query point is needed. In particular,
the exact face-center result is also the limit from either side.

The face differences are computed in double precision. Near a symmetry
plane, absolute accuracy in units of mu/side**2 is the relevant measure;
relative accuracy is not guaranteed for arbitrarily small components.
"""

import math

import numpy as np


def _primitive(u, v, d):
    """Corner primitive for 1/sqrt(u*u + v*v + d*d), with d >= 0."""
    # u*asinh(v/hypot(u,d)) has continuous limit zero at u=0,
    # including when d=0. The analogous rule holds for the v term.
    a = u * math.asinh(v / math.hypot(u, d)) if u != 0.0 else 0.0
    b = v * math.asinh(u / math.hypot(v, d)) if v != 0.0 else 0.0
    # Using a nonnegative distance keeps the angle on the correct branch
    # regardless of which side of the face contains the query point.
    c = d * math.atan2(u * v, d * math.hypot(u, v, d)) if d != 0.0 else 0.0
    return math.fsum((a, b, -c))


def _face(d, u0, u1, v0, v1):
    return math.fsum((
        _primitive(u1, v1, d), -_primitive(u0, v1, d),
        -_primitive(u1, v0, d), _primitive(u0, v0, d),
    ))


def cube_accel(position, side, mu):
    """Return the acceleration at one point as a float64 array of shape (3,).

    position is a finite Cartesian triple in km, side is the finite positive
    side length in km, and mu is the finite positive total GM in km**3/s**2.
    Interior, exterior and exact boundary points are accepted within the
    normalized domain above. Invalid shapes, values or positions outside
    that domain raise ValueError. A component too large for float64 raises
    OverflowError; sufficiently small components may underflow to zero.
    """
    p = np.asarray(position, dtype=np.float64)
    side, mu = float(side), float(mu)
    if p.shape != (3,) or not np.all(np.isfinite(p)):
        raise ValueError('position must be a finite Cartesian triple')
    if not (math.isfinite(side) and side > 0.0 and math.isfinite(mu) and mu > 0.0):
        raise ValueError('side and mu must be finite and positive')
    q = tuple(float(x) / side for x in p)
    if any(abs(x) > 4.0 for x in q):
        raise ValueError('each normalized position coordinate must be in [-4, 4]')

    # Separate mantissas from exponents so mu/side**2 need not itself
    # be representable, even when the requested components are.
    sm, se = math.frexp(side)
    mm, me = math.frexp(mu)
    scale_m = mm / (sm * sm)
    scale_e = me - 2 * se
    out = np.zeros(3, dtype=np.float64)
    for i, x in enumerate(q):
        if x == 0.0:
            continue  # Exact reflection symmetry for this component.
        y, z = q[(i + 1) % 3], q[(i + 2) % 3]
        bounds = (-0.5 - y, 0.5 - y, -0.5 - z, 0.5 - z)
        lower = _face(abs(-0.5 - x), *bounds)
        upper = _face(abs(0.5 - x), *bounds)
        out[i] = math.ldexp((lower - upper) * scale_m, scale_e)
    return out
