"""Fixed-orbit planar spin bridge for OEL revision 76ad5422bb476dd6cd3e162fc6346b3f7af6c069.

Time is in seconds, angles in radians, and inertial spin rates in rad/s.
The orbit is prograde in the inertial xy plane, with periapsis on +x at
t=0. The orbit is prescribed: neither spin torque changes it.

Let psi = theta - f, with primes denoting d/df and dots denoting d/dt.
Then ddot(theta) = f_dot**2 * psi'' + f_ddot * (1 + psi'). At a fixed
state the incremental tide enters only through the first term. Keplerian
kinematics gives f_dot = n * sqrt(1-e**2) * (a/r)**2 > 0, hence
sign(psi') = sign(theta_dot - f_dot). The supplied anomaly-domain tide
therefore gives the inertial angular acceleration

    alpha_tide = -D * n**2 * (a/r)**6 * sign(theta_dot - f_dot).

D and a/r are dimensionless, and n = sqrt(mu/a**3) is in s^-1.
The result is in rad/s**2; multiplying by C in kg m**2 gives torque
in N m. There is no smoothing or sticking rule: the tide is zero at
exact equality and changes sign there, not at the 3:2 spin rate.
"""

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.optimize import brentq

from sim.dynamics.attitude.disturbances import gravity_gradient_torque_body_nm
from sim.dynamics.orbit.elements import coe_to_rv_eci


_INV_SQRT_2 = math.sqrt(0.5)
_C_BN_ZERO = np.eye(3)
_C_BN_QUARTER = np.array([
    [_INV_SQRT_2, _INV_SQRT_2, 0.0],
    [-_INV_SQRT_2, _INV_SQRT_2, 0.0],
    [0.0, 0.0, 1.0],
])


@dataclass(frozen=True)
class OrbitSample:
    """One periodic orbital state and its OEL axial torque basis."""

    true_anomaly_rad: float  # wrapped to [0, 2*pi)
    f_dot_rad_s: float
    a_over_r: float
    torque_at_zero_nm: float
    torque_at_pi_over_4_nm: float


@dataclass(frozen=True)
class SpinBridge:
    """Study parameters and scalar-time, array-capable spin kernels.

    The default Sun/1-au orbit sets the conversion to seconds; it is not
    an extra dimensionless study parameter. The inertia normalization
    uses A/C = 0.5 and B/C = 0.5 + triaxiality. Only B-A enters the
    planar axial gravity-gradient torque.
    """

    eccentricity: float = 0.2
    triaxiality: float = 1e-4  # (B-A)/C
    D: float = 1e-5
    a_km: float = 149_597_870.7
    mu_km3_s2: float = 132_712_440_041.27942
    C_kg_m2: float = 1.0
    tide_enabled: bool = True

    def __post_init__(self):
        if not 0.0 <= self.eccentricity < 1.0:
            raise ValueError("eccentricity must lie in [0, 1)")
        if not 0.0 <= self.triaxiality <= 0.5:
            raise ValueError("triaxiality must lie in [0, 0.5] for this inertia normalization")
        if not math.isfinite(self.D) or self.D < 0.0:
            raise ValueError("D must be finite and nonnegative")
        for name in ("a_km", "mu_km3_s2", "C_kg_m2"):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")

    @property
    def mean_motion_rad_s(self):
        return math.sqrt(self.mu_km3_s2 / self.a_km**3)

    @property
    def period_s(self):
        return math.tau / self.mean_motion_rad_s

    def orbit_sample(self, t_s):
        """Sample the orbit and torque basis at a scalar inertial time."""
        t_s = float(t_s)
        if not math.isfinite(t_s):
            raise ValueError("t_s must be finite")
        return self._sample_at_phase(t_s % self.period_s)

    @lru_cache(maxsize=512)
    def _sample_at_phase(self, phase_s):
        # Cache exact requests, without rounding or interpolating phase.
        e = self.eccentricity
        mean_anomaly = self.mean_motion_rad_s * phase_s
        if e == 0.0 or mean_anomaly == 0.0:
            E = mean_anomaly
        else:
            E = brentq(
                lambda x: x - e * math.sin(x) - mean_anomaly,
                mean_anomaly - e,
                mean_anomaly + e,
                xtol=1e-15,
            )
        f = math.atan2(
            math.sqrt((1.0 - e) * (1.0 + e)) * math.sin(E),
            math.cos(E) - e,
        ) % math.tau
        r, v = coe_to_rv_eci(
            a_km=self.a_km, ecc=e, inc_deg=0.0, raan_deg=0.0,
            argp_deg=0.0, true_anomaly_deg=math.degrees(f) % 360.0,
            mu_km3_s2=self.mu_km3_s2,
        )
        radius_km = float(np.linalg.norm(r))
        h_z = float(r[0] * v[1] - r[1] * v[0])
        C = self.C_kg_m2
        inertia = np.diag([0.5 * C, (0.5 + self.triaxiality) * C, C])

        def axial_torque(c_bn):
            torque = gravity_gradient_torque_body_nm(
                mu_km3_s2=self.mu_km3_s2, inertia_kg_m2=inertia,
                position_km=r, c_bn=c_bn,
            )
            return float(torque[2])

        return OrbitSample(
            true_anomaly_rad=f,
            f_dot_rad_s=h_z / radius_km**2,
            a_over_r=self.a_km / radius_km,
            torque_at_zero_nm=axial_torque(_C_BN_ZERO),
            torque_at_pi_over_4_nm=axial_torque(_C_BN_QUARTER),
        )

    def conservative_accel_rad_s2(self, sample, theta_rad):
        """OEL gravity-gradient acceleration at inertial body angle theta."""
        # With diagonal inertia, the axial torque is proportional to the
        # product of the two in-plane body-frame radial components. A
        # rotation by theta makes that product a pure second harmonic.
        # OEL evaluations at theta=0 and pi/4 supply its two coefficients.
        twice_theta = 2.0 * np.asarray(theta_rad, dtype=float)
        torque_nm = (
            sample.torque_at_zero_nm * np.cos(twice_theta)
            + sample.torque_at_pi_over_4_nm * np.sin(twice_theta)
        )
        return torque_nm / self.C_kg_m2

    def tide_accel_rad_s2(self, sample, theta_dot_rad_s):
        """Constant angular-lag acceleration given the inertial spin rate."""
        spin = np.asarray(theta_dot_rad_s, dtype=float)
        if not self.tide_enabled or self.D == 0.0:
            return np.zeros_like(spin)
        scale = self.D * self.mean_motion_rad_s**2 * sample.a_over_r**6
        return -scale * np.sign(spin - sample.f_dot_rad_s)

    def angular_accel_rad_s2(self, t_s, theta_rad, theta_dot_rad_s):
        """Sum both accelerations, broadcasting spins at one orbital time."""
        sample = self.orbit_sample(t_s)
        return (
            self.conservative_accel_rad_s2(sample, theta_rad)
            + self.tide_accel_rad_s2(sample, theta_dot_rad_s)
        )
