"""Offline Poynting–Robertson inspiral study. See METHOD.md."""
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from dust_bridge import acceleration


MU = 132712440000.0  # km^3/s^2, IAU nominal solar mass parameter
C = 299792.458  # km/s, exact
AU = 149597870.700  # km, exact
YEAR = 31557600.0  # s, Julian year
SUN = 695700.0  # km, nominal photospheric radius
NORMALIZATION = C * AU**2 / (4 * MU * YEAR)
DIRECT_RTOLS = (2e-9, 2e-11)
POLAR_RTOLS = (2e-9, 2e-11)
FULL_CASES = (
    (0.05, 0.60, 'circular'),
    (0.20, 0.60, 'circular'),
    (0.20, 0.85, 'circular'),
    (1.00, 0.85, 'circular'),
    (1.00, 0.9999, 'circular'),
    (1.00, 0.999999, 'circular'),
    (0.05, 0.35, 'gravity_only_release'),
)


def scales(a_au, beta):
    r0 = a_au * AU
    v0 = np.sqrt(MU * (1 - beta) / r0)
    return r0, v0, r0 / v0, beta * MU / (C * r0 * v0)


def secular_years(a_au, beta, q):
    return NORMALIZATION * a_au**2 * (1 - q**2) / beta


def bridge_quadrature(a_au, beta, order, phases):
    """Numerically integrate inverse drift reconstructed from bridge power."""
    r0, _, _, _ = scales(a_au, beta)
    m = MU * (1 - beta)
    k = beta * MU / C
    z, w = leggauss(order)
    radii = SUN + (r0 - SUN) * (z + 1) / 2
    weights = (r0 - SUN) * w / 2
    # A tilted plane; avoid depending on one Cartesian orientation.
    p = np.array([1., 0., 0.])
    b = np.array([0., 0.6, 0.8])
    rows = []
    for r, weight in zip(radii, weights):
        rates = []
        for phi in 2 * np.pi * np.arange(phases) / phases:
            er = np.cos(phi) * p + np.sin(phi) * b
            et = -np.sin(phi) * p + np.cos(phi) * b
            v = np.sqrt(m / r) * et
            force = acceleration(r * er, v, MU, beta, C)
            power = np.dot(v, force)
            rates.append(2 * r**2 * power / m)
        drift = float(np.mean(rates))
        if not drift < 0:
            raise RuntimeError('Expected inward circular energy drift')
        rows.append({'radius_km': r, 'weight_km': weight,
                     'phase_da_dt_km_s': rates, 'mean_da_dt_km_s': drift,
                     'contribution_s': -weight / drift,
                     'relative_error_vs_minus_2k_over_r': drift / (-2 * k / r) - 1})
    return {'order': order, 'phases': phases, 'nodes': rows,
            'time_years': sum(row['contribution_s'] for row in rows) / YEAR}


def make_events(cuts, radius_of):
    events = []
    for i, q in enumerate(cuts):
        def event(t, y, q=q):
            return radius_of(y) - q
        event.direction = -1
        event.terminal = i == len(cuts) - 1
        events.append(event)
    return events


def cartesian(a_au, beta, launch, rtol, cuts=(), stop=None):
    """Unaveraged Cartesian bridge; independent variable is physical time/t0."""
    r0, v0, t0, lam = scales(a_au, beta)
    s = 1.0 if launch == 'circular' else 1 / np.sqrt(1 - beta)
    y0 = [1., 0., 0., 0., s, 0., 0.]

    def rhs(t, y):
        r, v = y[:3], y[3:6]
        acc = acceleration(r0 * r, v0 * v, MU, beta, C)
        phi_dot = (r[0] * v[1] - r[1] * v[0]) / np.dot(r, r)
        return np.r_[v, (t0 / v0) * acc, phi_dot]

    atol = rtol / 100
    events = make_events(cuts, lambda y: np.linalg.norm(y[:3])) if cuts else None
    upper = 32 / (4 * lam) if stop is None else stop
    sol = solve_ivp(rhs, (0., upper), y0, method='DOP853', rtol=rtol,
                    atol=atol, events=events)
    if sol.status != (1 if cuts else 0):
        raise RuntimeError(f'Cartesian integration did not reach endpoint: {sol.message}')
    h_z = sol.y[0] * sol.y[4] - sol.y[1] * sol.y[3]
    torque_residual = abs(h_z - (s - lam * sol.y[6]))
    result = {'method': 'DOP853', 'rtol': rtol, 'atol': atol,
              'nfev': sol.nfev, 'accepted_steps': len(sol.t) - 1,
              'max_step_angle_rad': float(np.max(np.diff(sol.y[6]))),
              'max_sampled_speed_over_c': v0 / C * np.max(np.linalg.norm(sol.y[3:6], axis=0)),
              'max_scaled_angular_momentum_identity_residual': float(np.max(torque_residual))}
    if not cuts:
        r, v = sol.y[:3, -1], sol.y[3:6, -1]
        e = np.dot(v, v) / 2 - 1 / np.linalg.norm(r)
        a_ratio = -1 / (2 * e)
        result.update({'end_state_scaled': sol.y[:, -1].tolist(),
                       'duration_years': sol.t[-1] * t0 / YEAR,
                       'end_a_over_r0': a_ratio,
                       'negative_change_in_a_squared_over_r0_squared': 1 - a_ratio**2,
                       'expected_negative_change': 4 * lam * sol.t[-1],
                       'measured_to_secular_decay_ratio': (1 - a_ratio**2) / (4 * lam * sol.t[-1])
                       if lam else None})
        return result
    result['events'] = []
    for q, times, states in zip(cuts, sol.t_events, sol.y_events):
        if not len(times):
            raise RuntimeError(f'Missing Cartesian crossing at q={q}')
        y = states[0]
        result['events'].append({'q': q, 'radius_km': q * r0, 'crossing_count': len(times),
                                 'time_years': times[0] * t0 / YEAR,
                                 'state_scaled': y.tolist(),
                                 'radial_velocity_km_s': np.dot(y[:3], y[3:6]) * v0 / q})
    return result


def polar_reference(a_au, beta, launch, rtol, cuts):
    """Unaveraged Binet equations in angle; no bridge or Cartesian RHS calls."""
    r0, v0, t0, lam = scales(a_au, beta)
    s = 1.0 if launch == 'circular' else 1 / np.sqrt(1 - beta)

    def rhs(theta, y):
        u, du, _ = y
        j = s - lam * theta
        return du, 1 / j**2 - u - lam * du / j, 1 / (j * u**2)

    events = make_events(cuts, lambda y: 1 / y[0])
    atol = rtol / 100
    max_step = np.pi / 3
    sol = solve_ivp(rhs, (0., (s - 1e-7) / lam), (1., 0., 0.),
                    method='RK45', rtol=rtol, atol=atol,
                    max_step=max_step, events=events)
    if sol.status != 1:
        raise RuntimeError(f'Polar integration did not reach endpoint: {sol.message}')
    rows = []
    for q, angles, states in zip(cuts, sol.t_events, sol.y_events):
        if not len(angles):
            raise RuntimeError(f'Missing polar crossing at q={q}')
        theta = angles[0]
        j = s - lam * theta
        y = states[0]
        rows.append({'q': q, 'crossing_count': len(angles), 'angle_rad': theta,
                     'scaled_u_dU_dtheta_time': y.tolist(), 'scaled_h': j,
                     'radial_velocity_km_s': -v0 * j * y[1],
                     'time_years': y[2] * t0 / YEAR})
    return {'method': 'RK45', 'rtol': rtol, 'atol': atol, 'max_step_rad': max_step,
            'nfev': sol.nfev, 'accepted_steps': len(sol.t) - 1, 'events': rows}


def full_case(a_au, beta, launch):
    r0, v0, t0, lam = scales(a_au, beta)
    q_sun = SUN / r0
    cuts = tuple(q for q in (0.7, 0.3, 0.1, 0.03) if q > q_sun) + (q_sun,)
    direct = [cartesian(a_au, beta, launch, tol, cuts=cuts) for tol in DIRECT_RTOLS]
    polar = [polar_reference(a_au, beta, launch, tol, cuts) for tol in POLAR_RTOLS]
    summaries = []
    for q, d0, d1, p0, p1 in zip(cuts, direct[0]['events'], direct[1]['events'],
                                polar[0]['events'], polar[1]['events']):
        t, tc, ref, refc = (row['time_years'] for row in (d1, d0, p1, p0))
        comp = secular_years(a_au, beta, q) if launch == 'circular' else None
        summaries.append({'q': q, 'direct_refinement_relative': tc / t - 1,
                          'polar_refinement_relative': refc / ref - 1,
                          'direct_vs_polar_relative': t / ref - 1,
                          'circular_secular_time_years': comp,
                          'direct_vs_secular_relative': t / comp - 1 if comp is not None else None,
                          'direct_minus_secular_years': t - comp if comp is not None else None})
    return {'r0_au': a_au, 'beta': beta, 'launch': launch, 'v0_km_s': v0, 't0_s': t0,
            'lambda0': lam, 'lambda_circular_at_solar_surface': lam / np.sqrt(q_sun),
            'omitted_fraction_if_surface_time_is_approximated_by_zero_cutoff_law': q_sun**2,
            'direct': direct, 'polar': polar, 'comparisons': summaries}


def fit_law(rows):
    # Fit only the prefactor and exponents after removing the explicit cutoff term.
    design = np.array([[1., np.log(a), np.log(b)] for a, b, _ in rows])
    values = np.log([t for _, _, t in rows])
    coeff, _, rank, singular = np.linalg.lstsq(design, values, rcond=None)
    return {'n': len(rows), 'rank': int(rank), 'singular_values': singular.tolist(),
            'K_years': np.exp(coeff[0]), 'r0_au_exponent': coeff[1], 'beta_exponent': coeff[2],
            'max_abs_log_residual': np.max(abs(design @ coeff - values))}


def main():
    radii = (0.03, 0.1, 1.0, 5.0, 10.0)
    betas = (0.01, 0.1, 0.5, 0.9)
    grid, fit_rows = [], []
    for a in radii:
        q = SUN / (a * AU)
        for beta in betas:
            _, _, _, lam = scales(a, beta)
            quad = [bridge_quadrature(a, beta, n, p) for n, p in ((4, 4), (8, 8))]
            coarse, fine = (row['time_years'] for row in quad)
            predicted = secular_years(a, beta, q)
            stop = min(16 * 2 * np.pi, 0.08 / (4 * lam))
            windows = [cartesian(a, beta, 'circular', tol, stop=stop) for tol in DIRECT_RTOLS]
            grid.append({'r0_au': a, 'beta': beta, 'q_sun': q,
                         'secular_years': predicted, 'quadrature': quad,
                         'quad_coarse_vs_fine_relative': coarse / fine - 1,
                         'quad_fine_vs_secular_relative': fine / predicted - 1,
                         'cartesian_windows': windows})
            fit_rows.append((a, beta, fine / (1 - q**2)))

    fit = fit_law(fit_rows)
    control = [cartesian(1.0, 0.0, 'circular', tol, stop=16 * 2 * np.pi)
               for tol in DIRECT_RTOLS]
    full = []
    for a, beta, launch in FULL_CASES:
        print(f'Full trajectory: {a:g} AU, beta={beta:g}, {launch}', flush=True)
        full.append(full_case(a, beta, launch))
    result = {'environment': {'python': platform.python_version(), 'numpy': np.__version__,
                              'scipy': scipy.__version__},
              'analytic_K_years': NORMALIZATION, 'quadrature_fit': fit,
              'grid': grid, 'zero_drag_control': control, 'full_cases': full}
    output = Path('results.json')
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(f'Wrote {output}; quadrature fit: K={fit["K_years"]:.6f} yr, '
          f'radius exponent={fit["r0_au_exponent"]:.6f}, '
          f'beta exponent={fit["beta_exponent"]:.6f}')


if __name__ == '__main__':
    main()
