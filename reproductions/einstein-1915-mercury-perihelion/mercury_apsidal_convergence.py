#!/usr/bin/env python3
"""Native OEL Sun–Mercury apsidal-rate study; Python does not integrate the orbit.

Usage: python mercury_apsidal_convergence.py --output result.json

Both branches use point-mass solar gravity and explicitly set the documented
Schwarzschild switch. The script checks the normalized configuration and the
effective solar GM. The replay launcher verifies checkout and package identity.
"""

import argparse
from collections.abc import Mapping
from copy import deepcopy
from dataclasses import asdict, is_dataclass
import json
import math
from pathlib import Path
import sys


AU_M = 149597870700
C_M_S = 299792458
DAY_S = 86400
JULIAN_YEAR_S = 36525 * DAY_S // 100
CENTURY_S = 100 * JULIAN_YEAR_S
INITIAL_JD_UTC = 2451545.0
BASELINE_YEARS = 20
BASELINE_S = BASELINE_YEARS * JULIAN_YEAR_S
END_JD_UTC = INITIAL_JD_UTC + BASELINE_S / DAY_S
A_AU, E = 0.38709927, 0.20563593
A_M = A_AU * AU_M
ARCSEC_PER_RAD = 180 * 3600 / math.pi
RATE_SCALE = ARCSEC_PER_RAD * CENTURY_S
NATIVE_ROOT_XTOL_S = 1e-6
HERMITE_ROOT_XTOL_U = 1e-13
ROOT_RTOL = 4 * sys.float_info.epsilon
ROOT_MAXITER = 100
LIMIT = 0.001
VERSION = '1.2.0+preview.1'
DEVELOPMENT_COMMIT = 'b619c3b790206ef629b3e6121e45ee5a20186857'
PUBLIC_PREVIEW_COMMIT = 'e3b0e7620d4ae616d01da10b7ee23ff282d17800'
DOC_URL = ('https://github.com/adamcohen8/orbital-engineering-lab/blob/'
           + PUBLIC_PREVIEW_COMMIT + '/docs/planetary-onp.md')
SWITCH_PATH = ('simulator', 'dynamics', 'orbit', 'schwarzschild', 'enabled')


class Blocker(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise Blocker(message)


def at(value, path):
    for key in path:
        if isinstance(value, (Mapping, list, tuple)):
            value = value[key]
        else:
            value = getattr(value, key) if isinstance(key, str) else value[key]
    return value


def resolved_config(cfg):
    if hasattr(cfg, 'model_dump'):
        raw = cfg.model_dump(mode='json')
    elif is_dataclass(cfg):
        raw = asdict(cfg)
    elif hasattr(cfg, 'json'):
        raw = json.loads(cfg.json())
    elif isinstance(cfg, Mapping):
        raw = dict(cfg)
    else:
        raise Blocker('No supported public serialization for the resolved configuration.')
    # Round-trip so comparison is on the same values that are archived.
    try:
        return json.loads(json.dumps(raw, allow_nan=False))
    except (TypeError, ValueError) as exc:
        raise Blocker('Resolved configuration is not JSON serializable.') from exc


def write_report(path, report):
    staging = path.with_name(path.name + '.tmp')
    with staging.open('w', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    staging.replace(path)


def config_input(x0, n, tol, enabled):
    return {
        'scenario_name': 'mercury_apsidal_convergence',
        'scenario_description': 'Isolated heliocentric point-Sun test particle',
        'objects': {'probe': {'kind': 'satellite', 'runtime_profile': 'trajectory_only',
                              'specs': {'mass_kg': 1.0}, 'initial_state': {
                                  'position_inertial_km': x0[:3].tolist(),
                                  'velocity_inertial_km_s': x0[3:].tolist()}}},
        'simulator': {'central_body': 'sun', 'initial_jd_utc': INITIAL_JD_UTC,
                      'duration_s': float(BASELINE_S), 'dt_s': BASELINE_S / n,
                      'planetary': {'gravity': {'model': 'point_mass'}},
                      'dynamics': {'orbit': {
                          'integrator': 'rkf78', 'numeric_backend': 'rust',
                          'adaptive_atol': tol, 'adaptive_rtol': tol,
                          'schwarzschild': {'enabled': enabled}}}},
        'outputs': {'stats': {'save_full_log': False, 'save_csv': False},
                    'output_dir': 'outputs/mercury_apsidal_convergence'},
    }


def radial(x):
    return float(x[0]*x[3] + x[1]*x[4] + x[2]*x[5])


def angle(x):
    return math.atan2(float(x[1]), float(x[0]))


def event_fit(events, np):
    require(len(events) >= 4, 'Insufficient perihelia for the full-baseline and split-half fits.')
    times, angles = np.asarray(events, dtype=float).T
    require(bool(np.all((0 < times) & (times < BASELINE_S))),
            'Perihelion times must lie strictly inside the integration baseline.')
    require(bool(np.all(np.diff(times) > 0)), 'Perihelion times are not strictly increasing.')
    unwrapped = np.unwrap(angles)

    def fit(t, y):
        require(len(t) >= 2, 'Too few events for a split-half fit.')
        centered = t - t.mean()
        dy = y - y.mean()
        slope = float(np.dot(centered, dy) / np.dot(centered, centered))
        return slope * RATE_SCALE, dy - slope * centered

    slope, residual = fit(times, unwrapped)
    first = times < BASELINE_S / 2
    s1, _ = fit(times[first], unwrapped[first])
    s2, _ = fit(times[~first], unwrapped[~first])
    return {'rate_arcsec_per_julian_century': slope, 'event_count': len(events),
            'first_event_s': float(times[0]), 'last_event_s': float(times[-1]),
            'event_time_s': times.tolist(), 'event_angle_unwrapped_rad': unwrapped.tolist(),
            'residual_rms_arcsec': float(np.sqrt(np.mean(residual**2))) * ARCSEC_PER_RAD,
            'residual_max_abs_arcsec': float(np.max(np.abs(residual))) * ARCSEC_PER_RAD,
            'first_half_event_count': int(np.count_nonzero(first)),
            'second_half_event_count': int(np.count_nonzero(~first)),
            'first_half_rate_arcsec_per_julian_century': s1,
            'second_half_rate_arcsec_per_julian_century': s2,
            'second_minus_first_arcsec_per_julian_century': s2 - s1}


def bracketed_root(brentq, g, upper, left_value, right_value, xtol):
    require(math.isfinite(upper) and math.isfinite(xtol) and 0 < xtol < upper,
            'Root tolerance must be positive and smaller than the local bracket span.')
    require(math.isfinite(left_value) and math.isfinite(right_value)
            and left_value < 0 <= right_value, 'Invalid perihelion root bracket.')
    root, info = brentq(g, 0.0, upper, xtol=xtol, rtol=ROOT_RTOL,
                       maxiter=ROOT_MAXITER, full_output=True, disp=False)
    require(info.converged and 0 <= root <= upper,
            'Perihelion root did not converge within its bracket and iteration limit.')
    return root


def sampled_events(times, states, stride, brentq):
    events = []
    for j in range(stride, len(times), stride):
        i = j - stride
        x, y = states[i], states[j]
        g0, g1 = radial(x), radial(y)
        if not g0 < 0 <= g1:
            continue
        h = times[j] - times[i]
        require(0 <= times[i] < times[j] <= BASELINE_S,
                'Hermite root bracket lies outside the integration baseline.')
        # Cubic in u=(t-t_i)/h. Velocity is the derivative of this same cubic.
        p, q = x[:3], y[:3]
        b, v = h*x[3:], h*y[3:]
        d = 3*(q-p) - 2*b - v
        f = 2*(p-q) + b + v

        def position(u):
            return p + u*(b + u*(d + u*f))

        def g(u):
            if u == 0:
                return g0
            if u == 1:
                return g1
            return float(position(u) @ (b + u*(2*d + 3*u*f))) / h

        u = bracketed_root(brentq, g, 1.0, g0, g1, HERMITE_ROOT_XTOL_U)
        t = float(times[j] if u == 1 else times[i] + h*u)
        if 0 < t < BASELINE_S:
            events.append((t, angle(position(u))))
    return events


def propagate_branch(prop, x0, n, save_samples, np, brentq):
    def step(x, h, t):
        result = np.asarray(prop.propagate(x.copy(), h, t), dtype=float)
        require(result.shape == (6,) and bool(np.isfinite(result).all()),
                'propagate did not return a finite primary-relative inertial six-vector.')
        return result.copy()

    times = np.empty(n + 1) if save_samples else None
    states = np.empty((n + 1, 6)) if save_samples else None
    x, t, g0 = x0.copy(), 0.0, radial(x0)
    events = []
    if save_samples:
        times[0], states[0] = t, x
    for j in range(1, n + 1):
        end = (j * BASELINE_S) / n
        h = end - t
        require(0 <= t < end <= BASELINE_S, 'Native step lies outside the integration baseline.')
        y = step(x, h, t)
        g1 = radial(y)
        if g0 < 0 <= g1:
            # Use native partial steps without replacing the accepted endpoint.
            def g(s):
                if s == 0:
                    return g0
                if s == h:
                    return g1
                return radial(step(x, s, t))

            s = bracketed_root(brentq, g, h, g0, g1, NATIVE_ROOT_XTOL_S)
            te = end if s == h else t + s
            if 0 < te < BASELINE_S:
                events.append((te, angle(y if s == h else step(x, s, t))))
        x, t, g0 = y, end, g1
        if save_samples:
            times[j], states[j] = t, x
    result = event_fit(events, np)
    outputs = {}
    if save_samples:
        for stride in (4, 2, 1):
            outputs[stride] = event_fit(sampled_events(times, states, stride, brentq), np)
    return result, outputs


def rate(fit):
    return fit['rate_arcsec_per_julian_century']


def pair(fits):
    return {'branches': fits, 'difference_arcsec_per_julian_century':
            rate(fits['schwarzschild']) - rate(fits['newtonian'])}


def delta(pair_record):
    return pair_record['difference_arcsec_per_julian_century']


def run(report, checkpoint):
    import numpy as np
    from scipy.optimize import brentq
    from sim.config.scenario_yaml import scenario_config_from_dict
    from sim.dynamics.planetary.definitions import BODIES
    from sim.dynamics.planetary.preparation import prepare_environment
    from sim.dynamics.planetary.propagation import PlanetaryForcePlan, PlanetaryPropagator

    require(JULIAN_YEAR_S == 31557600 and CENTURY_S == 3155760000
            and BASELINE_S == 631152000 and END_JD_UTC - INITIAL_JD_UTC == 7305,
            'Inconsistent conversion between elapsed seconds, Julian years, and Julian days.')
    require(math.isclose(math.tau * ARCSEC_PER_RAD, 1296000, rel_tol=1e-14)
            and math.isclose(RATE_SCALE / 100, ARCSEC_PER_RAD * JULIAN_YEAR_S,
                             rel_tol=1e-14),
            'Inconsistent angular conversion or annualization factor.')
    require(NATIVE_ROOT_XTOL_S >= math.ulp(float(BASELINE_S)),
            'Native root absolute tolerance is finer than elapsed-time resolution at the baseline end.')
    mu = float(BODIES['sun'].mu_km3_s2)
    require(math.isfinite(mu) and mu > 0, 'Solar catalog GM must be finite and positive.')
    mu_m3_s2 = mu * 1e9
    report['body_catalog_solar_mu_km3_s2'] = mu
    report['constants']['solar_mu_m3_s2'] = mu_m3_s2
    period_s = 2 * math.pi * math.sqrt(A_M**3 / mu_m3_s2)
    rp = A_M * (1 - E)
    x0 = np.array([rp / 1000, 0, 0, 0,
                   math.sqrt(mu_m3_s2 * (1 + E) / rp) / 1000, 0])
    report['initial_state_km_km_s'] = x0.tolist()
    report['invocation'] = {
        'setup': ['cfg = scenario_config_from_dict(input_config)',
                  'env = prepare_environment(cfg)', 'plan = PlanetaryForcePlan(env, mass_kg=1.0)',
                  "prop = PlanetaryPropagator(plan, cfg.simulator.dynamics['orbit'])"],
        'advance': 'prop.propagate(x.copy(), t_next - t, t)',
        'root_trial': 'prop.propagate(x.copy(), s, t), 0 < s < t_next - t',
        'state_frame': 'Sun-relative inertial',
        'state_units': ['km', 'km', 'km', 'km/s', 'km/s', 'km/s'],
        'time_origin': 'scenario start', 'time_unit': 's',
        'public_step_limit': 'Maximum outer native call span; not a claim about an internal max_step knob.',
        'integration_events': 'Native partial-step roots; no independent output interpolation on this sweep.',
        'output_events': 'Cubic Hermite roots on decimations of the same finest native trajectory.',
    }
    expected_halves = (math.floor((BASELINE_S / 2) / period_s),
                       math.floor(BASELINE_S / period_s))
    expected_halves = (expected_halves[0], expected_halves[1] - expected_halves[0])
    report['baseline']['kepler_expected_perihelia_by_half'] = list(expected_halves)
    n0 = math.ceil(64 * BASELINE_S / period_s)
    divisors = (64, 128, 256, 512, 1024, 2048, 4096)
    levels = [(k, 1e-11) for k in divisors]
    levels += [(divisors[-1], 1e-12), (divisors[-1], 1e-13)]
    sweep = report['integration_sweep'] = []
    sampled = {4: {}, 2: {}, 1: {}}
    for index, (divisor, tol) in enumerate(levels):
        n = n0 * (divisor // 64)
        fine = index == len(levels)-1
        entry = {'kepler_period_divisor': divisor, 'outer_step_count': n,
                 'outer_step_s': BASELINE_S/n, 'adaptive_atol': tol, 'adaptive_rtol': tol,
                 'inputs': {}, 'resolved': {}, 'effective_mu_km3_s2': {}}
        sweep.append(entry)
        setups = {}
        for branch, enabled in (('newtonian', False), ('schwarzschild', True)):
            raw = config_input(x0, n, tol, enabled)
            entry['inputs'][branch] = deepcopy(raw)
            cfg = scenario_config_from_dict(raw)
            resolved = resolved_config(cfg)
            entry['resolved'][branch] = resolved
            require(at(resolved, SWITCH_PATH) is enabled,
                    'Resolved Schwarzschild switch does not match the branch.')
            sim = resolved['simulator']
            require(sim['planetary']['gravity']['model'] == 'point_mass',
                    'Resolved gravity model is not point_mass.')
            require(sim['central_body'] == 'sun' and sim['initial_jd_utc'] == INITIAL_JD_UTC
                    and sim['duration_s'] == BASELINE_S
                    and math.isclose(sim['dt_s'], BASELINE_S/n, rel_tol=1e-14),
                    'Resolved central body, epoch, duration, or step differs from the request.')
            orbit = sim['dynamics']['orbit']
            require(orbit['integrator'] == 'rkf78' and orbit['numeric_backend'] == 'rust'
                    and orbit['adaptive_atol'] == tol and orbit['adaptive_rtol'] == tol,
                    'Resolved integration controls differ from the request.')
            env = prepare_environment(cfg)
            effective_mu = float(env.gravity.mu_km3_s2)
            entry['effective_mu_km3_s2'][branch] = effective_mu
            require(math.isclose(effective_mu, mu, rel_tol=1e-14),
                    'Resolved environment solar GM differs from the catalog GM used for initialization.')
            setups[branch] = PlanetaryPropagator(PlanetaryForcePlan(env, mass_kg=1.0),
                                                 cfg.simulator.dynamics['orbit'])
        check = deepcopy(entry['resolved']['newtonian'])
        at(check, SWITCH_PATH[:-1])[SWITCH_PATH[-1]] = True
        require(check == entry['resolved']['schwarzschild'],
                'Paired resolved configurations differ beyond the Schwarzschild Boolean switch.')
        checkpoint()
        fits = {}
        for branch, prop in setups.items():
            fits[branch], diagnostics = propagate_branch(prop, x0, n, fine, np, brentq)
            for stride, fit in diagnostics.items():
                sampled[stride][branch] = fit
        entry.update(pair(fits))
        if index:
            entry['difference_shift_from_previous'] = delta(entry) - delta(sweep[index-1])
        checkpoint()

    finest = sweep[-1]
    output = report['output_sweep'] = []
    for stride in (4, 2, 1):
        entry = {'stride': stride, 'sample_interval_s': stride * finest['outer_step_s'],
                 **pair(sampled[stride])}
        if output:
            entry['difference_shift_from_previous'] = delta(entry) - delta(output[-1])
        output.append(entry)

    limit = LIMIT
    step_shifts = [delta(sweep[j]) - delta(sweep[j-1])
                   for j in range(len(divisors) - 2, len(divisors))]
    tolerance_shifts = [delta(sweep[j]) - delta(sweep[j-1])
                        for j in range(len(divisors), len(sweep))]
    output_shift = delta(output[-1]) - delta(output[-2])
    native_rate = delta(finest)
    reconstruction_shift = delta(output[-1]) - native_rate
    counts_match = all(
        fit['event_count'] == finest['branches'][branch]['event_count']
        for record in sweep + output for branch, fit in record['branches'].items()
    )
    checks = {
        'last_two_step_refinements': all(abs(s) < limit for s in step_shifts),
        'last_two_tolerance_refinements': all(abs(s) < limit for s in tolerance_shifts),
        'finest_output_refinement': abs(output_shift) < limit,
        'finest_sampled_vs_native_root_agreement': abs(reconstruction_shift) < limit,
        'finest_newtonian_absolute_rate': abs(rate(finest['branches']['newtonian'])) < limit,
        'event_counts_match': counts_match,
        'split_half_counts_match_kepler_expectation': all(
            (fit['first_half_event_count'], fit['second_half_event_count']) == expected_halves
            for record in sweep + output for fit in record['branches'].values()
        ),
    }
    converged = all(checks.values())
    observed_scale = max([abs(s) for s in step_shifts + tolerance_shifts] +
                         [abs(output_shift), abs(reconstruction_shift),
                          abs(rate(finest['branches']['newtonian']))])
    report['convergence'] = {
        'passed': converged, 'checks': checks,
        'rate_unit': 'arcsec per Julian century', 'threshold': limit,
        'last_two_step_shifts': step_shifts, 'last_two_tolerance_shifts': tolerance_shifts,
        'finest_output_shift': output_shift,
        'finest_sampled_minus_native': reconstruction_shift,
        'observed_numerical_scale_not_an_error_bound': observed_scale,
    }
    places = (max(0, min(3, math.floor(-math.log10(observed_scale))))
              if observed_scale > 0 else 3)
    report['native_result'] = {
        **pair(finest['branches']),
        'status': 'converged' if converged else 'nonconverged',
        'reported_difference_arcsec_per_julian_century': (
            round(native_rate, places) if converged else None),
        'reported_decimal_places': places if converged else None,
        'difference_arcsec_per_julian_year': native_rate / 100,
        'reported_difference_arcsec_per_julian_year': (
            round(native_rate, places) / 100 if converged else None),
        'sign': 'Schwarzschild minus Newtonian; positive is prograde',
    }
    report['status'] = report['native_result']['status']
    checkpoint()

    advance = 6 * math.pi * mu_m3_s2 / (A_M * C_M_S**2 * (1 - E**2))
    prediction = advance * RATE_SCALE / period_s
    report['independent_check'] = {
        'c_m_s': C_M_S, 'advance_rad_per_orbit': advance, 'kepler_period_s': period_s,
        'weak_field_rate_arcsec_per_julian_century': prediction,
        'native_minus_weak_field_arcsec_per_julian_century': native_rate - prediction,
        'interpretation': ('Leading-order comparison using Newtonian elements and period; '
                           'not an exact reference for evolution under a truncated force law.'),
    }
    checkpoint()
    return converged


def main():
    parser = argparse.ArgumentParser(description='Paired native Mercury perihelion study')
    parser.add_argument('--output', '-o', type=Path, default=Path('mercury_apsidal_convergence.json'))
    args = parser.parse_args()

    report = {
        'schema': 'mercury_apsidal_convergence.v1', 'status': 'running',
        'requested_build': {'version': VERSION, 'revision': DEVELOPMENT_COMMIT,
                            'public_preview_revision': PUBLIC_PREVIEW_COMMIT,
                            'identity_verification': 'Replay launcher, not this script'},
        'force_model_provenance': {'configuration_documentation': DOC_URL,
                                   'schwarzschild_switch': '.'.join(SWITCH_PATH),
                                   'native_source_path': 'rust/oel-orbit/src/precision_forces.rs',
                                   'native_c_km_s': C_M_S / 1000},
        'constants': {'au_m': AU_M, 'c_m_s': C_M_S, 'day_s': DAY_S,
                      'julian_year_s': JULIAN_YEAR_S, 'julian_century_s': CENTURY_S,
                      'arcsec_per_rad': ARCSEC_PER_RAD},
        'baseline': {
            'initial_jd_utc': INITIAL_JD_UTC, 'end_jd_utc': END_JD_UTC,
            'duration_s': BASELINE_S, 'duration_julian_days': BASELINE_S / DAY_S,
            'duration_julian_years': BASELINE_YEARS,
            'first_second_half_boundary_s': BASELINE_S / 2,
            'split_definition': ('Two 10-Julian-year halves of the 20-Julian-year integration, '
                                 'not half centuries. An event at the midpoint belongs to the second half.'),
            'coverage': ('Both endpoints are interior to the published 1972–2030 Sun rotation '
                         'implementation envelope. Normal OEL environment preparation and validation apply.'),
            'rate_scale_from_rad_per_s': RATE_SCALE,
            'rate_scaling': ('Fit unwrapped perihelion angle in radians against elapsed seconds. '
                             'Multiply each slope by arcseconds per radian and seconds per Julian century. '
                             'Residuals are angles, so convert them only by arcseconds per radian.'),
            'annualization': ('Divide the per-Julian-century rate by 100 for arcseconds per Julian year. '
                              'The reported annualized value uses the rounded century-unit result; '
                              'it is not a separate annual fit.'),
            'limitations': ('The century unit expresses a slope measured over 20 Julian years; '
                            'no century-long trajectory or long-term constancy is tested. '
                            'The shorter baseline is more sensitive to endpoint selection and numerical drift '
                            'in the fitted slope. Split halves and refinements are diagnostics, not error bounds.'),
        },
        'root_finding': {
            'native_bracket': 'Elapsed seconds within the current native step: [0, h].',
            'native_xtol_s': NATIVE_ROOT_XTOL_S,
            'hermite_bracket': 'Dimensionless sample-interval coordinate: [0, 1], with t = t_i + h*u.',
            'hermite_xtol_u': HERMITE_ROOT_XTOL_U,
            'rtol': ROOT_RTOL, 'max_iterations': ROOT_MAXITER,
            'accepted_events': ('Negative-to-nonnegative radial sign brackets only; converged roots '
                                'must remain inside their bracket. Exclude both integration endpoints. '
                                'Root tolerances govern localization, not trajectory or rate accuracy.'),
        },
        'initial_elements': {'a_au': A_AU, 'eccentricity': E, 'a_m': A_M,
                             'source': 'NASA JPL Solar System Dynamics, J2000 Mercury mean elements, Table 1',
                             'url': 'https://ssd.jpl.nasa.gov/planets/approx_pos.html',
                             'use': 'initial-state inputs only; no trajectory or apsidal rate fitted to JPL data'},
        'scope': {
            'included': ['stationary point-Sun Newtonian monopole',
                         'leading nonspinning solar Schwarzschild correction in the enabled branch'],
            'excluded': ['planetary perturbations', 'solar oblateness', 'solar spin and frame dragging',
                         'finite Mercury mass and solar recoil', 'higher-order relativistic force terms',
                         'reference-frame precession'],
            'interpretation': 'Isolated test-particle difference, not the total observed motion or an ephemeris fit.',
        },
    }

    def checkpoint():
        args.output.parent.mkdir(parents=True, exist_ok=True)
        write_report(args.output, report)

    try:
        ok = run(report, checkpoint)
    except Exception as exc:
        report['status'] = 'blocked'
        report['error'] = {'type': type(exc).__name__, 'message': str(exc)}
        checkpoint()
        print(f'{type(exc).__name__}: {exc}', file=sys.stderr)
        return 2
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
