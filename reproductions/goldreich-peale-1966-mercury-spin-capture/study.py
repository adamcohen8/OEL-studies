"""Run from the directory containing spin_bridge.py: python study.py

Results go to study.json. Use --help for refinement and horizon controls.

Preregistered design: the six psi values include both endpoints. The primary
launch is at periapsis with spin/n = 1.525: the offset is small relative to
1.5 but exceeds the nominal resonant half-width. Change one assumption at
a time: use 1.522 at periapsis, or 1.525 at apoapsis. These are sensitivity
checks, not a probability distribution for the unspecified launch conditions.
No trajectory outcome was used to choose them.

Prior recognition: this is the classical spin-orbit separatrix energy-budget
problem (Goldreich-Peale type). Use the independently checked spin_bridge
for the fixed orbit and native conservative torque. No published probability
or six-trial result is used as input. The probability calculation below is
an averaged single-passage model, not a fit to the six trajectories.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from spin_bridge import SpinBridge


TAU = math.tau
RESONANCE = 1.5
LAUNCHES = ((0.0, 0.025), (0.0, 0.022), (math.pi, 0.025))
PSI = np.linspace(0, math.pi / 2, 6)
TABLE_SIZE = 2048
SPIN_LO, SPIN_HI = 1.468, 1.529  # certification band, not a clamp
CAPTURE_ORBITS = 256
ESCAPE_ORBITS = 32
EVENT_RESIDUAL_LIMIT = 2e-7


def orbit_row(bridge, M):
    n = bridge.mean_motion_rad_s
    sample = bridge.orbit_sample(M / n)
    scale = bridge.C_kg_m2 * n**2
    c = sample.torque_at_zero_nm / scale
    s = sample.torque_at_pi_over_4_nm / scale
    cm, sm = math.cos(3 * M), math.sin(3 * M)
    # q = theta - 3 M / 2; independent variable is n*t, i.e. unwrapped M.
    return (c * cm + s * sm, s * cm - c * sm,
            bridge.D * sample.a_over_r**6, sample.f_dot_rad_s / n)


def make_orbit(bridge):
    grid = np.linspace(0, TAU, TABLE_SIZE + 1)
    rows = np.array([orbit_row(bridge, M) for M in grid[:-1]])
    spline = CubicSpline(grid, np.vstack((rows, rows[0])),
                         bc_type='periodic', extrapolate='periodic')
    off_grid = TAU * (np.arange(257) + 0.37) / 257
    direct = np.array([orbit_row(bridge, M) for M in off_grid])
    delta = np.max(np.abs(spline(off_grid) - direct), axis=0)
    if np.any(delta > np.array([2e-11, 2e-11, 2e-12, 2e-9])):
        raise RuntimeError(f'periodic orbit cache check failed: {delta}')
    return spline, rows, delta


def averaged_tide(bridge):
    # Independently average the historical tide at constant inertial spin.
    # In true anomaly the measure (a/r)^6 dM is proportional to
    # (1+e*cos(f))^4 df. The negative-lag interval is |f| < f_switch.
    e = bridge.eccentricity
    b = 1 - e * e

    def weight(f):
        return (1 + e * math.cos(f))**4

    def integral(f):
        return quad(weight, 0, f, epsabs=1e-12, epsrel=1e-12)[0]

    total = integral(math.pi)
    scale = bridge.D / (math.pi * b**4.5)

    def tide(spin):
        z = (math.sqrt(spin * b**1.5) - 1) / e
        f_switch = math.acos(min(1.0, max(-1.0, z)))
        return -scale * (total - 2 * integral(f_switch))

    return tide, (1 - e)**2 / b**1.5, (1 + e)**2 / b**1.5


def prediction(bridge, coeff):
    a, b = coeff
    K = math.hypot(a, b)
    center = math.atan2(a, -b) / 2
    half_width = math.sqrt(2 * K)
    tide, fmin, fmax = averaged_tide(bridge)
    bias = tide(RESONANCE)
    if abs(bias) >= K:
        raise RuntimeError('no equilibrium in the averaged torque model')

    # Relative angle x is measured from the untilted center. Along the
    # conservative separatrix v = +/-sqrt(2K)*cos(x). Energy work is
    # integral T(v) dx, positive x on entry and negative x on return.
    # Uniform incident energy overshoot on [0,L_up] is a phase-mixing
    # hypothesis. The part of that interval which cannot leave on the
    # lower branch has length clip(L_up + L_down, 0, L_up).
    def branch(sign):
        breaks = []
        for edge in (fmin, fmax):
            ratio = (edge - RESONANCE) / (sign * half_width)
            if 0 < ratio < 1:
                breaks.append(math.acos(ratio))
        value, error = quad(
            lambda x: tide(RESONANCE + sign * half_width * math.cos(x)),
            0, math.pi / 2, points=breaks, epsabs=1e-14, epsrel=1e-10,
        )
        return -2 * sign * value, 2 * error

    up, up_error = branch(1)
    down, down_error = branch(-1)
    if up <= 0:
        raise RuntimeError('unexpected incident branch work sign')
    return dict(K=K, untilted_center_rad=center,
                tilted_center_rad=center + math.asin(bias / K) / 2,
                conservative_half_width=half_width,
                small_libration_period_orbits=1 / math.sqrt(2 * K),
                tide_at_resonance=bias, relative_bias=abs(bias) / K,
                lag_switch_spin_range=[fmin, fmax],
                energy_losses=[up, down], quadrature_error_estimates=[up_error, down_error],
                single_passage_probability=float(np.clip((up + down) / up, 0, 1)))


def acceleration(row, q, sign):
    a, b, d, _ = row
    return a * np.cos(2 * q) + b * np.sin(2 * q) - d * sign


def rk4(y, h, r0, rm, r1, sign):
    q, v = y
    a0 = acceleration(r0, q, sign)
    v2 = v + h / 2 * a0
    a2 = acceleration(rm, q + h / 2 * v, sign)
    v3 = v + h / 2 * a2
    a3 = acceleration(rm, q + h / 2 * v2, sign)
    v4 = v + h * a3
    a4 = acceleration(r1, q + h * v3, sign)
    return np.array((q + h / 6 * (v + 2 * (v2 + v3) + v4),
                     v + h / 6 * (a0 + 2 * (a2 + a3) + a4)))


def switch_step(orbit, t, h, initial, trial, r0, r1, sign):
    # Continue the incoming smooth branch for a provisional full step.
    # Hermite interpolate v on that branch to bracket/locate the switch;
    # redo RK4 to the root, then integrate the outgoing branch. The point
    # value of acceleration at exact equality has zero measure here.
    v0, v1 = initial[1], trial[1]
    m0 = h * acceleration(r0, initial[0], sign)
    m1 = h * acceleration(r1, trial[0], sign)

    def gap(u):
        v = v0 + u * (m0 + u * (3 * (v1 - v0) - 2 * m0 - m1
                                 + u * (2 * (v0 - v1) + m0 + m1)))
        return RESONANCE + v - orbit(t + h * u)[3]

    u = brentq(gap, 0, 1, xtol=2e-13)
    dt = u * h
    r_switch = orbit(t + dt)
    at_switch = rk4(initial, dt, r0, orbit(t + dt / 2), r_switch, sign)
    final = rk4(at_switch, h - dt, r_switch, orbit(t + (dt + h) / 2), r1, -sign)
    residual = abs(RESONANCE + at_switch[1] - r_switch[3])
    slope = abs(acceleration(r_switch, at_switch[0], sign) - orbit(t + dt, 1)[3])
    return final, float(residual), float(slope)


def crossing_certificate(bridge, steps):
    # Bounds in M time. Between consecutive apses f_dot/n is monotone.
    # In the only band where lag equality can occur, its slope beats both
    # smooth spin branches. Thus a step not straddling an apsis cannot
    # hide a pair of switches. Give endpoint speeds a full-step margin.
    e = bridge.eccentricity
    b = 1 - e**2
    low, high = (math.sqrt(w * b**1.5) for w in (SPIN_LO, SPIN_HI))
    z = max(abs((g - 1) / e) for g in (low, high))
    f_slope_lower = 2 * e * math.sqrt(1 - z**2) * low**3 / b**3
    spin_slope_upper = (1.5 * bridge.triaxiality / (1 - e)**3
                        + bridge.D / (1 - e)**6)
    if f_slope_lower < 10 * spin_slope_upper:
        raise RuntimeError('transversality margin too small')
    margin = TAU / steps * spin_slope_upper
    return dict(f_slope_lower=f_slope_lower, spin_slope_upper=spin_slope_upper,
                step_start_spin_range=[SPIN_LO + margin, SPIN_HI - margin])


def run_group(bridge, orbit, averaged, steps, horizon, M0, offset):
    h = TAU / steps
    cert = crossing_certificate(bridge, steps)
    lo, hi = cert['step_start_spin_range']
    times = np.linspace(M0, M0 + TAU, steps + 1)
    nodes = orbit(times)
    middle = orbit((times[:-1] + times[1:]) / 2)
    f0 = bridge.orbit_sample(M0 / bridge.mean_motion_rad_s).true_anomaly_rad
    y = np.array([PSI + f0 - RESONANCE * M0, np.full(PSI.size, offset)])
    active = np.ones(PSI.size, dtype=bool)
    count = np.zeros(PSI.size, dtype=int)
    residual = np.zeros(PSI.size)
    slope = np.full(PSI.size, np.inf)
    dwell = np.zeros(PSI.size, dtype=int)
    last_well = np.full(PSI.size, np.nan)
    first_below = [None] * PSI.size
    status = ['unresolved'] * PSI.size
    reason = ['horizon'] * PSI.size
    stop = np.zeros(PSI.size)
    trace = [[[0, float(q), float(v)]] for q, v in y.T]
    history = np.zeros((ESCAPE_ORBITS + 1, PSI.size))
    history[0] = y[0]
    center = averaged['tilted_center_rad']
    untilted = averaged['untilted_center_rad']
    K = averaged['K']

    def finish(i, why, age, label='unresolved'):
        active[i] = False
        status[i], reason[i], stop[i] = label, why, age
        row = [float(age), *map(float, y[:, i])]
        if trace[i][-1] != row:
            trace[i].append(row)

    for k in range(1, horizon + 1):
        start_count = count.copy()
        ids = np.flatnonzero(active)
        for j in range(steps):
            if not ids.size:
                break
            age = k - 1 + j / steps
            spin = RESONANCE + y[1, ids]
            bad = (spin < lo) | (spin > hi) | ~np.isfinite(spin)
            if np.any(bad):
                for i in ids[bad]:
                    finish(i, 'left_switch_certification_band', age)
                ids = ids[~bad]
                if not ids.size:
                    break
            r0, rm, r1 = nodes[j], middle[j], nodes[j + 1]
            old = y[:, ids]
            signs = np.where(RESONANCE + old[1] >= r0[3], 1.0, -1.0)
            trial = rk4(old, h, r0, rm, r1, signs)
            crossed = np.flatnonzero(signs * (RESONANCE + trial[1] - r1[3]) < 0)
            failed = []
            for p in crossed:
                i = ids[p]
                trial[:, p], err, derivative = switch_step(
                    orbit, times[j], h, old[:, p], trial[:, p], r0, r1, signs[p])
                count[i] += 1
                residual[i] = max(residual[i], err)
                slope[i] = min(slope[i], derivative)
                if (err > EVENT_RESIDUAL_LIMIT or derivative < cert['f_slope_lower'] / 2
                        or signs[p] * (RESONANCE + trial[1, p] - r1[3]) >= 0):
                    failed.append(i)
            y[:, ids] = trial
            if failed:
                for i in failed:
                    finish(i, 'switch_audit_failed', age + 1 / steps)
                ids = np.flatnonzero(active)

        history[k % len(history)] = y[0]
        for i in np.flatnonzero(active):
            q, v = y[:, i]
            if count[i] - start_count[i] != 2:
                finish(i, 'expected_two_switches_per_completed_orbit', k)
                continue
            if v < 0 and first_below[i] is None:
                first_below[i] = k
            if k % 32 == 0:
                trace[i].append([k, float(q), float(v)])
            # Strobe in one fixed well, away from the conservative saddle.
            # Deep sustained libration, not a temporary separatrix pause.
            well = round((q - center) / math.pi)
            energy = v*v / 2 + K * math.sin(q - untilted)**2
            if abs(q - center - well * math.pi) < 1.10 and energy < 0.85 * K:
                dwell[i] = dwell[i] + 1 if well == last_well[i] else 1
                last_well[i] = well
            else:
                dwell[i] = 0
            if dwell[i] >= CAPTURE_ORBITS:
                finish(i, 'deep_same_well_dwell', k, 'capture')
            elif (k >= ESCAPE_ORBITS and v < -0.023
                  and q - history[(k - ESCAPE_ORBITS) % len(history), i] < -TAU):
                finish(i, 'subresonant_and_backward_winding', k, 'escape')
        if not np.any(active):
            break

    for i in np.flatnonzero(active):
        finish(i, 'horizon', horizon)
    return dict(initial_M_rad=M0, initial_spin_over_n=RESONANCE + offset,
                steps_per_orbit=steps, certificate=cert,
                trials=[dict(psi_rad=float(PSI[i]), classification=status[i], reason=reason[i],
                             stop_orbits=float(stop[i]), first_strobe_below_3_2=first_below[i],
                             switches=int(count[i]), maximum_switch_gap=float(residual[i]),
                             minimum_switch_slope=float(slope[i]) if count[i] else None,
                             final_dwell_orbits=int(dwell[i]), trace=trace[i])
                        for i in range(PSI.size)])


def compare(coarse, fine):
    trials = []
    for a, b in zip(coarse['trials'], fine['trials'], strict=True):
        label = a['classification'] if a['classification'] == b['classification'] else 'unresolved'
        complete = a['first_strobe_below_3_2'], b['first_strobe_below_3_2']
        trials.append(dict(psi_rad=a['psi_rad'], classification=label,
                           first_below_difference_orbits=(complete[1] - complete[0]
                                                          if None not in complete else None),
                           stop_difference_orbits=b['stop_orbits'] - a['stop_orbits']))
    counts = {key: sum(t['classification'] == key for t in trials)
              for key in ('capture', 'escape', 'unresolved')}
    n = len(trials)
    return dict(initial_M_rad=fine['initial_M_rad'],
                initial_spin_over_n=fine['initial_spin_over_n'], trials=trials, counts=counts,
                uniform_discrete_psi_capture_fraction_bounds=[counts['capture'] / n,
                    (counts['capture'] + counts['unresolved']) / n])


def main():
    cli = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    cli.add_argument('--output', type=Path, default=Path('study.json'))
    cli.add_argument('--steps', type=int, default=64, help='coarse steps/orbit; also runs at twice this')
    cli.add_argument('--orbits', type=int, default=4000, help='maximum orbits per trajectory')
    args = cli.parse_args()
    if args.steps < 32 or args.steps % 4:
        cli.error('--steps must be a multiple of four and at least 32')
    if args.orbits < 1:
        cli.error('--orbits must be positive')

    bridge = SpinBridge()
    orbit, rows, delta = make_orbit(bridge)
    orbit_audit = dict(maximum_absolute_off_grid_error=delta)
    # Uniform mean-anomaly projection of the two native torque columns.
    averaged = prediction(bridge, rows[:, :2].mean(axis=0))

    settings = LAUNCHES
    experiments = []
    for M0, offset in settings:
        coarse = run_group(bridge, orbit, averaged, args.steps, args.orbits, M0, offset)
        fine = run_group(bridge, orbit, averaged, 2 * args.steps, args.orbits, M0, offset)
        experiments.append(dict(comparison=compare(coarse, fine), coarse=coarse, fine=fine))

    result = dict(
        model=dict(eccentricity=bridge.eccentricity, triaxiality=bridge.triaxiality,
                   D=bridge.D, orbit='fixed planar prograde Kepler ellipse',
                   time_unit='1/n', phase_origin='periapsis',
                   lag_switch='theta_dot = f_dot'),
        settings=dict(horizon_orbits=args.orbits, steps_per_orbit=[args.steps, 2 * args.steps],
                      capture_dwell_orbits=CAPTURE_ORBITS, escape_window_orbits=ESCAPE_ORBITS,
                      event_residual_limit=EVENT_RESIDUAL_LIMIT,
                      psi_rad=PSI.tolist(), primary_experiment_index=0,
                      initial_M_and_spin_offset=settings),
        interpretation=dict(
            prior='Goldreich–Peale spin-orbit separatrix-flux framework; use the supplied lag-sign law, not a substituted tidal law.',
            single_passage='Averaged leading-order separatrix estimate with uniform incoming energy phase in the narrow flux interval. Not an exact probability for the unaveraged dynamics or for six equally spaced initial angles.',
            unknown_initial_conditions='The source gives neither the exact positive spin offset nor the starting anomaly. The primary choice and the two sensitivity choices are assumptions.',
            six_state_measure='At each fixed epoch and offset, equal mass 1/6 at each specified psi. The reported fraction interval assigns unresolved states first all to escape and then all to capture. It is not a binomial confidence interval.',
            convergence='Agreement under step halving is a numerical check, not an error bound; disagreements and unfinished or failed audits are unresolved. No pooling across epochs or spin offsets.',
            capture='Operational finite-time classification from deep same-well stroboscopic dwell, not a proof of infinite-time trapping.',
            escape='Subresonant spin and sustained net backward winding through more than two resonant wells.',
            trace_columns=['orbits_since_start', 'unwrapped_theta_minus_1.5_M_rad', 'spin_over_n_minus_1.5']),
        orbit_cache_audit=orbit_audit,
        single_passage_averaged=averaged,
        experiments=experiments,
    )

    def numpy_json(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f'cannot encode {type(value).__name__}')

    args.output.write_text(json.dumps(result, ensure_ascii=False, allow_nan=False,
                                      separators=(',', ':'), default=numpy_json) + '\n',
                           encoding='utf-8')


if __name__ == '__main__':
    main()
