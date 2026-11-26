"""Source-informed follow-up; run `python paper_targets.py` from the workspace.

The Wisdom, Peale & Mignard Figure 2 comparison supplied the three targets.
This is not an independent discovery or a revision of the sealed blind search.
All trajectories use the source-locked OEL periapse map in replay.py.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

import replay as r


STEPS = (128, 256)
TARGETS = (('9/4', 9), ('3', 12), ('7/2', 14))  # pi of lift in TWO orbits
CLOSE = 1e-8
FIT_TOL = 1e-10
PERIOD_GAP = 1e-5
MATCH = 1e-3
DUP = 2e-6
FD = (1e-5, 5e-6)
SCAN_OFFSETS = (-0.48, -0.32, -0.16, 0., 0.16, 0.32, 0.48)
MAX_BRACKETS_PER_LINE = 2
NEAR_ORBITS = 16


class BudgetExpired(Exception):
    pass


def distance(x, y):
    d = np.asarray(x) - np.asarray(y)
    return float(np.hypot((d[0] + np.pi / 2) % np.pi - np.pi / 2, d[1]))


def section(x):
    return [float(x[0] % np.pi), float(x[1])]


class NativeMap:
    def __init__(self, engine, deadline):
        self.engine = engine
        self.deadline = deadline

    def flow(self, x, p):
        if time.process_time() >= self.deadline:
            raise BudgetExpired
        y, lift = r.flow(self.engine, *x, p)
        return np.asarray(y, dtype=float), float(lift)

    def residual(self, x, p, m):
        y, lift = self.flow(x, p)
        return np.array([np.pi * (lift - m), y[1] - x[1]])


def correct(native, seed, m):
    seed = np.asarray(seed, dtype=float)
    fun = lambda x: native.residual(x, 2, m)
    initial = fun(seed)
    x, err, status, nfev = seed, initial, None, 0
    if np.linalg.norm(initial) > FIT_TOL:
        span = np.array([0.24, 0.10])
        fit = least_squares(fun, seed, bounds=(seed - span, seed + span),
                            max_nfev=8, ftol=1e-12, xtol=1e-12, gtol=1e-13)
        x, err, status, nfev = fit.x, fit.fun, int(fit.status), int(fit.nfev)
    return x, dict(theta_w=section(x), seed_residual=initial.tolist(),
                   residual=err.tolist(), solver_status=status, nfev=nfev,
                   closed=bool(np.linalg.norm(err) < CLOSE))


def examine(native, x, m):
    y, lift1 = native.flow(x, 1)
    gap = distance(y, x)
    short_lift_error = abs(np.pi * (lift1 - m / 2))
    if gap < CLOSE and short_lift_error < CLOSE and m % 2 == 0:
        period = 1
    elif gap > PERIOD_GAP:
        period = 2
    else:
        period = None
    samples = []
    if period is not None:
        for h in FD:
            cols = []
            for dx in h * np.eye(2):
                plus, lp = native.flow(x + dx, period)
                minus, lm = native.flow(x - dx, period)
                # Differentiate the lifted endpoint, not its wrapped chart.
                cols.append([2 * dx[0] + np.pi * (lp - lm), plus[1] - minus[1]])
            jac = np.array(cols).T / (2 * h)
            eig = np.linalg.eigvals(jac)
            samples.append(dict(h=h, trace=float(np.trace(jac)), determinant=float(np.linalg.det(jac)),
                                multipliers_re_im=[[float(v.real), float(v.imag)] for v in eig]))
    return dict(fundamental_period=period, first_return_gap=gap, first_return_lift_pi=lift1,
                first_return_half_target_error_rad=short_lift_error,
                cycle_theta_w=[section(x), section(y)][:period], linear_maps=samples)


def classify(a, b):
    if a['fundamental_period'] is None or a['fundamental_period'] != b['fundamental_period']:
        return 'period_unresolved'
    maps = a['linear_maps'] + b['linear_maps']
    tr = np.array([v['trace'] for v in maps])
    det_error = max(abs(v['determinant'] - 1) for v in maps)
    spread = float(np.ptp(tr))
    if det_error > 1e-3 or spread > 0.01 * max(1., abs(tr[-1])):
        return 'linear_type_unresolved'
    margin = max(1e-3, 4 * spread, 4 * det_error)
    if abs(tr[-1]) < 2 - margin:
        return 'elliptic_linearization'
    if abs(tr[-1]) > 2 + margin:
        return 'hyperbolic_linearization'
    return 'near_marginal_or_unresolved'


def neighborhood(native, cycle, row):
    row['neighborhood_columns'] = ['orbit', 'theta_mod_pi_rad', 'omega_over_n', 'distance_to_center_at_same_phase']
    row['neighborhood'] = []
    for delta in ((0.002, 0.), (0., 0.002)):
        x = np.asarray(cycle[0]) + delta
        samples = []
        path = dict(initial_offset=list(delta), samples=samples, completed=False)
        row['neighborhood'].append(path)
        for k in range(NEAR_ORBITS + 1):
            gap = distance(x, cycle[k % 2])
            samples.append([k, *section(x), gap])
            path['max_phase_matched_distance'] = max(s[3] for s in samples)
            if k < NEAR_ORBITS:
                x, _ = native.flow(x, 1)
        path['completed'] = True


def trial(coarse, fine, target, seed, kind, extra=None):
    m = target['target_lift_pi_two_orbits']
    rows = target['trials']
    row = dict(seed_theta_w=section(seed), method=kind, status='started')
    if extra is not None:
        row.update(extra)
    rows.append(row)
    x, row['coarse_fit'] = correct(coarse, seed, m)
    if not row['coarse_fit']['closed']:
        row['status'] = 'coarse_closure_not_found'
        return
    # Both section points of a two-cycle represent the same orbit.
    for i, old in enumerate(rows[:-1]):
        if old.get('status') not in ('matched_two_cycle', 'short_return_rejected'):
            continue
        if min(distance(x, point) for point in old['coarse_cycle']['cycle_theta_w']) < DUP:
            row.update(status='duplicate', duplicate_of_trial=i)
            return
    y, row['fine_fit'] = correct(fine, x, m)
    row['root_shift'] = distance(x, y)
    if not row['fine_fit']['closed'] or row['root_shift'] > MATCH:
        row['status'] = 'step_halved_root_not_matched'
        return
    row['coarse_at_fine_residual'] = coarse.residual(y, 2, m).tolist()
    a = row['coarse_cycle'] = examine(coarse, x, m)
    b = row['fine_cycle'] = examine(fine, y, m)
    row['classification'] = classify(a, b)
    pa, pb = a['fundamental_period'], b['fundamental_period']
    if pa == pb == 1:
        row['status'] = 'short_return_rejected'
    elif pa == pb == 2:
        row['status'] = 'matched_two_cycle'
    else:
        row['status'] = 'minimal_period_unresolved'


def symmetry_search(coarse, fine, target):
    # At a periapse the reversible section involution is (theta,w)->(-theta,w).
    # An exact step from one fixed line to either fixed line implies a two-return.
    # We use that only to seed an unconstrained two-coordinate closure test.
    m = target['target_lift_pi_two_orbits']
    wind = m / 4
    for theta in (0., np.pi / 2):
        scan = dict(theta_rad=float(theta), samples_w_lift_error_pi=[], brackets=[], shots=[])
        target['scans'].append(scan)
        cache = {}

        def shot(w):
            w = float(w)
            if w not in cache:
                _, lift = coarse.flow((theta, w), 1)
                cache[w] = lift - m / 2
            return cache[w]

        grid = [wind + d for d in SCAN_OFFSETS]
        for w in grid:
            scan['samples_w_lift_error_pi'].append([w, shot(w)])
        brackets = [(w, w) for w in grid if cache[w] == 0]
        brackets += [(a, b) for a, b in zip(grid, grid[1:]) if cache[a] * cache[b] < 0]
        brackets.sort(key=lambda ab: abs(sum(ab) / 2 - wind))
        scan['brackets'] = brackets
        for a, b in brackets[:MAX_BRACKETS_PER_LINE]:
            rec = dict(bracket=[a, b])
            scan['shots'].append(rec)
            try:
                w = a if a == b else brentq(shot, a, b, xtol=1e-10, maxiter=16)
            except (RuntimeError, ValueError) as exc:
                rec['failure'] = type(exc).__name__
                continue
            rec['w'] = float(w)
            rec['lift_error_pi'] = shot(w)
            trial(coarse, fine, target, (theta, w), 'symmetry_shot', dict(bracket=[a, b]))


def add_neighborhoods(fine, target):
    for row in target['trials']:
        if (row['status'] == 'matched_two_cycle' and row['classification'] == 'elliptic_linearization'
                and 'neighborhood' not in row):
            neighborhood(fine, row['fine_cycle']['cycle_theta_w'], row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lock', type=Path, default=Path('oel-revision.lock.json'))
    parser.add_argument('--out', type=Path, default=Path('outputs/paper_targets.json'))
    parser.add_argument('--cpu-seconds', type=float, default=165., help='soft process CPU limit, checked before each map call')
    args = parser.parse_args()
    if not np.isfinite(args.cpu_seconds) or args.cpu_seconds <= 0:
        parser.error('--cpu-seconds must be positive and finite')
    start = time.process_time()
    api, source = r.replay_api(args.lock)
    r.CHUNKS = 16
    engines, calibration = zip(*(r.checks(api, r.replay_config(4, r.T / n, 'rkf78')) for n in STEPS))
    coarse, fine = [NativeMap(e, start + args.cpu_seconds) for e in engines]
    targets = [dict(winding_turns_per_orbit=label, target_lift_pi_two_orbits=m, scans=[], trials=[])
               for label, m in TARGETS]
    report = dict(study='source_informed_followup_not_independent_discovery', source=source,
                  prior_recognition='The supplied Figure 2 comparison identifies 9/4, 3 and 7/2 as targets.',
                  evidence='Native OEL Saturn point-mass orbit and passive attitude, via validated locked replay.',
                  orbit_integrator='rkf78', chunks=r.CHUNKS, attitude_steps_per_orbit=STEPS,
                  calibration=calibration, closure_tol=CLOSE, period_gap=PERIOD_GAP, match_tol=MATCH,
                  unique_cycle_tol=DUP, max_brackets_per_line=MAX_BRACKETS_PER_LINE,
                  cpu_limit_s=args.cpu_seconds, targets=targets,
                  section='Nominal periapse; theta from +X modulo pi, radians; w = body z spin / mean motion.',
                  interpretation=[
                      'The labels are net 2pi turns per orbit, not the instantaneous w at periapse.',
                      'Symmetry is a seed heuristic only. Roots require both coordinates and the specified lift.',
                      'One-orbit closures are rejected as new two-cycles; their linear maps are one-orbit maps.',
                      'Stability is local. Determinant near one is a consistency check, not a classification by itself.',
                      'Neighborhoods iterate the fixed-orbit one-return map at fine resolution; each call restarts the same Kepler phase. The samples test only finite-time behavior.',
                      'A finite search, including failed, unmatched or unattempted trials, cannot exclude other cycles.'
                  ])
    phase = 'symmetry_scans'
    try:
        for target in targets:
            symmetry_search(coarse, fine, target)
        island = targets[0]
        phase = 'first_9_over_4_neighborhoods'
        add_neighborhoods(fine, island)
        phase = 'off_axis_9_over_4_trials'
        shots = [shot['w'] for scan in island['scans'] for shot in scan['shots'] if 'w' in shot]
        w = float(np.median(shots)) if shots else 9 / 4
        for theta in (np.pi / 4, -np.pi / 4):
            trial(coarse, fine, island, (theta, w), 'off_axis')
        phase = 'remaining_9_over_4_neighborhoods'
        add_neighborhoods(fine, island)
        report['termination'] = 'bounded_design_completed'
    except BudgetExpired:
        report['termination'] = 'cpu_budget_reached'
        report['interrupted_phase'] = phase
    report['native_step_calls'] = [int(e.calls) for e in engines]
    report['process_cpu_s'] = round(time.process_time() - start, 3)
    text = json.dumps(report, separators=(',', ':'), allow_nan=False)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text + '\n', encoding='utf-8')
    print(text)


if __name__ == '__main__':
    main()
