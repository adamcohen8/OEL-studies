"""From the workspace: python multireturn.py

Bounded searches of the nominal periapse map. The solve uses both coordinates
and a specified lift, without imposing reflection symmetry. All trajectories
and finite-difference maps are from the source-locked native replay.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

import replay as r

CLOSE = 1e-8
PERIOD_GAP = 1e-5
FD = (1e-5, 5e-6)


def residual(eng, x, p, m):
    y, lift_pi = r.flow(eng, *x, p)
    return np.array([np.pi * (lift_pi - m), y[1] - x[1]])


def solve(eng, start, center, p, m):
    fun = lambda x: residual(eng, x, p, m)
    span = np.array([0.4, 0.30])
    fit = least_squares(fun, start, bounds=(center - span, center + span),
                        max_nfev=20, ftol=1e-11, xtol=1e-11, gtol=1e-11)
    return fit.x, {'x': fit.x.tolist(), 'residual': fit.fun.tolist(),
                   'closed': bool(np.linalg.norm(fit.fun) < CLOSE),
                   'solver_status': int(fit.status), 'nfev': fit.nfev}


def examine(eng, x, p):
    data, missing = r.track(eng, x, p, 1, False)
    peri = data['peri']
    gaps = r.section_delta(peri[1:], peri[0])
    divisors = [k for k in range(1, p) if p % k == 0]
    period = next((k for k in divisors if gaps[k - 1] < CLOSE), p)
    if any(CLOSE <= gaps[k - 1] <= PERIOD_GAP for k in divisors if k < period):
        period = None
    maps = []
    if period is not None:
        for h in FD:
            cols = []
            for dx in h * np.eye(2):
                plus, lp = r.flow(eng, *(x + dx), period)
                minus, lm = r.flow(eng, *(x - dx), period)
                cols.append([2 * dx[0] + np.pi * (lp - lm), plus[1] - minus[1]])
            jac = np.array(cols).T / (2 * h)
            maps.append([float(np.trace(jac)), float(np.linalg.det(jac))])
    return {'cycle_theta_w': peri[:-1, :2].tolist(), 'return_distances_1_up': gaps.tolist(),
            'turns_per_orbit': float(peri[-1, 2] / p), 'fundamental_period': period,
            'trace_det': maps, 'unbracketed_events': missing,
            'event_columns': r.EVENT_COLUMNS, 'events': data['events'].tolist()}


def compare(a, b, drift):
    if drift > 1e-3 or a['fundamental_period'] != b['fundamental_period']:
        return 'root_match_unresolved'
    maps = np.array(a['trace_det'] + b['trace_det'])
    if len(maps) != 4:
        return 'fundamental_period_unresolved'
    tr, det = maps.T
    spread = np.ptp(tr) / max(1., abs(tr[-1]))
    if max(abs(det - 1)) > 1e-3 or spread > 0.01:
        return 'linear_type_unresolved'
    margin = max(1e-3, 4 * np.ptp(tr), 4 * max(abs(det - 1)))
    if abs(tr[-1]) < 2 - margin:
        return 'elliptic_linearization'
    if abs(tr[-1]) > 2 + margin:
        return 'hyperbolic_linearization'
    return 'near_marginal_or_unresolved'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=Path('outputs/multireturn.json'))
    args = parser.parse_args()
    api, source = r.replay_api(Path('oel-revision.lock.json'))
    r.CHUNKS = 16
    steps = [128, 256]
    (coarse, cc), (fine, fc) = [r.checks(api, r.replay_config(4, r.T / n, 'rkf78')) for n in steps]
    rows = []
    for p, m, w in [(3, 2, 0.), (4, 8, 0.5)]:
        center = np.array([0., w])
        seen = []
        for offset in [(0, 0), (0.14, 0.08), (-0.14, 0.08), (0.09, -0.10), (-0.09, -0.10)]:
            start = center + offset
            x, fit = solve(coarse, start, center, p, m)
            row = dict(p=p, target_pi_winding=m, start=start.tolist(), coarse_fit=fit)
            rows.append(row)
            if not fit['closed']:
                continue
            duplicate = next((i for i, c in seen if np.min(r.section_delta(c, x)) < 1e-6), None)
            if duplicate is not None:
                row['duplicate_of_row_zero_based'] = duplicate
                continue
            a = row['coarse_cycle'] = examine(coarse, x, p)
            y, row['fine_fit'] = solve(fine, start, center, p, m)
            if not row['fine_fit']['closed'] or r.section_delta(x, y) > 1e-3:
                y, row['fine_fit_from_coarse'] = solve(fine, x, center, p, m)
            row['fine_at_coarse_residual'] = residual(fine, x, p, m).tolist()
            if not row.get('fine_fit_from_coarse', row['fine_fit'])['closed']:
                continue
            b = row['fine_cycle'] = examine(fine, y, p)
            row['coarse_at_fine_residual'] = residual(coarse, y, p, m).tolist()
            row['root_shift'] = float(r.section_delta(x, y))
            row['classification'] = compare(a, b, row['root_shift'])
            if row['root_shift'] <= 1e-3:
                seen.append((len(rows) - 1, np.array(a['cycle_theta_w'])))
    report = dict(source=source, chunks=r.CHUNKS, orbit_integrator='rkf78',
                  attitude_steps_per_orbit=steps, calibration=[cc, fc], trials=rows,
                  native_step_calls=[coarse.calls, fine.calls], closure_tol=CLOSE,
                  interpretation='Winding is in 2pi turns per orbit. A shorter closure is not a new cycle. '
                  'Linear type is local, not a claim of long-time stability. Failed or mismatched '
                  'solves are inconclusive: finite starts, bounds, lift and work limits do not exclude roots.')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
