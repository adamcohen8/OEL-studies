"""Reconstruct exterior attraction-error contours and scan their radial extent."""

import argparse
import csv
import math
import platform
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from sim.dynamics.orbit.accelerations import accel_two_body

from cube_bridge import cube_accel


LEVELS = (0.1, 1.0, 10.0, 40.0)
METRICS = ('magnitude_cube', 'vector_cube', 'magnitude_mascon')
RV = (3.0 / (4.0 * math.pi)) ** (1.0 / 3.0)


def fields(x, y):
    if max(abs(x), abs(y)) < 0.5:
        raise ValueError('query must be on or outside the cube')
    p = np.array([x, y, 0.0])
    return cube_accel(p, 1.0, 1.0), np.asarray(accel_two_body(p, 1.0))


def measures(c, m):
    cn, mn = math.hypot(*c), math.hypot(*m)
    diff = math.hypot(*(m - c))
    return cn, mn, 100.0 * (mn - cn) / cn, 100.0 * diff / cn, 100.0 * abs(mn - cn) / mn


def errors(x, y):
    _, _, signed, vector, other = measures(*fields(x, y))
    return abs(signed), vector, other


def write_csv(path, header, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def controls():
    for x, y in [(0.5, 0.0), (0.5, 0.5), (0.74, 0.17), (1.8, -1.2), (4.0, 0.0)]:
        c, m = fields(x, y)
        p = np.array([x, y, 0.0])
        np.testing.assert_allclose(m, -p / math.hypot(x, y)**3, rtol=2e-14, atol=1e-15)
        reflected_c, reflected_m = fields(-y, x)
        for a, b in ((c, reflected_c), (m, reflected_m)):
            np.testing.assert_allclose(b, [-a[1], a[0], a[2]], rtol=0.0, atol=2e-14)
            if not np.dot(a, p) < 0.0:
                raise AssertionError('expected inward radial projection')
    print('Controls passed: native inverse-square values, rotational symmetry, inward projection.', flush=True)


def near_table(out):
    cases = [('face_normal', 0.5 + d, 0.0) for d in
             (0.0, 1e-8, 1e-5, 0.001, 0.01, 0.025, 0.05, 0.1, RV-0.5, 0.15, 0.25, 0.5)]
    cases += [('face_surface', 0.5, y) for y in (0.1, 0.2, 0.3, 0.4, 0.49, 0.5)]
    cases += [('edge_diagonal', 0.5+d, 0.5+d) for d in (0.001, 0.01, 0.05, 0.1, 0.25)]
    rows = []
    for label, x, y in cases:
        c, m = fields(x, y)
        cn, mn, signed, vector, other = measures(c, m)
        r = math.hypot(x, y)
        angle = math.degrees(math.atan2(abs(c[0]*m[1]-c[1]*m[0]), np.dot(c, m)))
        rows.append((label, x, y, r, r >= RV, cn, mn, signed, vector, other, angle))
    write_csv(out/'near.csv', ('track', 'x', 'y', 'r', 'equal_volume_exterior_or_limit',
              'cube_mag', 'mascon_mag', 'signed_magnitude_error_pct', 'vector_error_pct',
              'mascon_normalized_magnitude_error_pct', 'direction_delta_deg'), rows)
    print(f'equal_volume_radius={RV:.12g}; saved {len(rows)} near-field rows', flush=True)


def contours(out, intervals, save_plot):
    half = intervals // 2
    axis = np.linspace(-2.0, 2.0, intervals + 1)
    quadrant = np.full((half+1, half+1), np.nan)
    for i in range(half+1):
        x = axis[half+i]
        if x < 0.5:
            continue
        for j in range(i+1):
            z = measures(*fields(x, axis[half+j]))[2]
            quadrant[i, j] = quadrant[j, i] = z
    indices = np.abs(np.arange(intervals+1)-half)
    grid = np.ma.masked_invalid(quadrant[np.ix_(indices, indices)])
    fig, ax = plt.subplots(figsize=(7, 6))
    palette = dict(zip(LEVELS, ('#313695', '#4575b4', '#d95f02', '#a50026')))
    # Interpolate the signed field before taking the union of the +/- levels.
    # Sampling its absolute value first can lose both crossings in one cell.
    signed_levels = [-v for v in reversed(LEVELS)] + list(LEVELS)
    cs = ax.contour(axis, axis, grid, levels=signed_levels, corner_mask=False,
                    colors=[palette[abs(v)] for v in signed_levels], linestyles='solid')
    rows = []
    segment_counts = dict.fromkeys(LEVELS, 0)
    for signed_level, segments in zip(cs.levels, cs.allsegs):
        level = abs(signed_level)
        radii, residuals = [], []
        nonempty = [seg for seg in segments if len(seg)]
        panel_ends, cube_ends = 0, 0
        for segment in nonempty:
            sid = segment_counts[level]
            segment_counts[level] += 1
            if not np.array_equal(segment[0], segment[-1]):
                for x, y in (segment[0], segment[-1]):
                    extent = max(abs(x), abs(y))
                    panel_ends += abs(extent - 2.0) <= 1e-12
                    cube_ends += abs(extent - 0.5) <= 1e-12
            for k, (x, y) in enumerate(segment):
                rows.append((level, sid, k, x, y, signed_level))
                radii.append(math.hypot(x, y))
                residuals.append(abs(measures(*fields(x, y))[2] - signed_level))
        print(f'contour n={intervals} signed_level_pct={signed_level:g} segments={len(nonempty)} '
              f'vertices={len(radii)} r_min={min(radii, default=math.nan):.9g} '
              f'r_max={max(radii, default=math.nan):.9g} '
              f'max_vertex_residual_pp={max(residuals, default=math.nan):.6g} '
              f'panel_ends={panel_ends} cube_ends={cube_ends}', flush=True)
    if not all(segment_counts.values()):
        raise AssertionError(f'missing absolute contour: {segment_counts}')
    write_csv(out/f'contours_{intervals}.csv',
              ('level_pct', 'segment', 'vertex', 'x', 'y', 'signed_level_pct'), rows)
    if save_plot:
        handles = [Line2D([0], [0], color=palette[v], label=f'{v:g}%') for v in LEVELS]
        ax.legend(handles=handles, title='Absolute error', loc='upper left',
                  bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
        ax.add_patch(plt.Rectangle((-0.5, -0.5), 1, 1, facecolor='0.88', edgecolor='black'))
        ax.set(xlim=(-2, 2), ylim=(-2, 2), aspect='equal', xlabel='x / side', ylabel='y / side',
               title='Absolute attraction-magnitude error (cube denominator)')
        fig.tight_layout()
        fig.savefig(out/'contours.png', dpi=180)
    plt.close(fig)


def ray_errors(r, c, s):
    return errors(max(0.5, r*c), r*s)


def refine(lo, hi, level, metric, c, s):
    # lo is not below level and hi is below. Equality is assigned to lo.
    for _ in range(32):
        mid = (lo + hi) * 0.5
        if ray_errors(mid, c, s)[metric] < level:
            hi = mid
        else:
            lo = mid
    return (lo + hi) * 0.5


def scan(out, na, nr):
    rows = []
    grouped = {(m, v): [] for m in METRICS for v in LEVELS}
    for theta in np.linspace(0.0, math.pi/4, na+1):
        c, s = math.cos(theta), math.sin(theta)
        rr = np.linspace(0.5/c, 4.0, nr+1)
        data = np.array([ray_errors(r, c, s) for r in rr])
        for j, metric in enumerate(METRICS):
            for level in LEVELS:
                high = data[:, j] >= level
                down = np.flatnonzero(high[:-1] & ~high[1:])
                up = np.flatnonzero(~high[:-1] & high[1:])
                root = None
                if len(down):
                    k = down[-1]
                    root = refine(rr[k], rr[k+1], level, j, c, s)
                row = (math.degrees(theta), metric, level, root, len(down), len(up),
                       bool(high[-1]), data[0, j], data[-1, j])
                rows.append(row)
                grouped[(metric, level)].append(row)
    write_csv(out/f'rays_{na}.csv', ('theta_deg', 'metric', 'level_pct', 'last_down_r',
              'down_steps', 'up_steps', 'unresolved_at_r4', 'boundary_pct', 'r4_pct'), rows)
    summary = {}
    for (metric, level), group in grouped.items():
        crossings = [row for row in group if row[3] is not None]
        last = max(crossings, key=lambda row: row[3]) if crossings else None
        earliest = min((row[3] for row in crossings), default=None)
        unresolved = sum(row[6] for row in group)
        summary[(metric, level)] = None if last is None else last[3]
        print(f'rays na={na} nr={nr} metric={metric} level_pct={level:g} '
              f'min_last_down_r={earliest} max_last_down_r={summary[(metric, level)]} '
              f'theta_at_max_deg={None if last is None else last[0]} '
              f'no_down={len(group)-len(crossings)} rays_with_up={sum(row[5] > 0 for row in group)} '
              f'max_transitions={max(row[4]+row[5] for row in group)} unresolved={unresolved}', flush=True)
    return summary


def shells(out, threshold_r):
    theta = np.linspace(0.0, math.pi/4, 721)
    xy = [(math.cos(t), math.sin(t)) for t in theta]
    radii = [0.75, 1, 1.5, 2, 2.5, 3, 3.25, 3.5, 4]
    if threshold_r is not None:
        radii.extend(r for r in (threshold_r-0.01, threshold_r, threshold_r+0.01)
                     if math.sqrt(0.5) <= r <= 4)
    rows = []
    for r in sorted(set(radii)):
        values = np.array([ray_errors(r, c, s) for c, s in xy])
        for j, metric in enumerate(METRICS):
            k = np.argmax(values[:, j])
            rows.append((r, metric, float(values[:, j].min()), values[k, j], math.degrees(theta[k])))
    write_csv(out/'shells.csv', ('r', 'metric', 'min_pct', 'max_pct', 'theta_at_max_deg'), rows)


def tail_bound(r):
    if r < 4.0:
        raise ValueError('tail check is specified for r >= 4')
    moment6 = 3/448 + 18/(80*12) + 6/1728
    b4 = (5*(7/480) + 7/96) / r**4
    b6plus = (28*moment6 / r**6) / (1 - (45/28)*0.75/r**2)
    return b4, b6plus


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=Path('results'))
    out = parser.parse_args().out
    out.mkdir(parents=True, exist_ok=True)
    print(f'Python={platform.python_version()} NumPy={np.__version__} '
          f'Matplotlib={matplotlib.__version__} native={accel_two_body.__module__}', flush=True)
    controls()
    near_table(out)
    contours(out, 200, False)
    contours(out, 400, False)
    contours(out, 800, True)
    coarse = scan(out, 90, 200)
    fine = scan(out, 180, 400)
    for metric in METRICS:
        key = (metric, 0.1)
        if coarse[key] is not None and fine[key] is not None:
            print(f'crossing_refinement metric={metric} fine_minus_coarse={fine[key]-coarse[key]:.9g}', flush=True)
    shells(out, fine[('magnitude_cube', 0.1)])
    b4, tail = tail_bound(4.0)
    b = b4 + tail
    print(f'tail r>=4 B4={b4:.12g} B6plus={tail:.12g} '
          f'mascon_normalized_bound_pct={100*b:.12g} cube_normalized_bound_pct={100*b/(1-b):.12g}', flush=True)
    if not 100*b/(1-b) < 0.1:
        raise AssertionError('analytic bound does not certify the tail')


if __name__ == '__main__':
    main()
