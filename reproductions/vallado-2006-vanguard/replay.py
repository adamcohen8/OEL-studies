import csv
import math
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import sim
from sim.dynamics.orbit.tle import parse_tle_lines
from sim.dynamics.orbit.ogp import (
    ogp_propagate_teme, ogp_propagator_name_for_elements,
)

ROOT = Path(__file__).resolve().parent
EXPECTED = [
    '1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753',
    '2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667',
]
TIMES = tuple(range(0, 4321, 360))
HEADER = 't_min,x_km,y_km,z_km,vx_km_s,vy_km_s,vz_km_s'.split(',')


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def main():
    out = ROOT / 'states.csv'
    out.unlink(missing_ok=True)
    lines = (ROOT / 'input.tle').read_text(encoding='ascii').splitlines()
    check(lines == EXPECTED, 'input.tle differs from the requested TLE')
    elements = parse_tle_lines(*lines, require_checksum=True)
    print(f'Python {sys.version.split()[0]}; sim source={sim.__file__}')
    for dist in ('orbital-engineering-lab', 'oel'):
        try:
            print(f'{dist} distribution version={version(dist)}')
        except PackageNotFoundError:
            pass
    name = ogp_propagator_name_for_elements(elements)
    print(f'Exact TLE and parser checksums: PASS; selected propagator={name!r}; backend=python')
    check(name == 'OGP-SGP4', f'Expected OGP-SGP4, got {name!r}')

    def sample(t, backend='python'):
        state = ogp_propagate_teme(elements, t, backend=backend)
        check(state.error is None, f'{backend} propagation at {t} min returned {state.error!r}')
        r = tuple(map(float, state.position_teme_km))
        v = tuple(map(float, state.velocity_teme_km_s))
        check(len(r) == len(v) == 3 and all(map(math.isfinite, r + v)), f'Invalid state at {t} min')
        return r + v

    states = {t: sample(t) for t in TIMES}
    rust_states = {t: sample(t, backend='rust') for t in TIMES}
    dr = {t: math.dist(states[t][:3], rust_states[t][:3]) for t in TIMES}
    dv = {t: math.dist(states[t][3:], rust_states[t][3:]) for t in TIMES}
    tr, tv = max(dr, key=dr.get), max(dv, key=dv.get)
    print(f'Max Python/Rust position difference: {dr[tr]:.9g} km at {tr} min')
    print(f'Max Python/Rust velocity difference: {dv[tv]:.9g} km/s at {tv} min')
    cross_ok = dr[tr] <= 1e-5 and dv[tv] <= 1e-8
    print(f'Cross-backend agreement at 1e-5 km and 1e-8 km/s: {"PASS" if cross_ok else "FAIL"}')
    residuals, step_changes = [], []
    for t in TIMES:
        def derivative(h_s):
            a, b = sample(t - h_s / 60), sample(t + h_s / 60)
            return tuple((b[i] - a[i]) / (2 * h_s) for i in range(3))
        coarse, fine = derivative(2.0), derivative(1.0)
        residuals.append(math.dist(fine, states[t][3:]))
        step_changes.append(math.dist(coarse, fine))
    worst = max(range(len(TIMES)), key=residuals.__getitem__)
    print(f'Max centered-difference velocity residual: {residuals[worst]:.9g} km/s at {TIMES[worst]} min')
    print(f'Max change in derivative on halving step: {max(step_changes):.9g} km/s')
    check(max(step_changes) < 1e-5, 'Finite difference is too sensitive to step size')
    check(max(residuals) < 2e-3, 'Position and velocity fail the kinematic consistency check')
    check(all(sample(t) == states[t] for t in reversed(TIMES)), 'Reverse-order replay differs')
    print('Reverse-order replay: bit-for-bit numeric equality')
    rows = [(t, *states[t]) for t in TIMES]
    temp = out.with_suffix('.tmp')
    with temp.open('w', newline='', encoding='ascii') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(HEADER)
        w.writerows(rows)
    with temp.open(newline='', encoding='ascii') as f:
        reader = csv.reader(f)
        check(next(reader) == HEADER, 'Header mismatch')
        check([(int(t), *map(float, vals)) for t, *vals in reader] == rows, 'CSV round-trip mismatch')
    temp.replace(out)
    radii = [math.hypot(*s[:3]) for s in states.values()]
    speeds = [math.hypot(*s[3:]) for s in states.values()]
    print(f'Wrote {len(rows)} checked rows to {out.name}; range {TIMES[0]}..{TIMES[-1]} min')
    print(f'Radius range {min(radii):.6f}..{max(radii):.6f} km; speed range {min(speeds):.6f}..{max(speeds):.6f} km/s')


if __name__ == '__main__':
    main()
