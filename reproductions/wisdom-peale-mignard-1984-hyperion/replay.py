#!/usr/bin/env python3
"""Blinded Saturn planar spin–orbit study; Python 3, NumPy/SciPy/Matplotlib.

Recognition: this is the planar gravity-gradient spin–orbit problem; the
parameters remind me of Hyperion. No reference result was used, and no
particular island or chaotic component was assumed in designing the study.

From the study workspace, with OEL_SOURCE set to the pinned source root:
    python replay.py pilot
    python replay.py production --help
    python replay.py refine --reference outputs/production

See the workspace README for the exact commands, initial states and work limits.
All output goes in a new outputs/<name> directory. A run is complete only
when summary.json exists. NPZ files contain only numeric section data; PNGs
are headless. All times and states plotted or classified come from OEL.

Coordinates are long-axis inertial longitude modulo pi, in [-pi/2,pi/2),
and positive body-Z spin divided by mean motion. The longitude is also the
relative angle at the initial periapse. At apoapse it agrees modulo pi up
to orbital error. Nominal T sampling is appropriate for the analytic fixed
Kepler orbit; refinement also finds radial-velocity zeroes in the *native*
orbit without feeding those event states back into the main trajectory.

The map has no dissipation. 'Supported' means periodic points and/or bounded
resonant or circulating structures, not attractors. A near return or a
diffuse section alone is not a stability or chaos diagnosis. Use the native
periodic-map Jacobian, same-seed refinement, rotational rates over separate
windows, and finite-time separation together. A finite survey cannot exclude
small islands, sticky trajectories, or resonances outside its search domain.
"""

import argparse
import dataclasses
import importlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
from scipy.optimize import brentq


PIN = 'c1724805b97937d87f0ee26ea8b77a74ebc7678d'
VERSION = '1.2.0+preview.1'
MU = 37931206.23436167
AORB, ECC, SIGMA = 150000.0, 0.1, 0.89
N = math.sqrt(MU / AORB**3)
T = 2 * math.pi / N
RP = AORB * (1 - ECC)
VP = math.sqrt(MU * (1 + ECC) / RP)
INERTIA = np.diag([0.55, 0.55 + SIGMA**2 / 3, 1.0]).tolist()
CHUNKS = 16  # native calls per orbit; endpoints also track quaternion winding
DEFAULTS = {'pilot': (56, 64), 'production': (176, 96), 'periodic': (0, 128)}


def wrap(x, period=math.pi):
    return (x + period / 2) % period - period / 2


def update(state, **changes):
    if dataclasses.is_dataclass(state):
        return dataclasses.replace(state, **changes)
    if hasattr(state, 'model_copy'):
        return state.model_copy(update=changes)
    if hasattr(state, '_replace'):
        return state._replace(**changes)
    import copy
    result = copy.copy(state)
    for key, val in changes.items():
        setattr(result, key, val)
    return result


def typed_like(old, values):
    if isinstance(old, tuple):
        return tuple(float(x) for x in values)
    if isinstance(old, list):
        return [float(x) for x in values]
    return np.array(values, dtype=float, copy=True)


class Engine:
    def __init__(self, api, config, sign=1):
        validate, prepare, create = api
        cfg = validate(config)
        rt = create('probe', cfg.objects['probe'], cfg, prepare(cfg))
        self.initial = rt.truth
        self.dyn = rt.dynamics
        self.sign = sign
        self.calls = 0

    def step_to(self, state, t):
        dt = float(t - state.t_s)
        if dt == 0:
            return state
        if dt < 0:
            raise ValueError('Native step must advance time')
        self.calls += 1
        o, a = self.dyn.attitude.propagate(self.dyn.orbit_propagator, state, dt)
        o, a = np.asarray(o).reshape(6), np.asarray(a).reshape(7)
        if not (np.isfinite(o).all() and np.isfinite(a).all()):
            raise FloatingPointError('Nonfinite native state')
        return update(state, t_s=float(t),
                      position_eci_km=typed_like(state.position_eci_km, o[:3]),
                      velocity_eci_km_s=typed_like(state.velocity_eci_km_s, o[3:]),
                      attitude_quat_bn=typed_like(state.attitude_quat_bn, a[:4]),
                      angular_rate_body_rad_s=typed_like(state.angular_rate_body_rad_s, a[4:]))

    def inject(self, state, theta, w):
        p = self.sign * theta / 2
        return update(state,
                      attitude_quat_bn=typed_like(state.attitude_quat_bn,
                                                  [math.cos(p), 0, 0, math.sin(p)]),
                      angular_rate_body_rad_s=typed_like(state.angular_rate_body_rad_s,
                                                        [0, 0, w * N]))

    def coordinates(self, state):
        q0, _, _, qz = state.attitude_quat_bn
        return self.sign * 2 * math.atan2(qz, q0), state.angular_rate_body_rad_s[2] / N

    def grid_step(self, state, j):
        return self.step_to(state, self.initial.t_s + j * T / CHUNKS)


def checks(api, config):
    import copy
    cfg_free = copy.deepcopy(config)
    cfg_free['simulator']['dynamics']['attitude']['disturbance_torques']['gravity_gradient'] = False
    free = Engine(api, cfg_free)
    dt = T * 1e-5
    s = free.step_to(free.initial, free.initial.t_s + dt)
    measured = free.coordinates(s)[0] / (N * dt)
    if not math.isclose(abs(measured), 1, abs_tol=1e-5):
        raise RuntimeError('Spin/quaternion convention check failed')
    sign = 1 if measured > 0 else -1
    eng = Engine(api, config, sign)
    s0 = eng.inject(eng.initial, math.pi / 8, 0)
    s1 = eng.step_to(s0, s0.t_s + dt)
    accel = eng.coordinates(s1)[1] / (N * dt)
    # Independent instantaneous Euler-equation check, not a trajectory solver.
    expected = -0.5 * SIGMA**2 * (AORB / RP)**3 * math.sin(math.pi / 4)
    if not math.isclose(accel, expected, rel_tol=5e-4):
        raise RuntimeError(f'Gravity-gradient sign/scale check failed: {accel} vs {expected}')
    return eng, {'unsigned_convention_raw_angle_rate_over_n': measured, 'angle_sign': sign,
                 'native_small_step_accel_over_n2': accel,
                 'instantaneous_euler_accel_over_n2': expected}


def radial_product(s):
    return float(np.dot(s.position_eci_km, s.velocity_eci_km_s))


def event_sample(eng, before, nominal):
    # Solve on a short native branch on either side of the nominal sample.
    # Endpoint branches reuse the same preceding grid checkpoint.
    tn = nominal.t_s
    delta = T / 2048
    lo, hi = tn - delta, tn + delta
    left, right = eng.step_to(before, lo), eng.step_to(before, hi)
    if not radial_product(left) < 0 < radial_product(right):
        return None
    te = brentq(lambda t: radial_product(eng.step_to(before, t)), lo, hi,
                xtol=1e-7, rtol=1e-14)
    event = eng.step_to(before, te)
    th, w = eng.coordinates(event)
    th0, w0 = eng.coordinates(nominal)
    return [te - tn, wrap(th), w, wrap(th - th0), w - w0,
            radial_product(event) / (AORB * VP)]


DIAG_COLUMNS = ['peri_position_error_over_a', 'peri_velocity_error_over_vp',
                'relative_specific_energy_error', 'relative_angular_momentum_error',
                'out_of_plane_q_norm', 'out_of_plane_spin_over_n', 'quaternion_norm_error']
EVENT_COLUMNS = ['orbit_index', 'time_shift_s', 'theta', 'w', 'theta_shift', 'w_shift',
                 'radial_product_over_a_vp']


def diagnostics(s):
    r, v = np.asarray(s.position_eci_km), np.asarray(s.velocity_eci_km_s)
    q, spin = np.asarray(s.attitude_quat_bn), np.asarray(s.angular_rate_body_rad_s)
    energy = np.dot(v, v) / 2 - MU / np.linalg.norm(r)
    return [np.linalg.norm(r - [RP, 0, 0]) / AORB,
            np.linalg.norm(v - [0, VP, 0]) / VP,
            abs(energy / (-MU / (2 * AORB)) - 1),
            np.linalg.norm(np.cross(r, v) / (RP * VP) - [0, 0, 1]),
            np.linalg.norm(q[1:3]), np.linalg.norm(spin[:2]) / N,
            abs(np.linalg.norm(q) - 1)]


def increment(eng, old, new):
    a, wa = eng.coordinates(old)
    b, wb = eng.coordinates(new)
    dt = (new.t_s - old.t_s) * N
    # Bound the angular displacement using the maximum possible periapse torque.
    # If needed, inspect a native intermediate state to resolve the winding.
    max_accel = 0.5 * SIGMA**2 / (1 - ECC)**3
    if abs(wa) * dt + 0.5 * max_accel * dt**2 >= math.pi / 2:
        mid = eng.step_to(old, 0.5 * (old.t_s + new.t_s))
        return increment(eng, old, mid) + increment(eng, mid, new)
    return wrap(b - a, 2 * math.pi)


def track(eng, seed, count, event_every, lyap):
    base = eng.inject(eng.initial, *seed)
    eps = 1e-7
    shadow = eng.inject(base, seed[0] + eps, seed[1]) if lyap else None
    peri, apo, diag, events, growth = [[wrap(seed[0]), seed[1], 0.0]], [], [], [], []
    winding = 0.0
    missing = 0
    for j in range(1, count * CHUNKS + 1):
        before = base
        base = eng.grid_step(base, j)
        winding += increment(eng, before, base)
        if shadow is not None:
            shadow = eng.grid_step(shadow, j)
        th, w = eng.coordinates(base)
        if j % CHUNKS == CHUNKS // 2:
            apo.append([wrap(th), w])
        if j % CHUNKS:
            continue
        k = j // CHUNKS
        peri.append([wrap(th), w, winding / (2 * math.pi)])
        diag.append(diagnostics(base))
        if event_every and k % event_every == 0:
            row = event_sample(eng, before, base)
            if row is None:
                missing += 1
            else:
                events.append([k, *row])
        if shadow is not None:
            ts, ws = eng.coordinates(shadow)
            dtheta, dw = wrap(ts - th), ws - w
            norm = math.hypot(dtheta, dw)
            if not norm or not math.isfinite(norm):
                raise FloatingPointError('Tangent probe separation unresolved')
            growth.append(math.log(norm / eps))
            shadow = eng.inject(base, th + eps * dtheta / norm, w + eps * dw / norm)
    return dict(peri=np.asarray(peri), apo=np.asarray(apo), diag=np.asarray(diag),
                events=np.asarray(events).reshape(-1, len(EVENT_COLUMNS)),
                log_growth=np.asarray(growth)), missing


def section_delta(a, b):
    return np.hypot(wrap(a[..., 0] - b[..., 0]), a[..., 1] - b[..., 1])


def describe(data, missing):
    p = data['peri']
    size = len(p) - 1
    splits = [0, size // 2, size]
    rot = [(p[b, 2] - p[a, 2]) / (b - a) for a, b in zip(splits, splits[1:])]
    tail = p[size // 2:]
    recur = [float(np.median(section_delta(tail[k:], tail[:-k])))
             for k in range(1, min(len(tail), 9))]
    ev = data['events']
    log = data['log_growth']
    return {'turns_per_orbit_in_halves': rot, 'tail_spin_range': [tail[:, 1].min(), tail[:, 1].max()],
            'tail_median_return_distances_periods_1_up': recur,
            'max_diagnostics': dict(zip(DIAG_COLUMNS, data['diag'].max(axis=0))),
            'event_samples': len(ev), 'unbracketed_events': missing,
            'max_abs_event_minus_nominal': (dict(zip(EVENT_COLUMNS[1:], np.max(np.abs(ev[:, 1:]), axis=0)))
                                           if len(ev) else None),
            'ftle_per_orbit_in_halves': ([float(x.mean()) for x in np.array_split(log, 2)]
                                         if len(log) else None)}


def default_seeds(mode):
    spins = [-0.5, 0, 0.5, 1, 1.5, 2, 2.65]
    angles = [0, math.pi / 4]
    if mode == 'production':
        spins = [-1, -0.5, 0, 0.45, 0.85, 1, 1.3, 1.5, 1.85, 2.1, 2.7, 3.2]
        angles = [0, math.pi / 5, -math.pi / 3]
    return [[a, w] for a in angles for w in spins]


def choose_refinement(rows):
    def near(r):
        return min(r['tail_median_return_distances_periods_1_up'])

    def wander(r):
        return abs(np.diff(r['turns_per_orbit_in_halves'])[0])

    ordered = sorted(rows, key=lambda r: r['seed'][1])
    middle = [r for r in rows if 0 <= r['seed'][1] <= 2.2] or rows
    picks = [min(middle, key=near), max(middle, key=wander), max(middle, key=near),
             ordered[0], ordered[-1]]
    return [list(x) for x in dict.fromkeys(tuple(r['seed']) for r in picks)]


def compare(fine, reference):
    size = min(len(fine), len(reference))
    d = section_delta(fine[:size], reference[:size])
    mid = (size - 1) // 2
    return {'common_orbits': size - 1, 'first_8_orbits_max_distance': d[:9].max(),
            'common_tail_median_distance': np.median(d[mid:]),
            'first_orbit_distance_exceeds': {str(level): (int(i[0]) if len(i) else None)
                for level in [1e-3, 1e-2, 0.1] for i in [np.flatnonzero(d > level)]},
            'common_tail_turns_per_orbit_fine_reference':
                [(p[size - 1, 2] - p[mid, 2]) / (size - 1 - mid) for p in [fine, reference]],
            'common_tail_spin_quantiles_fine_reference':
                [np.quantile(p[mid:size, 1], [0.1, 0.5, 0.9]).tolist() for p in [fine, reference]]}


def flow(eng, theta, w, orbits):
    state = eng.inject(eng.initial, theta, w)
    winding = 0.0
    for j in range(1, round(orbits * CHUNKS) + 1):
        nx = eng.grid_step(state, j)
        winding += increment(eng, state, nx)
        state = nx
    th, wn = eng.coordinates(state)
    return np.array([wrap(th), wn]), winding / math.pi


def periodic_search(eng, args):
    # At an apsis the Kepler forcing is reflection symmetric. Starting with
    # a principal axis aligned, seek a second alignment (with either axis)
    # at apoapse. Test, rather than assume, closure of the full native map.
    grid = np.linspace(args.search_min, args.search_max, args.search_grid)
    rows = []
    for a in (0.0, math.pi / 2):
        def half(w):
            return a / math.pi + flow(eng, a, float(w), 0.5)[1]

        if CHUNKS % 2:
            raise RuntimeError('The half-period symmetry search requires even CHUNKS')
        values = [half(w) for w in grid]
        for j, (u, v) in enumerate(zip(values, values[1:])):
            for m in range(math.ceil(2 * min(u, v)), math.floor(2 * max(u, v)) + 1):
                w = float(brentq(lambda x: half(x) - m / 2, grid[j], grid[j + 1], xtol=1e-10))
                if any(abs(r['seed'][0] - a) < 1e-10 and abs(r['seed'][1] - w) < 1e-7
                       for r in rows):
                    continue
                p, winding = flow(eng, a, w, 1)
                error = float(section_delta(p, np.array([a, w])))
                row = {'seed': [a, w], 'apo_axis_index_unwrapped': m,
                       'full_map_closure_distance': error, 'turns_per_orbit': winding / 2,
                       'integer_pi_winding_residual': abs(winding - round(winding)),
                       'classification': 'closure_not_verified'}
                rows.append(row)
                if error > 1e-5 or row['integer_pi_winding_residual'] > 1e-5:
                    continue
                results = []
                for h in (1e-5, 5e-6):
                    jac = np.empty((2, 2))
                    for i, (da, dw) in enumerate(((h, 0), (0, h))):
                        plus = flow(eng, a + da, w + dw, 1)[0]
                        minus = flow(eng, a - da, w - dw, 1)[0]
                        jac[:, i] = [wrap(plus[0] - minus[0]) / (2 * h),
                                     (plus[1] - minus[1]) / (2 * h)]
                    results.append([float(np.trace(jac)), float(np.linalg.det(jac))])
                row['trace_det_at_two_fd_scales'] = results
                tr, det = results[-1]
                margin = max(1e-3, 4 * abs(results[0][0] - tr), 4 * abs(det - 1))
                if abs(det - 1) > 1e-3:
                    label = 'unresolved_area_check'
                elif abs(tr) < 2 - margin:
                    label = 'elliptic_linearization'
                elif abs(tr) > 2 + margin:
                    label = 'hyperbolic_linearization'
                else:
                    label = 'near_marginal_or_unresolved'
                row['classification'] = label
    return rows


def plot_sections(path, arrays, seeds, roots):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True, layout='constrained')
    cmap = plt.get_cmap('turbo')
    for i, d in enumerate(arrays):
        color = cmap(i / max(1, len(seeds) - 1))
        for ax, key in zip(axes, ('peri', 'apo')):
            p = d[key]
            ax.scatter(p[:, 0], p[:, 1], s=4, color=color, alpha=0.64, linewidths=0)
        axes[0].scatter(*seeds[i], marker='+', s=23, color=color)
    for root in roots:
        label = root['classification']
        if label in ('elliptic_linearization', 'hyperbolic_linearization'):
            x, y = root['seed']
            axes[0].scatter(wrap(x), y, marker='o' if label.startswith('elliptic') else 'x',
                            facecolors='none' if label.startswith('elliptic') else 'black',
                            edgecolors='black' if label.startswith('elliptic') else None, s=55)
    for ax, title in zip(axes, ('Nominal periapse (+ seeds; o/x verified linear types)',
                               'Nominal apoapse (different forcing phase)')):
        ax.set(xlabel='Inertial long-axis angle modulo π [rad]', title=title,
               xlim=(-math.pi / 2, math.pi / 2))
        ax.grid(alpha=0.2)
    axes[0].set_ylabel('Spin / mean motion')
    fig.savefig(path, dpi=160)
    plt.close(fig)


def replay_api(lockfile):
    import importlib
    import os
    from pathlib import Path
    import subprocess
    import sys

    lock = json.loads(lockfile.read_text(encoding='utf-8'))
    pin = 'c1724805b97937d87f0ee26ea8b77a74ebc7678d'
    if lock['development_commit'] != pin or lock['version'] != '1.2.0+preview.1':
        raise RuntimeError('The workspace lock does not match the study pin')
    if not os.environ.get('OEL_SOURCE'):
        raise RuntimeError('Set OEL_SOURCE to the pinned OEL checkout')
    source = Path(os.environ['OEL_SOURCE']).resolve(strict=True)

    def git(*args):
        return subprocess.check_output(['git', '-C', str(source), *args], text=True).strip()

    if Path(git('rev-parse', '--show-toplevel')).resolve() != source:
        raise RuntimeError('OEL_SOURCE must point to the checkout root')
    head = git('rev-parse', 'HEAD')
    public = lock.get('public_commit')
    if len(head) != 40 or head not in (pin, public):
        raise RuntimeError(f'Unexpected OEL git HEAD: {head}')
    if git('status', '--porcelain', '--untracked-files=no'):
        raise RuntimeError('Pinned OEL checkout has tracked modifications')
    sys.path.insert(0, str(source))
    from importlib.metadata import version
    if version('oel-rust-orbit') != lock['native_orbit_version']:
        raise RuntimeError('Installed native orbit version differs from the study lock')
    methods = []
    for mod, name in [('sim.config.scenario_yaml', 'scenario_config_from_dict'),
                      ('sim.dynamics.planetary.preparation', 'prepare_environment'),
                      ('sim.dynamics.planetary.runtime', 'create_runtime')]:
        module = importlib.import_module(mod)
        if not Path(module.__file__).resolve().is_relative_to(source):
            raise RuntimeError(f'{mod} was not loaded from OEL_SOURCE')
        methods.append(getattr(module, name))

    return tuple(methods), {'git_head': head, 'source_kind': 'development' if head == pin else 'public_export',
                                'lock': lock, 'python': sys.version.split()[0], 'numpy': np.__version__}


def replay_config(count, attitude_step, orbit_integrator):
    inertia = [[0.55, 0., 0.], [0., 0.55 + SIGMA**2 / 3, 0.], [0., 0., 1.]]
    return {'objects': {'probe': {'kind': 'satellite', 'runtime_profile': 'trajectory_only',
                                 'specs': {'mass_kg': 100., 'mass_properties': {
                                     'inertia_reference_point': 'center_of_mass', 'inertia_kg_m2': inertia}},
                                 'initial_state': {'position_inertial_km': [RP, 0., 0.],
                                                   'velocity_inertial_km_s': [0., VP, 0.],
                                                   'attitude_quat_bn': [1., 0., 0., 0.],
                                                   'angular_rate_body_rad_s': [0., 0., N]}}},
            'simulator': {'central_body': 'saturn', 'initial_jd_utc': 2451545.0,
                          'duration_s': (count + 1) * T, 'dt_s': T / CHUNKS,
                          'planetary': {'gravity': {'model': 'point_mass'}},
                          'dynamics': {'orbit': {'integrator': orbit_integrator, 'numeric_backend': 'rust'},
                                       'attitude': {'enabled': True, 'attitude_substep_s': attitude_step,
                                                    'disturbance_torques': {'gravity_gradient': True}}}}}


def replay_cli(argv=None):
    import argparse
    from pathlib import Path

    examples = ('python replay.py pilot\n'
                'python replay.py production --periodic\n'
                'python replay.py refine --reference outputs/production\n'
                'python replay.py refine --reference outputs/production --out outputs/check-dopri5 '
                '--orbit-integrator dopri5 --orbits 96\n'
                'python replay.py production --seed 0 1 --seed 0.4 1.5 --out outputs/targeted')
    parser = argparse.ArgumentParser(description=__doc__, epilog=examples,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    modes = ('pilot', 'production', 'refine', 'refinement')
    parser.add_argument('command', nargs='?', choices=modes)
    parser.add_argument('--mode', choices=modes)
    parser.add_argument('--out', '--output', '--output-dir', type=Path,
                        help='Fresh run directory under outputs/; defaults to outputs/<mode>.')
    parser.add_argument('--reference', '--from-run', '--baseline', type=Path,
                        help='Previous run directory (or its summary.json).')
    parser.add_argument('--orbits', type=int)
    parser.add_argument('--seed', nargs=2, type=float, action='append', metavar=('THETA_RAD', 'SPIN_OVER_N'))
    parser.add_argument('--steps-per-orbit', '--attitude-steps-per-orbit', type=int)
    parser.add_argument('--attitude-substep-s', type=float)
    parser.add_argument('--orbit-integrator', choices=('rkf78', 'rk4', 'dopri5'))
    parser.add_argument('--event-every', type=int, help='Periapse root sampling stride; 0 disables.')
    parser.add_argument('--lyap', action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument('--periodic', action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument('--search-min', type=float, default=-1.15)
    parser.add_argument('--search-max', type=float, default=2.9)
    parser.add_argument('--search-grid', type=int, default=33)
    args = parser.parse_args(argv)
    if args.mode and args.command and args.mode != args.command:
        parser.error('conflicting modes')
    args.mode = args.mode or args.command or 'pilot'
    if args.mode == 'refinement':
        args.mode = 'refine'
    if args.mode == 'refine' and args.reference is None:
        args.reference = Path('outputs/production')
    if args.out is None:
        args.out = Path('outputs') / args.mode
    if args.attitude_substep_s is not None and args.steps_per_orbit is not None:
        parser.error('choose only one attitude step control')
    if not (math.isfinite(args.search_min) and math.isfinite(args.search_max)
            and args.search_min < args.search_max and 3 <= args.search_grid <= 1001):
        parser.error('invalid periodic search range or grid')
    return args


def replay_main(argv=None):
    import hashlib
    import json
    from pathlib import Path
    import sys
    import time

    args = replay_cli(argv)
    source_file = Path(__file__).resolve()
    output_root = Path('outputs').resolve()
    destination = args.out.resolve()
    if destination == output_root or not destination.is_relative_to(output_root):
        raise ValueError('--out must name a directory under outputs/')
    if destination.exists():
        raise FileExistsError(f'{args.out} exists; choose a fresh --out')
    ref = None
    ref_dir = None
    ref_digest = None
    if args.reference is not None:
        ref_file = args.reference if args.reference.is_file() else args.reference / 'summary.json'
        raw = ref_file.read_bytes()
        ref = json.loads(raw)
        ref_dir = ref_file.parent
        ref_digest = hashlib.sha256(raw).hexdigest()
        if ref.get('study_schema') != 1:
            raise ValueError('Incompatible reference summary schema')
    opts = ref['settings'] if ref else {}
    seeds = args.seed if args.seed is not None else (
        choose_refinement(ref['tracks']) if args.mode == 'refine' else default_seeds(args.mode))
    seeds = [[float(wrap(a)), float(w)] for a, w in seeds]
    if not seeds or not all(math.isfinite(a) and math.isfinite(w) for a, w in seeds):
        raise ValueError('At least one finite seed is required')
    count = args.orbits if args.orbits is not None else (
        24 if args.mode == 'pilot' else (max(96, opts['orbits']) if args.mode == 'refine' else 96))
    stride = args.event_every if args.event_every is not None else (4 if args.mode == 'pilot' else 8)
    if not 4 <= count <= 5000 or stride < 0:
        raise ValueError('Require 4 <= --orbits <= 5000 and --event-every >= 0')
    step = args.attitude_substep_s
    if step is None:
        if args.steps_per_orbit is not None:
            if args.steps_per_orbit <= 0:
                raise ValueError('--steps-per-orbit must be positive')
            step = T / args.steps_per_orbit
        elif args.mode == 'refine':
            step = min(opts['attitude_substep_s'], T / CHUNKS) / 2
        else:
            step = T / max(128, 2 * CHUNKS)
    if not math.isfinite(step) or step <= 0:
        raise ValueError('Attitude step must be finite and positive')
    orbit_integrator = args.orbit_integrator or opts.get('orbit_integrator', 'rkf78')
    lyap = args.lyap if args.lyap is not None else args.mode != 'pilot'
    search = args.periodic if args.periodic is not None else (
        args.mode == 'production' or bool(opts.get('periodic')))

    # The API loader checks the checkout before importing any OEL module.
    api, provenance = replay_api(source_file.with_name('oel-revision.lock.json'))
    if ref and ref['provenance']['git_head'] != provenance['git_head']:
        raise RuntimeError('Reference and current checkout differ; do not conflate revision and step effects')
    config = replay_config(count, step, orbit_integrator)
    eng, calibration = checks(api, config)
    start = time.perf_counter()
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = destination.with_name(destination.name + '.in-progress')
    stage.mkdir()

    def dump(filename, obj):
        def convert(x):
            if isinstance(x, np.generic):
                return x.item()
            if isinstance(x, np.ndarray):
                return x.tolist()
            raise TypeError(type(x).__name__)

        (stage / filename).write_text(json.dumps(obj, default=convert, allow_nan=False,
                                                 indent=2, sort_keys=True) + '\n', encoding='utf-8')

    provenance['replay_sha256'] = hashlib.sha256(source_file.read_bytes()).hexdigest()
    dump('scenario.json', config)
    dump('provenance.json', provenance)
    rows, arrays = [], []
    for i, seed in enumerate(seeds):
        data, missing = track(eng, seed, count, stride, lyap)
        name = f'track_{i:03d}.npz'
        np.savez_compressed(stage / name, **data)
        row = {'seed': seed, 'data': name, **describe(data, missing)}
        if ref:
            for old in ref['tracks']:
                if float(section_delta(np.array(seed), np.array(old['seed']))) <= 1e-12:
                    with np.load(ref_dir / old['data'], allow_pickle=False) as native:
                        row['reference_comparison'] = compare(data['peri'], native['peri'])
                    break
        arrays.append(data)
        rows.append(row)
        print(f'{args.mode}: {i + 1}/{len(seeds)} seeds complete', file=sys.stderr, flush=True)

    periodic = periodic_search(eng, args) if search else []
    if search and ref and ref['periodic_orbits']:
        for root in periodic:
            candidates = [r for r in ref['periodic_orbits']
                          if r['apo_axis_index_unwrapped'] == root['apo_axis_index_unwrapped']
                          and abs(wrap(r['seed'][0] - root['seed'][0])) < 1e-10]
            if candidates:
                nearest = min(candidates, key=lambda r: abs(r['seed'][1] - root['seed'][1]))
                shift = root['seed'][1] - nearest['seed'][1]
                if abs(shift) < 0.025:
                    root['reference_match'] = {'spin_shift': shift, 'classification': nearest['classification'],
                                               'closure_distance': nearest['full_map_closure_distance']}
    plot_sections(stage / 'sections.png', arrays, seeds, periodic)
    summary = {
        'study_schema': 1, 'mode': args.mode, 'provenance': provenance,
        'prior_recognition': ('Planar gravity-gradient spin–orbit problem; parameters reminiscent of Hyperion. '
                              'No prior island or chaos identification, expected outcome, or consulted result.'),
        'physics': {'mu_km3_s2': MU, 'a_km': AORB, 'e': ECC, 'sigma': SIGMA,
                    'mean_motion_rad_s': N, 'period_s': T, 'apsidal_direction': '+J2000 X'},
        'settings': {'orbits': count, 'grid_chunks_per_orbit': CHUNKS, 'orbit_integrator': orbit_integrator,
                     'attitude_substep_s': step, 'event_every': stride, 'lyap': lyap, 'periodic': search,
                     'search_min_max_grid': [args.search_min, args.search_max, args.search_grid] if search else None},
        'reference_summary_sha256': ref_digest, 'reference_settings': opts if ref else None,
        'calibration': calibration, 'tracks': rows, 'periodic_orbits': periodic,
        'data_columns': {'peri': ['theta_rad_mod_pi', 'w_over_n', 'unwrapped_turns_from_seed'],
                         'apo': ['theta_rad_mod_pi', 'w_over_n'], 'diag': DIAG_COLUMNS,
                         'events': EVENT_COLUMNS, 'log_growth': ['natural_log_separation_gain_per_orbit']},
        'interpretation_guide': [
            'Every track, shadow, event branch and periodic-map evaluation uses native OEL orbit plus passive attitude.',
            'The Kepler state/period, orbital invariants, instantaneous torque and map area preservation are independent checks, not alternate trajectory integrations.',
            'The section coordinates are the inertial long-axis angle reduced modulo pi and the signed body-Z spin divided by mean motion. Periapse and apoapse sections are distinct drives of the same dynamics and should not be pooled.',
            'Compare the nominal-period section against radial events recorded on selected revolutions before interpreting fine structure. Changing event cadence changes the audit sample; it does not change the nominal-period section.',
            'A closed one-orbit map point with an elliptic linearization is a candidate island center; nearby bounded libration supplies additional evidence. A hyperbolic point is a periodic state, not a stable spin lock. Do not infer nonlinear stability from the linearization alone.',
            'The periodic search is limited to its stated symmetry lines, spin interval and mesh. It can miss roots, including asymmetric and longer-period states. Seed occupancy is not a phase-space measure; transient behavior or a missed island is possible.',
            'Use the unwrapped angle to estimate winding, rather than differencing wrapped section angles. Near-rational finite-time winding does not by itself establish a resonance.',
            'At matched seeds, compare early pointwise differences and then finite-time winding, section morphology and separation statistics. Late pointwise decorrelation alone does not show numerical failure or physical chaos.',
            'A sustained positive finite-time growth rate that survives step refinement and longer observation supports a chaotic interpretation; it is not an asymptotic exponent measurement. Conservative dynamics does not imply capture or attraction.',
        ],
        'wall_seconds_after_calibration': time.perf_counter() - start,
    }
    dump('summary.json', summary)
    stage.rename(destination)
    print(destination.relative_to(output_root.parent) / 'summary.json', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(replay_main())
