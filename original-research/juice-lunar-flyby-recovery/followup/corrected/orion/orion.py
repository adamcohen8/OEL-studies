"""Retrospective short-arc test using selected, timestamped Orion astrometry."""

from pathlib import Path
from datetime import datetime, timezone
import re
import json
import sys
import numpy as np
from scipy.optimize import least_squares, brentq
from scipy.spatial.distance import pdist, squareform

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
sys.path.insert(0, str(ROOT))
from sim.dynamics.orbit.propagator import OrbitPropagator, third_body_moon_plugin, third_body_sun_plugin, j2_plugin
from sim.dynamics.orbit.accelerations import OrbitContext
from sim.dynamics.orbit.environment import EARTH_MU_KM3_S2
from sim.dynamics.orbit.de440_hpop import hpop_de440_positions_km
from sim.dynamics.orbit.frames import eci_to_ecef_rotation_hpop_like

ASEC = 206264.80624709636
LIGHT_SPEED = 299792.458
DAY = 86400.0


def unit(x):
    return x / np.linalg.norm(x, axis=-1, keepdims=True)


def basis(u):
    east = unit(np.cross([0.0, 0.0, 1.0], u))
    return east, np.cross(u, east)


def distance(a, b):
    return ASEC * np.arctan2(np.linalg.norm(np.cross(a, b), axis=-1), np.sum(a * b, axis=-1))


def sexagesimal(s):
    a, b, c = map(float, s.split())
    return np.sign(a) * (abs(a) + b / 60 + c / 3600)


def load():
    records = []
    for m in re.finditer(r'<a name="o(\d+)"></a>(.*)', (HERE / "orion.html").read_text()):
        idx = int(m[1])
        if not 19 <= idx <= 43:
            continue
        line = re.sub("<[^>]*>", "", m[2])
        if line.endswith("U52"):
            mpc_time, ra, dec = line[15:56].split()
            stamp = datetime.strptime(mpc_time, "K%y%m%d:%H%M%S").replace(tzinfo=timezone.utc).timestamp()
            records.append((idx, "U52", stamp, np.deg2rad(float(ra)), np.deg2rad(float(dec))))
            continue
        year, month, day = line[15:32].split()
        stamp = (
            datetime(int(year), int(month), int(float(day)), tzinfo=timezone.utc).timestamp() + (float(day) % 1) * DAY
        )
        ra, dec = np.deg2rad([15 * sexagesimal(line[32:44]), sexagesimal(line[44:56])])
        records.append((idx, line[-3:], stamp, ra, dec))
    assert [r[1] for r in records] == ["204"] * 5 + ["I47"] * 3 + ["U52"] * 12 + ["J69"] * 5
    return records


records = load()
epoch = records[0][2]
JD = 2440587.5 + epoch / DAY
times = np.array([r[2] - epoch for r in records])
ras, decs = np.array([r[3:] for r in records]).T
obs = np.array([np.cos(ras) * np.cos(decs), np.sin(ras) * np.cos(decs), np.sin(decs)]).T
east, north = basis(obs)
eop = np.array(
    [
        [float(x[a:b]) for a, b in [(7, 15), (18, 27), (37, 46), (58, 68)]]
        for x in (HERE / "eop.txt").read_text().splitlines()
    ]
)
fixed = {}
for line in (HERE / "stations.txt").read_text().splitlines():
    lon = np.deg2rad(float(line[4:13]))
    cosphi, sinphi = float(line[13:21]), float(line[21:30])
    fixed[line[:3]] = 6378.14 * np.array([np.cos(lon) * cosphi, np.sin(lon) * cosphi, sinphi])


def site(code, t):
    xp, yp, ut1 = [np.interp(JD + t / DAY - 2400000.5, eop[:, 0], a) for a in eop[:, 1:].T]
    rot = eci_to_ecef_rotation_hpop_like(t, JD, dut1_s=ut1, xp_arcsec=xp, yp_arcsec=yp, dat_s=37)
    return rot.T @ fixed[code]


sites = np.array([site(r[1], t) for r, t in zip(records, times)])
ephem = dict(de440_coeff_path=str(HERE.parent / "de440-geocentric.npz"), de440_tai_utc_s=37.0)
earth_v = np.array(
    [
        (
            hpop_de440_positions_km(JD + (t + 60) / DAY, ephem)["earth"]
            - hpop_de440_positions_km(JD + (t - 60) / DAY, ephem)["earth"]
        )
        / 120
        for t in times
    ]
)
env = dict(ephem, ephemeris_mode="de440", jd_utc_start=JD)
ctx = OrbitContext(mu_km3_s2=EARTH_MU_KM3_S2, mass_kg=1)
zero = np.zeros(3)
engine = OrbitPropagator(
    integrator="rkf78",
    numeric_backend="rust",
    plugins=[j2_plugin, third_body_moon_plugin, third_body_sun_plugin],
    adaptive_atol=1e-11,
    adaptive_rtol=1e-12,
)


def sky(state, n):
    x = state.copy()
    result = []
    last = 0.0
    for i, t in enumerate(times[:n]):
        if t > last:
            x = engine.propagate(x, t - last, last, zero, env, ctx)
        delta = x[:3] - sites[i]
        vel = x[3:] + earth_v[i]
        lag = np.linalg.norm(delta) / LIGHT_SPEED
        for _ in range(2):
            lag = np.linalg.norm(delta - lag * vel) / LIGHT_SPEED
        result.append(unit(delta - lag * vel))
        last = t
    return np.array(result)


def residual(u, n):
    return (
        ASEC
        * np.column_stack([np.sum(u * east[:n], axis=1), np.sum(u * north[:n], axis=1)])
        / np.sum(u * obs[:n], axis=1)[:, None]
    )


# Approximate angular velocity from the first tracklet. Range and range rate are varied.
angular_slope = np.polyfit(times[:5], obs[:5], 1)[0]
site_v = (site("204", 1) - site("204", -1)) / 2


def start_state(r, v):
    # Convert the first received astrometric ray to a reception-epoch state.
    vel = site_v + v * obs[0] + r * angular_slope
    return np.r_[sites[0] + r * obs[0] + (earth_v[0] + vel) * (r / LIGHT_SPEED), vel]


# Local Cartesian corrections rotate into radial and two transverse components.
axes = np.array([obs[0], east[0], north[0]])


def correction(p):
    p = np.asarray(p)
    return np.r_[(p[:3] * [30000, 3, 3]) @ axes, (p[3:] * [0.35, 0.003, 0.003]) @ axes]


# Broad return-coast prior, measured through the end of the UK tracklet.
# Large margins between 15-minute samples cover motion to the nearest node.
coast_times = np.unique(np.r_[np.arange(0, times[-1], 900.0), times[-1]])
coast_moon = np.array([hpop_de440_positions_km(JD + t / DAY, ephem)["moon"] for t in coast_times])


def admissible_coast(state):
    last = 0.0
    for t, moon in zip(coast_times, coast_moon):
        if t > last:
            state = engine.propagate(state, t - last, last, zero, env, ctx)
        if not (
            201500 < np.linalg.norm(state[:3]) < 498500
            and np.linalg.norm(state[:3] - moon) > 102000
            and np.linalg.norm(state[3:]) < 2.99
        ):
            return False
        last = t
    return True


def run(n):
    fits = []
    # Same range/range-rate seeds for both cuts. No held-out measurement enters the fit.
    mids = [(r, v) for r in (185000.0, 285000.0, 385000.0) for v in (-1.5, -0.5, 0.5)]
    for r, v in mids:
        start = start_state(r, v)

        def fun(p):
            return residual(sky(start + correction(p), n), n).ravel()

        # Profile the four transverse coordinates with radial coordinates fixed.
        def at_rv(q):
            return fun([0, *q[:2], 0, *q[2:]])

        opt = least_squares(at_rv, np.zeros(4), jac="3-point", max_nfev=30)
        p = np.r_[0, opt.x[:2], 0, opt.x[2:]]
        fits.append((start + correction(p), "fixed range and range rate"))
        opt = least_squares(
            fun,
            np.clip(p, [-3, -100, -100, -6, -100, -100], [3, 100, 100, 6, 100, 100]),
            jac="3-point",
            max_nfev=60,
            bounds=([-3, -100, -100, -6, -100, -100], [3, 100, 100, 6, 100, 100]),
        )
        fits.append((start + correction(opt.x), "full least squares"))
    accepted = []
    for x, kind in fits:
        fit = distance(sky(x, n), obs[:n])
        if max(fit) < 2.0 and admissible_coast(x):
            accepted.append((x, kind))
    if not accepted:
        raise ValueError("No admissible fits")
    best = min(accepted, key=lambda xf: np.sum(distance(sky(xf[0], n), obs[:n]) ** 2))[0]

    # Expand a pair at fixed offsets in the two weak directions of the fit.
    def centered(p):
        return residual(sky(best + correction(p), n), n).ravel()

    local = least_squares(centered, np.zeros(6), jac="3-point", max_nfev=1)
    _, s, vt = np.linalg.svd(local.jac, full_matrices=False)
    for vec, sing in zip(vt[-2:], s[-2:]):
        for sign in (-1.0, 1.0):

            def excess(scale):
                x = best + correction(sign * scale * vec)
                return max(distance(sky(x, n), obs[:n])) - 1.97

            if excess(0) >= 0:
                continue
            hi = min(0.5 / max(sing, 1e-12), 0.1)
            while hi < 3 and excess(hi) < 0:
                hi = min(3.0, 2 * hi)
            scale = brentq(excess, 0, hi) if excess(hi) > 0 else hi
            for factor in (1.0, 0.5):
                x = best + correction(sign * scale * factor * vec)
                if max(distance(sky(x, n), obs[:n])) < 2 and admissible_coast(x):
                    accepted.append((x, "profile along local singular direction"))
    predictions = np.array([sky(x, 21)[20] for x, _ in accepted])
    dist = squareform(pdist(predictions))
    i, j = np.unravel_index(np.argmax(dist), dist.shape)

    def describe(pair):
        x, kind = pair
        u = sky(x, 25)
        errors = distance(u, obs)
        return dict(
            method=kind,
            state_eci_km_km_s=x.tolist(),
            fit_max_as=max(errors[:n]),
            fit_rms_as=np.sqrt(np.mean(errors[:n] ** 2)),
            first_holdout_los=u[20].tolist(),
            residuals_as=errors.tolist(),
        )

    # Linearized Gaussian baseline: illustrative, independent isotropic 1 arcsecond per coordinate.
    nominal = sky(best, 21)[20]
    tangent = np.array(basis(nominal))

    def holding(p):
        u = sky(best + correction(p), 21)[20]
        return ASEC * (tangent @ u) / (u @ nominal)

    h = 0.001
    derivative = np.column_stack([(holding(h * k) - holding(-h * k)) / (2 * h) for k in np.eye(6)])
    proj = derivative @ vt.T @ np.diag(1 / s)
    sigma_max = np.linalg.svd(proj, compute_uv=False)[0]
    return dict(
        training_rows=n,
        accepted=len(accepted),
        baseline=describe((best, "least squares")),
        pair=[describe(accepted[i]), describe(accepted[j])],
        separation_as=distance(predictions[i], predictions[j]),
        formal_major_axis_1sigma_as=sigma_max,
    )


if __name__ == "__main__":
    result = dict(
        rows=[
            dict(
                source_row=r[0],
                station=r[1],
                utc=datetime.fromtimestamp(r[2], timezone.utc).isoformat(),
                ra_deg=np.rad2deg(r[3]),
                dec_deg=np.rad2deg(r[4]),
            )
            for r in records
        ]
    )
    for key, n in [("single", 5), ("double", 8), ("triple", 20)]:
        result[key] = run(n)
        print(key, result[key]["separation_as"], result[key]["baseline"]["residuals_as"][20], flush=True)
    (HERE / "orion-results.json").write_text(json.dumps(result, indent=2) + "\n")
