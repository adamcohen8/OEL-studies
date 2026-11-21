"""Conservative kinematic enclosure for the two- and three-station Orion arcs.

The linear program is a relaxation of the equations of motion. Its result
requires the continuous state and measurement bounds specified below.
"""

import json
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
import orion as o

R_MAX = 500000.0  # Earth-centred km
V_MAX = 3.0  # Earth-centred km/s
A_MAX = 0.000012  # geocentric km/s^2
SITE_MAX = 6400.0
RAY_MAX = 507000.0  # km, includes barycentric translation during light time
# For angular caps of 2", ray error < 4.92 km. During light time < 1.70 s,
# geocentric object motion < 5.10 km. Earth velocity approximation, timescale
# rounding and numeric constants have > 1.9 km additional component margin.
TUBE = 12.0
# Entire reception-to-emission displacement of the barycentric object.
PREDICT_PAD = 60.0  # km, with Earth barycentric speed <32 km/s


def enclose(n, fit, refine=False, angular_cap_as=2.0):
    t = np.unique(np.r_[o.times[:n], np.arange(3600, o.times[20], 3600) if refine else [], o.times[20]])
    count = len(t)
    dim = 6 * count + n
    lower = np.r_[np.tile([-R_MAX] * 3 + [-V_MAX] * 3, count), np.zeros(n)]
    upper = np.r_[np.tile([R_MAX] * 3 + [V_MAX] * 3, count), np.full(n, RAY_MAX)]
    A = []
    b = []

    def enclosure(row, mid, half):
        A.extend((row, -row))
        b.extend((mid + half, -mid + half))

    for k in range(n):
        for axis in range(3):
            row = np.zeros(dim)
            row[6 * np.searchsorted(t, o.times[k]) + axis] = 1.0
            row[6 * count + k] = -(o.obs[k, axis] + o.earth_v[k, axis] / o.LIGHT_SPEED)
            if refine:
                # For true ray length in [193500,507000] km, approximate
                # (rho/c)*v by mid-lag*v. Error <= 1.57 km for speed <=3.
                row[6 * np.searchsorted(t, o.times[k]) + 3 + axis] = -(193500.0 + 507000.0) / (2 * o.LIGHT_SPEED)
                for sign in (-1.0, 1.0):
                    side = sign * row.copy()
                    # Euclidean chord error per component <= angular radius in radians.
                    side[6 * count + k] -= angular_cap_as / o.ASEC
                    A.append(side)
                    b.append(sign * o.sites[k, axis] + 1.70)
            else:
                enclosure(row, o.sites[k, axis], TUBE)
    for k, h in enumerate(np.diff(t)):
        for axis in range(3):
            x = 6 * k + axis
            row = np.zeros(dim)
            row[x + 6], row[x], row[x + 3] = 1.0, -1.0, -h
            enclosure(row, 0.0, A_MAX * h * h / 2 + 0.001)
            row = np.zeros(dim)
            row[x + 9], row[x + 3] = 1.0, -1.0
            enclosure(row, 0.0, A_MAX * h + 1e-8)
    base = fit["baseline"]["state_eci_km_km_s"]
    center = o.sky(np.array(base), 21)[20]
    east, north = o.basis(center)
    station = o.sites[20]

    def maximize(axis, k):
        c = np.zeros(dim)
        c[6 * k : 6 * k + 3] = axis
        result = linprog(-c, A_ub=matrix, b_ub=rhs, bounds=list(zip(lower, upper)), method="highs")
        if not result.success:
            raise ValueError(result.message)
        # Any nonnegative multipliers give an upper bound after adding the
        # box support of the dual residual. Evaluate with exact float ratios.
        f = Fraction.from_float
        lam = np.maximum(-result.ineqlin.marginals, 0)
        used = np.flatnonzero(lam)
        z = [f(float(val)) for val in c]
        answer = Fraction()
        for j in used:
            weight = f(float(lam[j]))
            answer += weight * f(rhs[j])
            for col in np.flatnonzero(matrix[j]):
                z[col] -= weight * f(matrix[j, col])
        answer += sum(max(val * f(lo), val * f(hi)) for val, lo, hi in zip(z, lower, upper))
        return float(answer) + 0.001

    def support(axis):
        return maximize(axis, count - 1) - axis @ station + PREDICT_PAD + 0.01

    def attraction(lo, hi, mu, floor):
        # Outward enclosure of -mu*x/|x|^3, using a physical distance floor.
        near = np.maximum(0, np.maximum(lo, -hi))
        rmin = max(floor, np.linalg.norm(near) - 0.01)
        rmax = np.linalg.norm(np.maximum(abs(lo), abs(hi))) + 0.01
        terms = np.array([-mu * x / r**3 for x in (lo, hi) for r in (rmin, rmax)])
        return terms.min(axis=0), terms.max(axis=0)

    matrix, rhs = np.array(A), np.array(b)
    for iteration in range(5 if refine else 0):
        maxima = np.array([[maximize(axis, k) for axis in np.eye(3)] for k in range(count)])
        minima = np.array([[-maximize(-axis, k) for axis in np.eye(3)] for k in range(count)])
        for k, h in enumerate(np.diff(t)):
            pad = V_MAX * h / 2 + 1.0
            lo = np.minimum(minima[k], minima[k + 1]) - pad
            hi = np.maximum(maxima[k], maxima[k + 1]) + pad
            alo, ahi = attraction(lo, hi, o.EARTH_MU_KM3_S2, 200000.0)
            midpoint = (t[k] + t[k + 1]) / 2
            moon = o.hpop_de440_positions_km(o.JD + midpoint / o.DAY, o.ephem)["moon"]
            lunar_pad = 1.4 * h / 2 + 1.0
            mlo, mhi = attraction(lo - moon - lunar_pad, hi - moon + lunar_pad, 4902.8001, 100000.0)
            # Uniform component margin for indirect lunar acceleration, solar tide,
            # J2, lunar-GM rounding and numeric enclosure error.
            alo += mlo - 1.2e-7
            ahi += mhi + 1.2e-7
            for axis in range(3):
                mid = (ahi[axis] + alo[axis]) / 2
                half = (ahi[axis] - alo[axis]) / 2
                col = 6 * k + axis
                row = np.zeros(dim)
                row[col + 6], row[col], row[col + 3] = 1.0, -1.0, -h
                enclosure(row, 0.5 * mid * h * h, 0.5 * half * h * h + 0.001)
                row = np.zeros(dim)
                row[col + 9], row[col + 3] = 1.0, -1.0
                enclosure(row, mid * h, half * h + 1e-8)
        matrix, rhs = np.array(A), np.array(b)
    e = max(support(east), support(-east))
    northabs = max(support(north), support(-north))
    forward = -support(-center)
    radius = np.rad2deg(np.arctan2(np.hypot(e, northabs), forward)) if forward > 0 else 180.0
    if forward > 0:
        # Physical radius floor minus site radius and barycentric light-time displacement.
        radius = min(radius, np.rad2deg(np.arcsin(min(1.0, np.hypot(e, northabs) / 193500.0))))
    return dict(
        fixed_center=center.tolist(),
        transverse_east_bound_km=e,
        transverse_north_bound_km=northabs,
        forward_lower_km=forward,
        sufficient_field_radius_deg=float(radius),
    )


if __name__ == "__main__":
    data = json.loads((o.HERE / "orion-results.json").read_text())
    answer = {}
    for key, n in (("double", 8), ("triple", 20)):
        answer[key] = enclose(n, data[key], refine=True)
        answer[key]["acceleration_only_radius_deg"] = enclose(n, data[key])["sufficient_field_radius_deg"]
    answer["triple"]["radius_with_extra_1as_bias_deg"] = enclose(20, data["triple"], refine=True, angular_cap_as=3.0)[
        "sufficient_field_radius_deg"
    ]
    print(answer)
    (o.HERE / "outer.json").write_text(json.dumps(answer, indent=2) + "\n")
