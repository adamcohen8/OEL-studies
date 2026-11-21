"""Check sensitivity of the archived least-squares boresight to the light-time fit.

Use the archived optimizer's equal-weight tangent residuals, iterated spacecraft
emission states, and a smooth local Earth displacement for each observation.
This is a local refit starting from the reported solution.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

STUDY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STUDY / "field-certificate"))
from common import Arc, o

arc = Arc("juice")
mission = STUDY / "followup/corrected/missions"
reference = np.array(json.loads((mission / "juice-forecast.json").read_text())["state"])
published = np.array(json.loads((mission / "juice-score.json").read_text())["replay"][0]["los"])


def earth(t):
    return o.hpop_de440_positions_km(arc.JD + t / 86400, arc.ephem)["earth"]


# Quadratic interpolation centered at reception avoids noise when subtracting
# Earth positions evaluated at trial-dependent, single-part Julian dates.
velocities, accelerations = [], []
for t in arc.times:
    before, at, after = earth(t - 64), earth(t), earth(t + 64)
    velocities.append((after - before) / 128)
    accelerations.append(((after - at) + (before - at)) / 64**2)


def directions(p, include_forecast=False):
    previous = -15.0
    state = reference + arc.delta(p)
    state[:3] += previous * state[3:]  # six free coordinates now specify the earlier epoch
    directions = []
    for i, t in enumerate(arc.times if include_forecast else arc.times[:-1]):
        pre = t - 15.0
        if pre > previous:
            state = arc.propagate(state, pre - previous, previous)
        previous = pre
        lag = 6.0
        for _ in range(4):
            pos = arc.propagate(state, t - lag - pre, pre)[:3]
            ray = pos - arc.sites[i] - lag * velocities[i] + 0.5 * lag**2 * accelerations[i]
            lag = np.linalg.norm(ray) / o.LIGHT_SPEED
        directions.append(o.unit(ray))
    return np.array(directions)


def residual(p):
    u = directions(p)
    tangents = np.column_stack(((u * arc.east).sum(1), (u * arc.north).sum(1)))
    return (o.ASEC * tangents / (u * arc.obs).sum(1)[:, None]).ravel()


def fit_once(stencil_scale):
    steps = stencil_scale * np.array([1e-3, 2e-2, 2e-2, 1e-4, 2e-3, 2e-3])

    def jac(p):
        return np.column_stack([(residual(p + dv) - residual(p - dv)) / (2 * h)
                                for h, dv in zip(steps, np.diag(steps))])

    fit = least_squares(residual, np.zeros(6), jac=jac, x_scale="jac", xtol=1e-10,
                        ftol=1e-12, max_nfev=40)
    assert fit.success, fit.message
    return dict(stencil_scale=stencil_scale,
                boresight_shift_arcsec=float(o.distance(published, directions(fit.x, True)[-1])),
                tangent_rms_arcsec=float(np.sqrt(np.sum(fit.fun**2) / len(arc.rows))))


if __name__ == "__main__":
    print(json.dumps([fit_once(scale) for scale in (1.0, 0.5, 2.0)], indent=2))
