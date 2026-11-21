"""Evaluate the 13 earlier observations with the paper's nominal OEL force model.

Run with the study's installed Python environment. Output is CSV on stdout.
The last row, 75, repeats the first training observation as a time/frame check.
"""
import csv
import html
import json
import re
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import root

STUDY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STUDY / "field-certificate"))
from common import Arc, forecast, o


def load(path):
    return json.loads((STUDY / path).read_text())


def earlier_observations(arc):
    # IERS finals2000A, 2024-08-14, columns MJD, xp, yp, UT1-UTC.
    eop = [60536.0, 0.187692, 0.466829, 0.0379662]
    arc.eop = np.insert(arc.eop, np.searchsorted(arc.eop[:, 0], eop[0]), eop, axis=0)
    # MPC observatory-code table: G17 longitude, rho*cos(phi'), rho*sin(phi').
    lon = np.deg2rad(11.20261)
    arc.fixed["G17"] = 6378.14 * np.array([np.cos(lon) * 0.723508, np.sin(lon) * 0.723508, 0.688023])
    raw = (STUDY / "followup/corrected/missions/juice.html").read_text()
    rows = []
    for match in re.finditer(r'<a name="o(\d+)"></a>([^\n]+)', raw):
        index = int(match[1])
        if not 62 <= index <= 75:
            continue
        line = html.unescape(re.sub("<[^>]+>", "", match[2]))
        stamp, code = line[15:32].strip(), line[-3:]
        t = forecast.parse_time(stamp) - arc.rows[0][0]
        ra, dec = np.deg2rad([15 * o.sexagesimal(line[32:44]), o.sexagesimal(line[44:56])])
        direction = np.array([np.cos(ra) * np.cos(dec), np.sin(ra) * np.cos(dec), np.sin(dec)])
        rows.append((index, code, stamp, t, arc.site(code, t), direction))
    assert [row[0] for row in rows] == list(range(62, 76))
    return rows


def residuals(arc, rows, state):
    state = np.asarray(state)
    start = rows[0][3] - 20.0
    seed = state.copy()
    seed[:3] += start * state[3:]
    scale = np.array([5000.0] * 3 + [0.05] * 3)

    def closure(p):
        return (arc.propagate(seed + scale * p, -start, start) - state) / scale

    fit = root(closure, np.zeros(6), options={"eps": 1e-9, "xtol": 1e-11})
    error = closure(fit.x) * scale
    assert np.linalg.norm(error[:3]) < 1e-5 and np.linalg.norm(error[3:]) < 1e-9
    early_state = seed + scale * fit.x

    def earth(t):
        return o.hpop_de440_positions_km(arc.JD + t / 86400, arc.ephem)["earth"]

    angles = []
    for _, _, _, t, site, direction in rows:
        lag = 5.0
        for _ in range(5):
            emission = t - lag
            position = arc.propagate(early_state, emission - start, start)[:3]
            ray = position - site + (earth(emission) - earth(t))
            updated_lag = np.linalg.norm(ray) / o.LIGHT_SPEED
            change = abs(lag - updated_lag)
            lag = updated_lag
        assert change < 1e-8
        angles.append(float(o.distance(o.unit(ray), direction)))
    return angles


if __name__ == "__main__":
    arc = Arc("juice")
    observations = earlier_observations(arc)
    trio = load("square-decision/threshold-certificate.json")["2.15"]["accepted"]
    states = {
        "LS": load("followup/corrected/missions/juice-forecast.json")["state"],
        "missed": load("square-decision/nominal-certificate.json")["2"]["accepted"][5]["state"],
        **{label: trio[i]["state"] for label, i in zip("ABC", (1, 3, 4))},
    }
    scores = [residuals(arc, observations, state) for state in states.values()]
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(["source_row", "station", "utc_fractional_day", *[key + "_arcsec" for key in states]])
    for row, values in zip(observations, zip(*scores)):
        writer.writerow([*row[:3], *[f"{v:.5f}" for v in values]])
