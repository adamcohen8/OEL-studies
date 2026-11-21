"""Regenerate four nominal witness enclosures before testing detector geometry."""
import json
import sys
from pathlib import Path
from types import FunctionType

import numpy as np

STUDY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(STUDY / "square-decision"), str(STUDY / "field-certificate")]
from common import Arc, o
from certify_pairs import certify
from interval_roll import run
from tubes import prepare


def load(path):
    return json.loads((STUDY / path).read_text())


def fine_knots(arc, state, scale, knots):
    return prepare(arc, state, scale / 3, knots)


def recompute(arc, state, cap):
    cert = certify(arc, state, cap)
    if cert["status"] != "bounded":
        # Same documented finer knot retry used by square-decision/verify.py.
        finer = FunctionType(certify.__code__, dict(certify.__globals__, prepare=fine_knots),
                             argdefs=certify.__defaults__)
        cert = finer(arc, state, cap)
    assert cert["accepted"] and cert["min_moon_floor_km"] > 1737.4
    assert cert["training_moon_floor_km"] > 1e5 and cert["min_earth_floor_km"] > 2e5
    assert cert["max_speed_upper"] < 6
    return cert


if __name__ == "__main__":
    arc = Arc("juice")
    low = load("square-decision/nominal-certificate.json")["2"]["accepted"]
    high = load("square-decision/threshold-certificate.json")["2.15"]["accepted"]
    missed = recompute(arc, low[5]["state"], 2.0)
    trio = [recompute(arc, high[i]["state"], 2.148) for i in (1, 3, 4)]
    ls = np.array(load("followup/corrected/missions/juice-score.json")["replay"][0]["los"])
    corner = np.rad2deg(np.arctan(np.sqrt(2) * np.tan(np.deg2rad(2.2))))
    sep = o.distance(ls, np.array(missed["forecast_center"])) - missed["forecast_radius_as"]
    miss_arcmin = (sep - corner * 3600) / 60
    roll = run([{"validation": cert} for cert in trio], [0, 1, 2])
    assert miss_arcmin > 0 and roll["success"]
    print(json.dumps(dict(miss_arcmin=miss_arcmin, roll=roll, witnesses=dict(zip(
        ["missed", "A", "B", "C"], [missed, *trio]))), indent=2))
