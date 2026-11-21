"""Independently score the native TEME states against the archived reference."""
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
criteria = json.loads((ROOT / "criteria.json").read_text())
reference = ROOT / "00005.e"
assert hashlib.sha256(reference.read_bytes()).hexdigest() == criteria["source_sha256"]
data = reference.read_text().split("EphemerisTimePosVel")[1].split("END Ephemeris")[0]
expected = [list(map(float, line.split())) for line in data.splitlines() if line.strip()]
states_path = Path(sys.argv[1]) if len(sys.argv) == 2 else ROOT / "states.csv"
with states_path.open() as stream:
    reader = csv.reader(stream)
    next(reader)
    states = [list(map(float, row)) for row in reader]
epochs = criteria["epochs_min"]
if [r[0] for r in states] != epochs or [r[0] / 60 for r in expected] != epochs:
    raise ValueError("The reference and calculated epoch lists must match the frozen grid.")
if any(len(row) != 7 or not all(map(math.isfinite, row)) for row in [*states, *expected]):
    raise ValueError("Expected six finite TEME coordinates and one time per row.")
dr = [math.dist(a[1:4], b[1:4]) for a,b in zip(states,expected)]
dv = [math.dist(a[4:], b[4:]) for a,b in zip(states,expected)]
passed = max(dr) <= criteria["max_position_norm_km"] and max(dv) <= criteria["max_velocity_norm_km_s"]
print(json.dumps(dict(epochs=len(states),max_position_error_km=max(dr),
    max_position_error_at_min=epochs[dr.index(max(dr))], max_velocity_error_km_s=max(dv),
    max_velocity_error_at_min=epochs[dv.index(max(dv))], passed=passed),indent=2))
raise SystemExit(0 if passed else 1)
