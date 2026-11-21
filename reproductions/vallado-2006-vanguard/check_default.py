"""Post-submission diagnostic: grade OEL 1.1.0's default propagation backend."""
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

from sim.dynamics.orbit.ogp import ogp_propagate_teme
from sim.dynamics.orbit.tle import parse_tle_lines

root = Path(__file__).resolve().parent
elements = parse_tle_lines(*(root / "input.tle").read_text().splitlines(), require_checksum=True)
with tempfile.TemporaryDirectory() as temporary:
    output = Path(temporary) / "default.csv"
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["t_min","x_km","y_km","z_km","vx_km_s","vy_km_s","vz_km_s"])
        for minutes in range(0, 4321, 360):
            state = ogp_propagate_teme(elements, minutes)
            if state.error is not None:
                raise ValueError(state.error)
            writer.writerow([minutes, *state.position_teme_km, *state.velocity_teme_km_s])
    raise SystemExit(subprocess.run([sys.executable, "-I", str(root / "grade.py"), str(output)]).returncode)
