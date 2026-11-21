"""Run the independent batch fits and check the inferred drift using pinned OEL."""

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from tempfile import mkdtemp

import numpy as np

ROOT = Path(__file__).resolve().parent
AU_KM = 149597870.7
DAY_S = 86400.0
MU = 0.01720209895**2 * AU_KM**3 / DAY_S**2


def energy_rate_drift(params, to_state):
    """Average d(a)/dt from v dot acceleration over one osculating ellipse."""
    a_km = params["q"] * AU_KM / (1 - params["e"])
    rates, weights = [], []
    for anomaly_deg in np.arange(0.0, 360.0, 360 / 512):
        r, v = to_state(
            a_km=a_km,
            ecc=params["e"],
            inc_deg=math.degrees(params["i"]),
            raan_deg=math.degrees(params["om"]),
            argp_deg=math.degrees(params["w"]),
            true_anomaly_deg=anomaly_deg,
            mu_km3_s2=MU,
        )
        radius = np.linalg.norm(r)
        h = np.cross(r, v)
        transverse = np.cross(h / np.linalg.norm(h), r / radius)
        accel = params["a2"] * AU_KM / DAY_S**2 * (radius / AU_KM) ** (-2.25)
        rates.append(2 * a_km**2 / MU * np.dot(v, transverse) * accel)
        weights.append(radius**2 / np.linalg.norm(h))
    return float(np.average(rates, weights=weights)) * 1000 * DAY_S * 365.25


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oel", default="oel")
    args = parser.parse_args()
    cmd = [args.oel, "--workspace", str(ROOT), "version"]
    version = json.loads(subprocess.check_output(cmd, text=True))
    if version["version"] != "1.1.0" or version["status"] != "ready":
        raise RuntimeError("The study requires a ready OEL 1.1.0 installation")
    inputs = ROOT / "inputs"
    provenance = json.loads((inputs / "provenance.json").read_text())
    for name, expected in provenance["normalized_sha256"].items():
        if hashlib.sha256((inputs / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Study input digest mismatch: {name}")
    sys.dont_write_bytecode = True
    source = Path(version["source"]).resolve()
    sys.path.insert(0, str(source))
    from sim.dynamics.orbit.elements import coe_to_rv_eci

    if not Path(sys.modules[coe_to_rv_eci.__module__].__file__).resolve().is_relative_to(source):
        raise RuntimeError("The state conversion must come from the selected OEL installation")
    outputs = ROOT / "outputs"
    outputs.mkdir(exist_ok=True)
    directory = Path(mkdtemp(prefix="replay-", dir=outputs))
    subprocess.run([sys.executable, str(ROOT / "bridge.py"), str(directory)], check=True)
    force = json.loads((directory / "force.json").read_text())
    null = json.loads((directory / "gravity.json").read_text())
    params = force["parameters"]
    drift = energy_rate_drift(params, coe_to_rv_eci)
    # Propagate covariance of eccentricity, perihelion distance and A2.
    # The orbit-averaged rate is independent of the orientation and perihelion epoch.
    gradient = np.zeros(len(force["parameter_order"]))
    for key, delta in [("e", 1e-6), ("q", 1e-6), ("a2", 1e-16)]:
        plus = energy_rate_drift({**params, key: params[key] + delta}, coe_to_rv_eci)
        minus = energy_rate_drift({**params, key: params[key] - delta}, coe_to_rv_eci)
        gradient[force["parameter_order"].index(key)] = (plus - minus) / (2 * delta)
    sigma = math.sqrt(float(gradient @ np.array(force["covariance"]) @ gradient))
    checks = {
        "orbit_average_agrees": abs(drift - force["drift_m_year"]) < 0.001,
        "replay_drift_agrees": abs(drift - (-282.9026)) < 0.1,
        "replay_uncertainty_agrees": abs(sigma - 1.48245) < 0.02,
        "gravity_only_fits_worse": null["chi2"] > force["chi2"],
        "publication_comparison": abs(drift - (-284.0)) < 3 * 1.5,
    }
    comparison = dict(
        oel=version["version"], grss="4.6.0", drift_m_year=drift, formal_sigma_m_year=sigma, checks=checks
    )
    (directory / "comparison.json").write_text(json.dumps(comparison, indent=2) + "\n")
    print(directory, json.dumps(comparison, indent=2))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
