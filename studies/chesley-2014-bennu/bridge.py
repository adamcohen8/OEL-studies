"""Independent reconstruction of Bennu's historical arc using pinned GRSS."""

import argparse
import json
import math
from pathlib import Path
import sysconfig
from unittest.mock import patch

import numpy as np
import pandas as pd
import requests
from astropy.time import Time
from astropy.utils import iers

ROOT = Path(__file__).resolve().parent


def fit_arc(directory, force, grss_fit, grss_columns):
    data = json.loads((ROOT / "inputs" / "observations.json").read_text())
    rows = []
    for optical in data["optical"]:
        if optical["use"]:
            rows.append(
                {
                    **optical,
                    "sigTime": 0.0,
                    "sigCorr": 0.0,
                    "selAst": "a",
                    "cosDec": math.cos(math.radians(optical["dec"])),
                }
            )
    for obs in data["radar"]:
        row = dict(
            permID="101955",
            mode="RAD",
            obsTime=obs["epoch"].replace(" ", "T") + "Z",
            rcv=obs["rcvr"],
            trx=obs["xmit"],
            frq=float(obs["freq"]),
            com=1,
            selAst="a",
        )
        value, sigma = float(obs["value"]), float(obs["sigma"])
        if obs["units"] == "us":
            row.update(delay=value / 1e6, sigDelay=sigma)
        else:
            row.update(doppler=value, sigDoppler=sigma)
        rows.append(row)
    df = pd.DataFrame(rows)
    for col in grss_columns:
        if col not in df:
            df[col] = np.nan
    times = Time(list(df.obsTime), format="isot", scale="utc")
    df["obsTimeMJD"], df["obsTimeMJDTDB"] = times.utc.mjd, times.tdb.mjd
    # Rounded initial geometry. The transverse coefficient starts at zero.
    seed = dict(
        t=55562.0,
        e=0.203745,
        q=0.896894,
        tp=Time("2010-08-30", scale="tdb").mjd + 0.642,
        om=math.radians(2.061),
        w=math.radians(66.223),
        i=math.radians(6.035),
    )
    if force:
        seed["a2"] = 0.0
    size = len(seed) - 1
    sim = grss_fit.FitSimulation(seed, df, np.eye(size), de_kernel=440, n_iter_max=12, nongrav_info={"m": 2.25})
    sim.analytic_partials = True
    sim.reject_outliers = False
    sim.filter_lsq()
    if not sim.converged:
        raise RuntimeError("GRSS fit did not converge")
    label = "force" if force else "gravity"
    res = sim.obs
    cov = np.asarray(sim.covariance)
    if not np.all(np.isfinite(cov)) or not np.all(np.diag(cov) > 0):
        raise ValueError("Invalid covariance")
    report = {
        "model": label,
        "parameters": sim.x_nom,
        "covariance": cov.tolist(),
        "parameter_order": list(sim.x_nom),
        "optical": len(sim.optical_idx),
        "delay": len(sim.delay_idx),
        "doppler": len(sim.doppler_idx),
    }
    # Count and validate every weighted component before comparing the nested models.
    normalized = []
    for field, indices, expected in [
        ("RA", sim.optical_idx, 561),
        ("Dec", sim.optical_idx, 561),
        ("Delay", sim.delay_idx, 22),
        ("Doppler", sim.doppler_idx, 7),
    ]:
        selected = res.iloc[indices]
        sigma = selected["sig" + field].to_numpy(dtype=float)
        residual = selected["res" + field].to_numpy(dtype=float)
        if len(selected) != expected or not (
            np.all(np.isfinite(sigma)) and np.all(sigma > 0) and np.all(np.isfinite(residual))
        ):
            raise ValueError("Invalid residual coverage for " + field)
        normalized.extend(residual / sigma)
    report["chi2"] = float(np.dot(normalized, normalized))
    if force:
        # Orbit average of Gauss's semimajor-axis equation, uniform in eccentric anomaly with time weight.
        e = sim.x_nom["e"]
        a = sim.x_nom["q"] / (1 - e)
        anomaly = 2 * np.pi * np.arange(4096) / 4096
        n = 0.01720209895 / a**1.5
        scale = (
            2
            * math.sqrt(1 - e * e)
            / n
            * a ** (-2.25)
            * np.mean((1 - e * np.cos(anomaly)) ** (-2.25))
            * 149597870700
            * 365.25
        )
        report.update(
            drift_m_year=scale * sim.x_nom["a2"], conditional_drift_sigma_m_year=scale * math.sqrt(cov[-1, -1])
        )
    (directory / f"{label}.json").write_text(json.dumps(report, indent=2) + "\n")
    res[
        [
            "mode",
            "stn",
            "obsTime",
            "selAst",
            "resRA",
            "resDec",
            "resDelay",
            "resDoppler",
            "sigRA",
            "sigDec",
            "sigDelay",
            "sigDoppler",
        ]
    ].to_csv(directory / f"{label}-residuals.csv", index=False)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    # This environment is dedicated to the study. Load the frozen MPC station catalog.
    import sys

    if sys.prefix == sys.base_prefix:
        raise SystemExit("Use the dedicated study virtual environment.")
    target = Path(sysconfig.get_path("purelib")) / "grss" / "fit"
    codes = json.loads((ROOT / "inputs" / "stations.json").read_text())
    special = {"gaia": ["258"], "occultation": [], "roving": [], "spacecraft": []}
    (target / "codes.json").write_text(json.dumps(codes))
    (target / "special_codes.json").write_text(json.dumps(special))
    iers.conf.auto_download = False
    with patch("requests.head", side_effect=requests.ConnectionError("offline pinned inputs")):
        from grss import fit
        from grss.fit.fit_ades import ades_column_types
    args.output.mkdir(parents=True, exist_ok=True)
    for force in (True, False):
        report = fit_arc(args.output, force, fit, ades_column_types)
        print(
            json.dumps({k: v for k, v in report.items() if k not in ("parameters", "covariance", "parameter_order")}),
            flush=True,
        )
