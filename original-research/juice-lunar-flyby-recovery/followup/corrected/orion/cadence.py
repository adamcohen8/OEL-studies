"""Locate the first sampled California data cut that passes the example field test."""

import json
import numpy as np
import orion as o
from outer import enclose

data = json.loads((o.HERE / "orion-results.json").read_text())
seed = np.array(data["double"]["baseline"]["state_eci_km_km_s"])
rows = []
for n in range(9, 20):
    opt = o.least_squares(
        lambda p: o.residual(o.sky(seed + o.correction(p), n), n).ravel(), np.zeros(6), jac="3-point", max_nfev=45
    )
    x = seed + o.correction(opt.x)
    fit = {"baseline": {"state_eci_km_km_s": x.tolist()}}
    upper = enclose(n, fit, refine=True)
    rows.append(
        dict(
            training_count=n,
            california_count=n - 8,
            last_row=data["rows"][n - 1],
            fit_max_as=float(max(o.distance(o.sky(x, n), o.obs[:n]))),
            first_holdout_as=float(o.distance(o.sky(x, 21)[20], o.obs[20])),
            upper_diameter_deg=2 * upper["sufficient_field_radius_deg"],
        )
    )
    print(rows[-1], flush=True)
(o.HERE / "cadence.json").write_text(json.dumps(rows, indent=2) + "\n")
