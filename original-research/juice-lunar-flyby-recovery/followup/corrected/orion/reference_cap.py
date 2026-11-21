"""Post-hoc sensitivity to caps spanning the maximum archive-to-Horizons difference."""
import json
import orion as o
from outer import enclose

fit = json.loads((o.HERE/"orion-results.json").read_text())["triple"]
results = {}
for cap in (6., 11.):
    result = enclose(20, fit, refine=True, angular_cap_as=cap)
    results[str(cap)] = 2 * result["sufficient_field_radius_deg"]
    print(cap, results[str(cap)], flush=True)
(o.HERE / "reference-cap.json").write_text(json.dumps(results, indent=2) + "\n")
