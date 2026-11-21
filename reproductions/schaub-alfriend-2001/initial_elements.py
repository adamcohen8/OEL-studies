"""Optional input audit. Run with Basilisk 2.12.1 to print the prepared elements."""
import json
from math import degrees, radians

from Basilisk.utilities import orbitalMotion


def initial_elements(perigee):
    states = {}
    for member in ("chief", "deputy"):
        delta = int(member == "deputy")
        elements = orbitalMotion.ClassicElements()
        elements.a = 7153 - delta * 0.000351765
        elements.e = 0.05 + delta * 0.0001
        elements.i = radians(48 + delta * 0.001035)
        elements.Omega = radians(delta * 0.005)
        elements.omega = radians(perigee + delta * 0.01)
        anomaly = orbitalMotion.M2E(radians(-delta * 0.01), elements.e)
        elements.f = orbitalMotion.E2f(anomaly, elements.e)
        mapped = orbitalMotion.ClassicElements()
        orbitalMotion.clMeanOscMap(6378.1366, 0.001082616, elements, mapped, 1)
        for method, state in (("direct", elements), ("mapped", mapped)):
            states[f"{method}_{member}"] = {
                "a_km": state.a, "ecc": state.e, "inc_deg": degrees(state.i),
                "raan_deg": degrees(state.Omega) % 360,
                "argp_deg": degrees(state.omega) % 360,
                "true_anomaly_deg": degrees(state.f) % 360,
            }
    return states


if __name__ == "__main__":
    print(json.dumps({angle: initial_elements(angle) for angle in (30, 75)}, indent=2))
