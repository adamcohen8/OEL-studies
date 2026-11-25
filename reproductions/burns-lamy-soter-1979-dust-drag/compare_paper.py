"""Compare the offline results to the preselected circular-orbit benchmark."""
import json
from pathlib import Path

def main():
    results = json.loads(Path("results.json").read_text())
    fit = results["quadrature_fit"]
    expected = {"K_years": 400.0, "r0_au_exponent": 2.0, "beta_exponent": -1.0}
    errors = {key: fit[key] - value for key, value in expected.items()}
    relative = errors["K_years"] / expected["K_years"]
    accepted = (abs(relative) <= 0.02 and abs(errors["r0_au_exponent"]) <= 0.03
                and abs(errors["beta_exponent"]) <= 0.03)
    report = {"source": "Burns, Lamy and Soter (1979), Eq. (50), p. 31",
              "expected_zero_cutoff": expected,
              "measured_fit": {key: fit[key] for key in expected},
              "coefficient_relative_difference": relative,
              "tolerances": {"coefficient_relative": 0.02, "exponents_absolute": 0.03},
              "within_preselected_tolerances": accepted}
    Path("comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Eq. (50) coefficient difference {relative:+.3%}. Passed: {accepted}")
    if not accepted:
        raise SystemExit(1)

if __name__ == "__main__":
    main()
