"""Check the editable public OEL installation and run the study from this folder."""
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys

import sim
from oel_rust_orbit import PlanetaryPlan


def main():
    workspace = Path(__file__).resolve().parent
    lock = json.loads((workspace / "oel-revision.lock.json").read_text())
    checkout = Path(sim.__file__).resolve().parent.parent
    revision = subprocess.check_output(["git", "-C", str(checkout), "rev-parse", "HEAD"], text=True).strip()
    if revision != lock["public_commit"] or version("orbital-engineering-lab") != lock["version"]:
        raise SystemExit("Install the exact public OEL checkout in editable mode as described in README.md")
    subprocess.run(["git", "-C", str(checkout), "diff", "--quiet", "HEAD", "--"], check=True)
    if not getattr(PlanetaryPlan, "supports_schwarzschild", False):
        raise SystemExit("Build the native orbit extension from the pinned checkout first")
    return subprocess.call([sys.executable, str(workspace / "mercury_apsidal_convergence.py"),
                            "--output", str(workspace / "outputs" / "mercury.json")], cwd=workspace)


if __name__ == "__main__":
    raise SystemExit(main())
