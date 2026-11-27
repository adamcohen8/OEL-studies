"""Verify the installed editable public OEL source revision used by this workspace."""
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
from sim.dynamics.orbit import elements
from sim.dynamics.attitude import disturbances

def verify():
    lock = json.loads(Path(__file__).with_name("oel-revision.lock.json").read_text())
    root = Path(elements.__file__).resolve().parents[3]
    def git(*args):
        return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()
    if Path(git("rev-parse", "--show-toplevel")).resolve() != root:
        raise RuntimeError("Install OEL in editable mode from the pinned public source checkout.")
    if not Path(disturbances.__file__).resolve().is_relative_to(root):
        raise RuntimeError("The OEL torque and orbit routines came from different source trees.")
    if git("rev-parse", "HEAD") != lock["public_commit"] or git("status", "--porcelain", "--untracked-files=no"):
        raise RuntimeError("The OEL source does not match the clean public study pin.")
    for name, expected in [("orbital-engineering-lab", lock["version"]), ("oel-rust-orbit", lock["native_orbit_version"])]:
        if version(name) != expected:
            raise RuntimeError(f"Expected {name} {expected}, found {version(name)}")
    return lock["public_commit"]

if __name__ == "__main__":
    print(verify())
