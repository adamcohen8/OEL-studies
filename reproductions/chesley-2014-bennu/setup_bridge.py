"""Prepare a dedicated CPython 3.11 environment for the verified Apple Silicon replay."""

import concurrent.futures
import hashlib
import platform
from pathlib import Path
import shutil
import subprocess
import sys
import sysconfig
import tarfile
import tempfile

from build_grss import __file__ as build_script

ROOT = Path(__file__).resolve().parent
SOURCE = "b0bc1124795493ad88589af2988ff0502cacc188"
PCK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/"
FILES = {
    "de440.bsp": (
        "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440.bsp",
        "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5",
    ),
    "sb441-n16s.bsp": (
        "https://ssd.jpl.nasa.gov/ftp/xfr/sb441-n16s.bsp",
        "815981128c8ecf06547c4562a9b2ae76e4ae17b9a8831743afdd74ea3f823ecf",
    ),
    "earth_historic.bpc": (
        PCK + "earth_620120_260806.bpc",
        "652b02c0b3a5a283d393497f7f96e1fe4aa47ee6dcc035a089eb9e6b5df86643",
    ),
    "earth_predict.bpc": (
        PCK + "earth_2026_260806_2126_predict.bpc",
        "9f565ab5aaccd33afe54f3a1b20cb8b136636c31a0bd4c0986fac722b931c280",
    ),
    "moon_pa_de440.bpc": (
        PCK + "moon_pa_de440_200625.bpc",
        "60cd55aa401ea2ea97360636f567554bfe4e37bb829f901b4460a455dfaf783f",
    ),
}


def download(path, url, expected):
    def digest(p):
        with p.open("rb") as stream:
            return hashlib.file_digest(stream, "sha256").hexdigest()

    if path.is_file() and digest(path) == expected:
        return
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as f:
        part = Path(f.name)
    try:
        subprocess.run(["curl", "-fLsS", "--retry", "2", "--max-time", "300", "-o", str(part), url], check=True)
        if digest(part) != expected:
            raise ValueError("Download digest mismatch for " + path.name)
        part.replace(path)
    finally:
        part.unlink(missing_ok=True)


def main():
    if (platform.system(), platform.machine(), sys.version_info[:2]) != ("Darwin", "arm64", (3, 11)):
        raise SystemExit("The tested bridge build requires macOS arm64 and CPython 3.11.")
    if sys.prefix == sys.base_prefix:
        raise SystemExit("Activate a dedicated study virtual environment first.")
    packages = Path(sysconfig.get_path("purelib"))
    kernels = packages / "grss" / "kernels"
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        tasks = [pool.submit(download, kernels / name, url, digest) for name, (url, digest) in FILES.items()]
        for task in tasks:
            task.result()
    # GRSS opens three Earth PCK handles. The historical file covers every date
    # in this arc, so keep its otherwise-unused "latest" handle immutable too.
    shutil.copyfile(kernels / "earth_historic.bpc", kernels / "earth_latest.bpc")
    with tempfile.TemporaryDirectory() as t:
        temp = Path(t)
        archive = temp / "openmp.tar.bz2"
        download(
            archive,
            "https://conda.anaconda.org/conda-forge/osx-arm64/llvm-openmp-15.0.5-h7cfbb63_0.tar.bz2",
            "df8d545338b3a895d5d010aa1d672d484006b4f3f406ff5146b18f51f3acd0ea",
        )
        with tarfile.open(archive) as tf:
            (temp / "omp.h").write_bytes(tf.extractfile("include/omp.h").read())
        src = temp / "grss"
        subprocess.run(["git", "init", "-q", str(src)], check=True)
        subprocess.run(
            ["git", "-C", str(src), "fetch", "--depth", "1", "https://github.com/rahil-makadia/grss.git", SOURCE],
            check=True,
        )
        subprocess.run(["git", "-C", str(src), "checkout", "-q", "--detach", "FETCH_HEAD"], check=True)
        if subprocess.check_output(["git", "-C", str(src), "rev-parse", "HEAD"], text=True).strip() != SOURCE:
            raise RuntimeError("GRSS source revision mismatch")
        subprocess.run(
            [sys.executable, build_script, str(src), str(temp), str(packages / "healpy" / ".dylibs" / "libomp.dylib")],
            check=True,
        )


if __name__ == "__main__":
    main()
