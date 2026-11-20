"""Build the pinned GRSS extension with Apple Clang using a local OpenMP installation."""

import argparse
import concurrent.futures
import subprocess
import sysconfig
import tempfile
from pathlib import Path
import pybind11

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("source", type=Path)
    p.add_argument("omp_header", type=Path)
    p.add_argument("omp_library", type=Path)
    args = p.parse_args()
    base = args.source.resolve()
    out = Path(sysconfig.get_path("purelib")) / "grss" / ("libgrss" + sysconfig.get_config_var("EXT_SUFFIX"))
    opts = [
        "clang++",
        "-std=c++17",
        "-O2",
        "-fPIC",
        "-Xpreprocessor",
        "-fopenmp",
        '-DGRSS_VERSION="4.6.0"',
        "-I" + str(base / "include"),
        "-I" + str(args.omp_header.resolve()),
        "-I" + sysconfig.get_path("include"),
        "-I" + pybind11.get_include(),
    ]
    with tempfile.TemporaryDirectory() as temp:

        def compile_one(src):
            obj = Path(temp) / (src.stem + ".o")
            subprocess.run([*opts, "-c", str(src), "-o", str(obj)], check=True)
            return str(obj)

        with concurrent.futures.ThreadPoolExecutor(4) as pool:
            objects = list(pool.map(compile_one, sorted((base / "src").glob("*.cpp"))))
        lib = args.omp_library.resolve()
        subprocess.run(
            [
                "clang++",
                "-bundle",
                "-undefined",
                "dynamic_lookup",
                *objects,
                str(lib),
                "-Wl,-rpath," + str(lib.parent),
                "-o",
                str(out),
            ],
            check=True,
        )
