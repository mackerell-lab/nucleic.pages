#!/usr/bin/env python3
"""Create or verify RNA-only dependencies without changing the DNA environment."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import venv

FR3D_URL = "https://github.com/BGSU-RNA/fr3d-python.git"
FR3D_COMMIT = "994e54ea8fcea1a0484ea8c082f4a41c5406191d"
HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[2]
DATA = WORKSPACE / "data" / "pure_rna"


def command(*args, **kwargs):
    return subprocess.check_output([str(arg) for arg in args], text=True, **kwargs).strip()


def verify(python, checkout):
    requirements = dict(line.split("==") for line in
                        (HERE / "requirements.txt").read_text().splitlines()
                        if line and not line.startswith("#"))
    versions = json.loads(command(python, "-c",
        "import importlib.metadata as m,json,sys; "
        "print(json.dumps({n:m.version(n) for n in sys.argv[1:]}))", *requirements))
    if versions != requirements:
        raise RuntimeError(f"RNA dependency versions differ: {versions}")
    revision = command("git", "-C", checkout, "rev-parse", "HEAD")
    if revision != FR3D_COMMIT:
        raise RuntimeError("FR3D revision differs from the tested pin")
    if command("git", "-C", checkout, "status", "--porcelain", "--untracked-files=no"):
        raise RuntimeError("FR3D tracked source contains changes")
    source_hashes = {str(path.relative_to(checkout)): hashlib.sha256(path.read_bytes()).hexdigest()
                     for path in sorted((checkout / "fr3d").rglob("*.py"))}
    return {"python": command(python, "--version"), "packages": versions,
            "requirements_sha256": hashlib.sha256((HERE / "requirements.txt").read_bytes()).hexdigest(),
            "fr3d": {"url": FR3D_URL, "commit": revision, "source_sha256": source_hashes}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Only verify existing dependencies")
    args = parser.parse_args()
    environment, checkout = DATA / "venv", DATA / "fr3d-python"
    python = environment / "bin" / "python"
    if not args.check:
        DATA.mkdir(parents=True, exist_ok=True)
        if not python.exists():
            venv.EnvBuilder(with_pip=True).create(environment)
        subprocess.check_call([str(python), "-m", "pip", "install", "-r", str(HERE / "requirements.txt")])
        if not checkout.exists():
            subprocess.check_call(["git", "clone", "--no-checkout", FR3D_URL, str(checkout)])
            subprocess.check_call(["git", "-C", str(checkout), "checkout", "--detach", FR3D_COMMIT])
    report = verify(python, checkout)
    target = DATA / "dependency_lock.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"ok": True, "python": report["python"], "packages": report["packages"],
                      "fr3d_commit": report["fr3d"]["commit"], "report": str(target)}, indent=2))


if __name__ == "__main__":
    main()
