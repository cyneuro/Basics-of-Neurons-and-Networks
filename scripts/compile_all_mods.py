#!/usr/bin/env python3
"""Compile NEURON .mod files in every tutorial folder that ships mechanisms."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MOD_DIRS = [
    ROOT / "S2_ActionPotential",
    ROOT / "S3_Burster",
    ROOT / "S4_Synapses",
    ROOT / "S5_CPG",
    ROOT / "S6_STM_WTA",
    ROOT / "B3_MoreProperties" / "mod",
]

BUILD_DIR_NAMES = ("x86_64", "arm64", "aarch64")


def clear_build_dirs(directory: Path) -> None:
    for name in BUILD_DIR_NAMES:
        build_dir = directory / name
        if build_dir.is_dir():
            shutil.rmtree(build_dir)
            print(f"  removed {build_dir.relative_to(ROOT)}")


def compile_mods(directory: Path) -> int:
    print(f"\n=== {directory.relative_to(ROOT)} ===")
    if not directory.is_dir():
        print("  ERROR: directory missing")
        return 1
    if not any(directory.glob("*.mod")):
        print("  ERROR: no .mod files found")
        return 1

    clear_build_dirs(directory)
    result = subprocess.run(
        ["nrnivmodl"],
        cwd=directory,
        capture_output=True,
        text=True,
    )
    if result.stdout.strip():
        print(result.stdout.rstrip())
    if result.returncode != 0:
        if result.stderr.strip():
            print(result.stderr.rstrip(), file=sys.stderr)
        print(f"  FAILED (return code {result.returncode})")
        return result.returncode

    print("  OK")
    return 0


def main() -> int:
    failures: list[str] = []
    for mod_dir in MOD_DIRS:
        code = compile_mods(mod_dir)
        if code != 0:
            failures.append(str(mod_dir.relative_to(ROOT)))

    print()
    if failures:
        print("Compilation failed in:")
        for path in failures:
            print(f"  - {path}")
        return 1

    print("All mechanism folders compiled successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
