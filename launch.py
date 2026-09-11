#!/usr/bin/env python3
"""Launch or validate the Graphical Bayesian Inference workbench.

Typical use after installing ``requirements.txt``::

    python launch.py

The default launch is detached and opens the local web interface. Use
``--foreground`` for a server that remains attached to the terminal and stops
with Ctrl+C, or ``--check`` to validate a newly extracted release without
starting a server.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent
GUI_DIR = PROJECT_ROOT / "gui"
CORE_IMPORTS = {
    "numpy": "numpy",
    "pandas": "pandas",
    "openpyxl": "openpyxl",
}
FEATURE_IMPORTS = {
    "scipy (parameter calibration)": "scipy",
    "Pillow": "PIL",
    "networkx": "networkx",
    "obonet": "obonet",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--no-browser", action="store_true")
    parser.add_argument(
        "--foreground",
        action="store_true",
        help="Keep the server attached to this terminal instead of launching it in the background.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Validate dependencies, registry entries, and packaged inputs without launching.",
    )
    return parser.parse_args()


def missing_dependencies(packages: dict[str, str]) -> list[str]:
    return [
        package
        for package, module in packages.items()
        if importlib.util.find_spec(module) is None
    ]


def validate_release() -> dict[str, object]:
    missing = missing_dependencies(CORE_IMPORTS)
    if missing:
        joined = ", ".join(missing)
        raise RuntimeError(
            f"Missing Python dependencies: {joined}. Run: "
            f"{sys.executable} -m pip install -r requirements.txt"
        )
    if str(GUI_DIR) not in sys.path:
        sys.path.insert(0, str(GUI_DIR))
    from workflow_engine import (  # noqa: PLC0415
        UNIVERSE_RELATIVE,
        default_configuration,
        load_node_factor_catalog,
        load_registry,
    )

    registry = load_registry(PROJECT_ROOT)
    config = default_configuration(registry)
    factors = load_node_factor_catalog(PROJECT_ROOT, registry)
    universe_path = PROJECT_ROOT / UNIVERSE_RELATIVE
    if not universe_path.is_file():
        raise FileNotFoundError(f"Seed universe is missing: {universe_path}")
    return {
        "status": "ok",
        "project_root": str(PROJECT_ROOT),
        "python": sys.version.split()[0],
        "missing_optional_dependencies": missing_dependencies(FEATURE_IMPORTS),
        "node_streams": len(registry["node_streams"]),
        "edge_streams": len(registry["edge_streams"]),
        "node_candidates": len(factors),
        "default_node_streams": sum(
            bool(state["enabled"]) for state in config["node_streams"].values()
        ),
        "default_edge_streams": sum(
            bool(state["enabled"]) for state in config["edge_streams"].values()
        ),
    }


def launch(args: argparse.Namespace) -> int:
    if args.foreground:
        command = [
            sys.executable,
            str(GUI_DIR / "server.py"),
            "--host",
            args.host,
            "--port",
            str(args.port),
        ]
        if args.no_browser:
            command.append("--no-browser")
    else:
        command = [
            sys.executable,
            str(GUI_DIR / "launch_workbench.py"),
            "--host",
            args.host,
            "--port",
            str(args.port),
            "--timeout",
            str(args.timeout),
        ]
        if args.no_browser:
            command.append("--no-browser")
    return subprocess.call(command, cwd=PROJECT_ROOT)


def main() -> int:
    args = parse_args()
    try:
        report = validate_release()
    except Exception as exc:  # noqa: BLE001 - command-line boundary
        print(f"Launch validation failed: {exc}", file=sys.stderr)
        return 2
    if args.check:
        print(json.dumps(report, indent=2))
        return 0
    return launch(args)


if __name__ == "__main__":
    raise SystemExit(main())
