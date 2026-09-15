#!/usr/bin/env python3
"""Launch the local workbench as a detached, health-checked process."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request
import webbrowser
from pathlib import Path


GUI_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = GUI_DIR.parent
SERVER_PATH = GUI_DIR / "server.py"
PID_PATH = GUI_DIR / ".workbench.pid"
STDOUT_PATH = GUI_DIR / "server.stdout.log"
STDERR_PATH = GUI_DIR / "server.stderr.log"


def workbench_url(host: str, port: int) -> str:
    return f"http://{host}:{port}"


def is_healthy(url: str, timeout: float = 5.0) -> bool:
    """Return true only when the workbench configuration API responds."""
    try:
        with urllib.request.urlopen(f"{url}/api/v2/config", timeout=timeout) as response:
            if response.status != 200:
                return False
            payload = json.load(response)
    except (OSError, ValueError, urllib.error.URLError):
        return False
    return isinstance(payload, dict) and "defaults" in payload and "registry" in payload


def detached_creation_options() -> dict[str, object]:
    if os.name == "nt":
        return {
            "creationflags": (
                subprocess.CREATE_NEW_PROCESS_GROUP
                | subprocess.CREATE_NO_WINDOW
                | subprocess.DETACHED_PROCESS
            )
        }
    return {"start_new_session": True}


def launch(host: str, port: int, timeout: float) -> tuple[int, bool]:
    """Return (pid, started_now); reuse an already healthy local server."""
    url = workbench_url(host, port)
    if is_healthy(url):
        stored_pid = PID_PATH.read_text(encoding="utf-8").strip() if PID_PATH.exists() else ""
        return (int(stored_pid) if stored_pid.isdigit() else 0), False

    command = [
        sys.executable,
        str(SERVER_PATH),
        "--host",
        host,
        "--port",
        str(port),
        "--no-browser",
    ]
    with STDOUT_PATH.open("ab", buffering=0) as stdout, STDERR_PATH.open(
        "ab", buffering=0
    ) as stderr:
        process = subprocess.Popen(  # noqa: S603 - fixed local command
            command,
            cwd=PROJECT_ROOT,
            stdin=subprocess.DEVNULL,
            stdout=stdout,
            stderr=stderr,
            close_fds=True,
            **detached_creation_options(),
        )
    PID_PATH.write_text(f"{process.pid}\n", encoding="utf-8")

    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        remaining = deadline - time.monotonic()
        if is_healthy(url, timeout=max(0.1, min(5.0, remaining))):
            return process.pid, True
        if process.poll() is not None:
            break
        time.sleep(0.1)
    raise RuntimeError(
        "The workbench server did not become healthy. Review "
        f"{STDERR_PATH} for details."
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--no-browser", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    url = workbench_url(args.host, args.port)
    pid, started_now = launch(args.host, args.port, args.timeout)
    status = "started" if started_now else "already running"
    pid_text = f" (PID {pid})" if pid else ""
    print(f"Graphical Bayesian Inference workbench {status}{pid_text}: {url}")
    if not args.no_browser:
        webbrowser.open(f"{url}/v2.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
