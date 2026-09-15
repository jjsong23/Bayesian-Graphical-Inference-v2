#!/usr/bin/env python3
"""Portable root entry point for the Version 2 backward frontier workflow."""

import sys
from pathlib import Path


CODE_DIRECTORY = Path(__file__).resolve().parent / "code"
if str(CODE_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CODE_DIRECTORY))

from backward_search.cli import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
