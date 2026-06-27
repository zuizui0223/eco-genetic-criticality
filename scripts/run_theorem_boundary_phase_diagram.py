#!/usr/bin/env python3
"""Run a theorem-boundary phase diagram from the repository root.

Example:
    python scripts/run_theorem_boundary_phase_diagram.py --profile standard --output-dir artifacts/standard
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from causal_model.theorem_boundary_cli import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
