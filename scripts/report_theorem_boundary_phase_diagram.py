#!/usr/bin/env python3
"""Generate figure-ready reports from a theorem-boundary CSV artifact.

Example:
    python scripts/report_theorem_boundary_phase_diagram.py artifacts/run.csv --output-dir artifacts/report
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from causal_model.theorem_boundary_report import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
