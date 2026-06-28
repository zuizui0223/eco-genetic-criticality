#!/usr/bin/env python3
"""Run paired trait/genetic/full baseline comparisons from the repository root.

Example:
    python scripts/run_paired_baseline_comparisons.py --profile standard --output-dir artifacts/baseline_comparison/standard
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from causal_model.paired_baseline_cli import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
