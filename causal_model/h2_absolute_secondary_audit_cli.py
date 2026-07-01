"""CLI for a no-resimulation H2-A fixed-threshold secondary audit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from causal_model.h2_absolute_secondary_audit import (
    audit_h2a_from_h2r_validation_payload,
    write_h2a_secondary_audit_artifacts,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Apply H2-A's existing fixed 0.20 thresholds to an H2-R validation JSON artifact."
    )
    parser.add_argument("--input", type=Path, required=True, help="independent_relative_warning_v1.json artifact path")
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/h2a_secondary_audit"))
    parser.add_argument("--prefix", default="h2a_fixed_threshold_secondary_audit_v1")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    with args.input.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    audit = audit_h2a_from_h2r_validation_payload(payload)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_h2a_secondary_audit_artifacts(
        audit,
        csv_path=args.output_dir / f"{args.prefix}.csv",
        json_path=args.output_dir / f"{args.prefix}.json",
    )
    print(f"Wrote H2-A secondary audit artifacts to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
