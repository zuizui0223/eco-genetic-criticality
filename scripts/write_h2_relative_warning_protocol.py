"""Write the frozen H2-R proposition and calibration contract as JSON."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from causal_model.h2_relative_warning_contract import h2r_protocol_manifest


def main() -> int:
    parser = argparse.ArgumentParser(description="Write the H2-R relative-warning protocol manifest.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/h2_relative_warning/h2r_protocol_v1.json"),
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump(h2r_protocol_manifest(), handle, indent=2, sort_keys=True)
    print(f"Wrote H2-R protocol manifest: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
