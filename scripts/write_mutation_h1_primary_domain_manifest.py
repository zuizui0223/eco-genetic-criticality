"""Write the frozen mutation-H1 primary-analysis domain manifest."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from causal_model.mutation_h1_primary_domain import domain_manifest


def main() -> int:
    parser = argparse.ArgumentParser(description="Write the independent-validation mutation-H1 primary-analysis domain.")
    parser.add_argument("--output", type=Path, default=Path("artifacts/mutation_h1_primary_domain/primary_domain_v1.json"))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump(domain_manifest(), handle, indent=2, sort_keys=True)
    print(f"Wrote primary domain manifest: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
