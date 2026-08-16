"""CLI for the preregistered H3 fragmentation-gradient sensitivity."""
from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path
from typing import Sequence

from causal_model.h3_fragmentation_gradient_sensitivity import (
    FRAGMENT_PATCH_COUNTS,
    FRESH_GRADIENT_MASTER_SEEDS,
    run_h3_fragmentation_gradient_sensitivity,
    write_fragmentation_gradient_artifacts,
)
from causal_model.multipatch_criticality_experiments import standard_profile


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description="Run one or more H3 fragmentation-gradient primary cells.")
    value.add_argument("--cell-index", action="append", type=int, dest="cell_indices")
    value.add_argument("--master-seed", action="append", type=int, dest="master_seeds")
    value.add_argument("--patch-count", action="append", type=int, dest="patch_counts")
    value.add_argument("--replicates", type=int, default=20)
    value.add_argument("--generations", type=int, default=30)
    value.add_argument("--output-dir", type=Path, default=Path("artifacts/h3_fragmentation_gradient"))
    value.add_argument("--prefix", default="h3_fragmentation_gradient")
    return value


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    spec = replace(standard_profile(), replicates=args.replicates, generations=args.generations)
    seeds = tuple(FRESH_GRADIENT_MASTER_SEEDS if args.master_seeds is None else args.master_seeds)
    counts = tuple(FRAGMENT_PATCH_COUNTS if args.patch_counts is None else args.patch_counts)
    artifacts = run_h3_fragmentation_gradient_sensitivity(
        spec,
        primary_cell_indices=args.cell_indices,
        master_seeds=seeds,
        patch_counts=counts,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_fragmentation_gradient_artifacts(
        artifacts,
        csv_path=args.output_dir / f"{args.prefix}.csv",
        json_path=args.output_dir / f"{args.prefix}.json",
    )
    prepared = sum(
        record.source_prepared
        for artifact in artifacts
        for record in artifact.records
        if record.patch_count == counts[0]
    )
    attempted = sum(
        1
        for artifact in artifacts
        for record in artifact.records
        if record.patch_count == counts[0]
    )
    print(
        f"Wrote {len(artifacts)} cell artifact(s); "
        f"prepared sources {prepared}/{attempted}; patch counts={counts}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
