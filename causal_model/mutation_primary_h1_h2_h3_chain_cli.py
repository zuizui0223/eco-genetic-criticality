"""Command-line runner for mutation-conditioned primary H1-H2-H3 chain."""
from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path
from typing import Sequence

from causal_model.multipatch_criticality_experiments import PROFILE_FULL, PROFILE_QUICK, PROFILE_STANDARD, full_profile, quick_profile, standard_profile
from causal_model.mutation_primary_h1_h2_h3_chain import (
    DEFAULT_PRIMARY_CHAIN_MASTER_SEEDS,
    run_mutation_primary_h1_h2_h3_chain,
    write_mutation_primary_h1_h2_h3_artifacts,
)

_FACTORIES = {PROFILE_QUICK: quick_profile, PROFILE_STANDARD: standard_profile, PROFILE_FULL: full_profile}


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description="Run full-state mutation-conditioned H1-H2-H3 primary chain.")
    value.add_argument("--profile", choices=tuple(_FACTORIES), default=PROFILE_STANDARD)
    value.add_argument("--replicates", type=int, default=None)
    value.add_argument("--generations", type=int, default=None)
    value.add_argument("--master-seed", action="append", type=int, dest="master_seeds")
    value.add_argument("--output-dir", type=Path, default=Path("artifacts/mutation_primary_h1_h2_h3"))
    value.add_argument("--prefix", default="primary_chain_v1")
    value.add_argument("--endpoint-padding-fraction", type=float, default=0.5)
    value.add_argument("--stage-generations", type=int, default=30)
    value.add_argument("--hold-generations", type=int, default=30)
    value.add_argument("--barrier-points", action="append", type=int, dest="nested_barrier_points")
    value.add_argument("--interaction-separation-threshold", type=float, default=0.05)
    value.add_argument("--maximum-normalized-bracket-width", type=float, default=0.03)
    return value


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    changes = {key: value for key, value in (("replicates", args.replicates), ("generations", args.generations)) if value is not None}
    spec = replace(_FACTORIES[args.profile](), **changes)
    seeds = tuple(DEFAULT_PRIMARY_CHAIN_MASTER_SEEDS if args.master_seeds is None else args.master_seeds)
    points = tuple((25, 49, 97) if args.nested_barrier_points is None else args.nested_barrier_points)
    cells = run_mutation_primary_h1_h2_h3_chain(
        spec, master_seeds=seeds, endpoint_padding_fraction=args.endpoint_padding_fraction,
        stage_generations=args.stage_generations, hold_generations=args.hold_generations,
        nested_barrier_points=points, interaction_separation_threshold=args.interaction_separation_threshold,
        maximum_normalized_bracket_width=args.maximum_normalized_bracket_width,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_mutation_primary_h1_h2_h3_artifacts(
        cells, csv_path=args.output_dir / f"{args.prefix}.csv", json_path=args.output_dir / f"{args.prefix}.json",
        manifest_path=args.output_dir / f"{args.prefix}.manifest.json",
    )
    print(f"Wrote mutation primary chain artifacts to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
