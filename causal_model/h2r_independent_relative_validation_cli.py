"""CLI for the frozen-domain independent H2-R relative-warning validation."""
from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path
from typing import Sequence

from causal_model.h2r_independent_relative_validation import (
    DEFAULT_H2R_VALIDATION_MASTER_SEEDS,
    run_h2r_independent_relative_validation,
    write_h2r_independent_validation_artifacts,
)
from causal_model.multipatch_criticality_experiments import (
    PROFILE_FULL,
    PROFILE_QUICK,
    PROFILE_STANDARD,
    full_profile,
    quick_profile,
    standard_profile,
)

_FACTORIES = {PROFILE_QUICK: quick_profile, PROFILE_STANDARD: standard_profile, PROFILE_FULL: full_profile}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run independent relative-warning validation in the frozen H2-R domain.")
    parser.add_argument("--profile", choices=tuple(_FACTORIES), default=PROFILE_STANDARD)
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--master-seed", action="append", type=int, dest="master_seeds")
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/h2r_independent_validation"))
    parser.add_argument("--prefix", default="independent_relative_warning_v1")
    parser.add_argument("--endpoint-padding-fraction", type=float, default=0.5)
    parser.add_argument("--stage-generations", type=int, default=30)
    parser.add_argument("--hold-generations", type=int, default=30)
    parser.add_argument("--barrier-points", action="append", type=int, dest="nested_barrier_points")
    parser.add_argument("--interaction-separation-threshold", type=float, default=0.05)
    parser.add_argument("--maximum-normalized-bracket-width", type=float, default=0.03)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    spec = replace(_FACTORIES[args.profile](), replicates=args.replicates)
    seeds = tuple(DEFAULT_H2R_VALIDATION_MASTER_SEEDS if args.master_seeds is None else args.master_seeds)
    points = tuple((25, 49, 97) if args.nested_barrier_points is None else args.nested_barrier_points)
    result = run_h2r_independent_relative_validation(
        spec,
        master_seeds=seeds,
        endpoint_padding_fraction=args.endpoint_padding_fraction,
        stage_generations=args.stage_generations,
        hold_generations=args.hold_generations,
        nested_barrier_points=points,
        interaction_separation_threshold=args.interaction_separation_threshold,
        maximum_normalized_bracket_width=args.maximum_normalized_bracket_width,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_h2r_independent_validation_artifacts(
        result,
        csv_path=args.output_dir / f"{args.prefix}.csv",
        json_path=args.output_dir / f"{args.prefix}.json",
        manifest_path=args.output_dir / f"{args.prefix}.manifest.json",
    )
    print(f"Wrote H2-R independent validation artifacts to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
