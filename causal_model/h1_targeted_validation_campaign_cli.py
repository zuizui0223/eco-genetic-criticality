"""Command-line runner for one-large canonical-H1-targeted validation campaigns."""
from __future__ import annotations

import argparse
import os
from dataclasses import replace
from pathlib import Path
from typing import Sequence

from causal_model.h1_targeted_validation_campaign import (
    DEFAULT_INSIDE_POSITIONS,
    run_h1_targeted_validation_campaign,
    targeted_output_paths,
    write_h1_targeted_validation_campaign,
)
from causal_model.multipatch_criticality_experiments import (
    ExperimentSpec,
    PROFILE_FULL,
    PROFILE_QUICK,
    PROFILE_STANDARD,
    full_profile,
    quick_profile,
    standard_profile,
)

PROFILE_FACTORIES = {
    PROFILE_QUICK: quick_profile,
    PROFILE_STANDARD: standard_profile,
    PROFILE_FULL: full_profile,
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run finite H1--H3 validation at barriers normalized to each one-large canonical H1 bistable interval."
        )
    )
    parser.add_argument("--profile", choices=tuple(PROFILE_FACTORIES), default=PROFILE_STANDARD)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/h1_targeted_validation"))
    parser.add_argument("--prefix", default=None, help="Output stem; defaults to h1_targeted_<profile>.")
    parser.add_argument("--replicates", type=int, default=None, help="Override profile replicate count.")
    parser.add_argument("--generations", type=int, default=None, help="Override profile generation count.")
    parser.add_argument("--master-seed", type=int, default=None, help="Override base master seed.")
    parser.add_argument(
        "--inside-position",
        action="append",
        type=float,
        dest="inside_positions",
        help="Fractional position inside each one-large bistable interval. Repeat; default is 0.25, 0.50, 0.75.",
    )
    parser.add_argument("--interaction-separation-threshold", type=float, default=0.05)
    parser.add_argument("--terminal-window", type=int, default=5)
    parser.add_argument("--hysteresis-barrier-points", type=int, default=7)
    parser.add_argument("--hysteresis-barrier-padding", type=float, default=0.1)
    parser.add_argument("--hysteresis-stage-generations", type=int, default=10)
    parser.add_argument("--hysteresis-low-state-threshold", type=float, default=0.25)
    parser.add_argument("--hysteresis-high-state-threshold", type=float, default=0.75)
    parser.add_argument("--coupled-tolerance", type=float, default=1e-12)
    parser.add_argument(
        "--code-revision",
        default=os.environ.get("GITHUB_SHA"),
        help="Optional source revision stored in the top-level manifest; defaults to GITHUB_SHA when available.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Allow replacement of all top-level and subcampaign artifacts.")
    return parser


def profile_spec(
    profile: str,
    *,
    replicates: int | None = None,
    generations: int | None = None,
    master_seed: int | None = None,
) -> ExperimentSpec:
    """Build a named profile with explicit validated run overrides."""
    if profile not in PROFILE_FACTORIES:
        raise ValueError(f"unknown profile: {profile}")
    updates = {
        name: value
        for name, value in (
            ("replicates", replicates),
            ("generations", generations),
            ("master_seed", master_seed),
        )
        if value is not None
    }
    return replace(PROFILE_FACTORIES[profile](), **updates)


def run_from_namespace(args: argparse.Namespace) -> dict[str, Path]:
    """Execute and write one H1-targeted campaign."""
    spec = profile_spec(
        args.profile,
        replicates=args.replicates,
        generations=args.generations,
        master_seed=args.master_seed,
    )
    positions = tuple(DEFAULT_INSIDE_POSITIONS if args.inside_positions is None else args.inside_positions)
    prefix = args.prefix or f"h1_targeted_{spec.profile}"
    paths = targeted_output_paths(args.output_dir, prefix)
    existing = tuple(path for path in paths.as_dict().values() if path.exists())
    if existing and not args.overwrite:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"refusing to overwrite existing artifacts: {names}; pass --overwrite to replace")
    audit_arguments = {
        "interaction_separation_threshold": args.interaction_separation_threshold,
        "terminal_window": args.terminal_window,
        "hysteresis_barrier_points": args.hysteresis_barrier_points,
        "hysteresis_barrier_padding": args.hysteresis_barrier_padding,
        "hysteresis_stage_generations": args.hysteresis_stage_generations,
        "hysteresis_low_state_threshold": args.hysteresis_low_state_threshold,
        "hysteresis_high_state_threshold": args.hysteresis_high_state_threshold,
        "coupled_tolerance": args.coupled_tolerance,
    }
    campaign = run_h1_targeted_validation_campaign(
        spec,
        normalized_positions=positions,
        **audit_arguments,
    )
    written = write_h1_targeted_validation_campaign(
        campaign,
        output_dir=args.output_dir,
        prefix=prefix,
        audit_arguments=audit_arguments,
        code_revision=args.code_revision,
    )
    return written.as_dict()


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        outputs = run_from_namespace(args)
    except (FileExistsError, ValueError) as error:
        build_parser().error(str(error))
    for name, path in sorted(outputs.items()):
        print(f"Wrote {name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
