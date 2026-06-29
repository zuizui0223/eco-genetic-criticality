"""Command-line runner for the unified finite H1--H3 validation campaign."""
from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict, replace
from pathlib import Path
from typing import Sequence

from causal_model.finite_validation_campaign import (
    campaign_output_paths,
    run_finite_validation_campaign,
    write_finite_validation_campaign_artifacts,
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
            "Run all finite H1--H3 validation audits on one shared parameter grid and write raw artifacts, a ledger, and a manifest."
        )
    )
    parser.add_argument("--profile", choices=tuple(PROFILE_FACTORIES), default=PROFILE_QUICK)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/finite_validation_campaign"))
    parser.add_argument("--prefix", default=None, help="Output stem; defaults to finite_validation_<profile>.")
    parser.add_argument("--replicates", type=int, default=None, help="Override profile replicate count.")
    parser.add_argument("--generations", type=int, default=None, help="Override profile generation count.")
    parser.add_argument("--master-seed", type=int, default=None, help="Override profile master random seed.")
    parser.add_argument(
        "--interaction-separation-threshold",
        type=float,
        default=0.05,
        help="Finite H1/H3 terminal interaction gap required for branch separation.",
    )
    parser.add_argument(
        "--terminal-window",
        type=int,
        default=5,
        help="Terminal snapshots averaged by the finite H1 branch audit and used for H2 conditioning.",
    )
    parser.add_argument("--hysteresis-barrier-points", type=int, default=7)
    parser.add_argument("--hysteresis-barrier-padding", type=float, default=0.1)
    parser.add_argument("--hysteresis-stage-generations", type=int, default=10)
    parser.add_argument("--hysteresis-low-state-threshold", type=float, default=0.25)
    parser.add_argument("--hysteresis-high-state-threshold", type=float, default=0.75)
    parser.add_argument("--coupled-tolerance", type=float, default=1e-12)
    parser.add_argument(
        "--code-revision",
        default=os.environ.get("GITHUB_SHA"),
        help="Optional source revision recorded in the manifest; defaults to GITHUB_SHA when available.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Allow replacement of all existing campaign outputs.")
    return parser


def profile_spec(
    profile: str,
    *,
    replicates: int | None = None,
    generations: int | None = None,
    master_seed: int | None = None,
) -> ExperimentSpec:
    """Build a named profile with explicit and validated run overrides."""
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


def manifest_path(output_dir: Path, prefix: str) -> Path:
    return output_dir / f"{prefix}.manifest.json"


def run_from_namespace(args: argparse.Namespace) -> dict[str, Path]:
    """Execute one campaign, write all artifacts, and return all output paths."""
    spec = profile_spec(
        args.profile,
        replicates=args.replicates,
        generations=args.generations,
        master_seed=args.master_seed,
    )
    prefix = args.prefix or f"finite_validation_{spec.profile}"
    paths = campaign_output_paths(args.output_dir, prefix)
    outputs = paths.as_dict()
    outputs["manifest"] = manifest_path(args.output_dir, prefix)
    existing = tuple(path for path in outputs.values() if path.exists())
    if existing and not args.overwrite:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"refusing to overwrite existing artifacts: {names}; pass --overwrite to replace")

    campaign = run_finite_validation_campaign(
        spec,
        interaction_separation_threshold=args.interaction_separation_threshold,
        terminal_window=args.terminal_window,
        hysteresis_barrier_points=args.hysteresis_barrier_points,
        hysteresis_barrier_padding=args.hysteresis_barrier_padding,
        hysteresis_stage_generations=args.hysteresis_stage_generations,
        hysteresis_low_state_threshold=args.hysteresis_low_state_threshold,
        hysteresis_high_state_threshold=args.hysteresis_high_state_threshold,
        coupled_tolerance=args.coupled_tolerance,
    )
    written = write_finite_validation_campaign_artifacts(
        campaign,
        output_dir=args.output_dir,
        prefix=prefix,
    ).as_dict()
    manifest = {
        "runner": "causal_model.finite_validation_campaign_cli",
        "campaign": "finite_validation_campaign_v1",
        "code_revision": args.code_revision,
        "profile": spec.profile,
        "spec": asdict(spec),
        "scenario_ids": list(campaign.scenario_ids),
        "audit_arguments": {
            "interaction_separation_threshold": args.interaction_separation_threshold,
            "terminal_window": args.terminal_window,
            "hysteresis_barrier_points": args.hysteresis_barrier_points,
            "hysteresis_barrier_padding": args.hysteresis_barrier_padding,
            "hysteresis_stage_generations": args.hysteresis_stage_generations,
            "hysteresis_low_state_threshold": args.hysteresis_low_state_threshold,
            "hysteresis_high_state_threshold": args.hysteresis_high_state_threshold,
            "coupled_tolerance": args.coupled_tolerance,
        },
        "cell_count": campaign.cell_count,
        "replicate_count_per_cell": spec.replicates,
        "selection_policy": (
            "All declared parameter cells and scenarios are retained. Canonical-H1-inapplicable cells, finite-precondition "
            "failures, censored first-passage comparisons, and counterexamples are reported rather than filtered."
        ),
        "outputs": {name: str(path) for name, path in written.items()},
    }
    outputs["manifest"].parent.mkdir(parents=True, exist_ok=True)
    with outputs["manifest"].open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
    written["manifest"] = outputs["manifest"]
    return written


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
