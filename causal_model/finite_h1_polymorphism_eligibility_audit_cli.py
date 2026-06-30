"""Command-line runner for H1 high-branch genetic eligibility audit."""
from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict, replace
from pathlib import Path
from typing import Sequence

from causal_model.finite_h1_polymorphism_eligibility_audit import (
    run_finite_h1_polymorphism_eligibility_audit,
    write_finite_h1_polymorphism_eligibility_artifacts,
)
from causal_model.finite_h1_fragment_projection_audit import DEFAULT_MASTER_SEEDS
from causal_model.multipatch_criticality_experiments import (
    PROFILE_FULL,
    PROFILE_QUICK,
    PROFILE_STANDARD,
    ExperimentSpec,
    full_profile,
    quick_profile,
    standard_profile,
)

_FACTORIES = {PROFILE_QUICK: quick_profile, PROFILE_STANDARD: standard_profile, PROFILE_FULL: full_profile}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Check whether full-state H1 high branches retain genetic variation for H2/H3 endpoints.")
    parser.add_argument("--profile", choices=tuple(_FACTORIES), default=PROFILE_STANDARD)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/h1_polymorphism_eligibility"))
    parser.add_argument("--prefix", default=None)
    parser.add_argument("--replicates", type=int, default=None)
    parser.add_argument("--generations", type=int, default=None)
    parser.add_argument("--master-seed", action="append", type=int, dest="master_seeds")
    parser.add_argument("--endpoint-padding-fraction", type=float, default=0.5)
    parser.add_argument("--stage-generations", type=int, default=30)
    parser.add_argument("--hold-generations", type=int, default=30)
    parser.add_argument("--barrier-points", action="append", type=int, dest="nested_barrier_points")
    parser.add_argument("--interaction-separation-threshold", type=float, default=0.05)
    parser.add_argument("--maximum-normalized-bracket-width", type=float, default=0.03)
    parser.add_argument("--polymorphism-epsilon", type=float, default=1e-12)
    parser.add_argument("--code-revision", default=os.environ.get("GITHUB_SHA"))
    parser.add_argument("--overwrite", action="store_true")
    return parser


def _profile_spec(profile: str, replicates: int | None, generations: int | None) -> ExperimentSpec:
    updates = {name: value for name, value in (("replicates", replicates), ("generations", generations)) if value is not None}
    return replace(_FACTORIES[profile](), **updates)


def run_from_namespace(args: argparse.Namespace) -> dict[str, Path]:
    spec = _profile_spec(args.profile, args.replicates, args.generations)
    seeds = tuple(DEFAULT_MASTER_SEEDS if args.master_seeds is None else args.master_seeds)
    points = tuple((25, 49, 97) if args.nested_barrier_points is None else args.nested_barrier_points)
    seed_tag = "_".join(str(seed) for seed in seeds)
    prefix = args.prefix or f"h1_polymorphism_eligibility_{spec.profile}_seeds_{seed_tag}"
    csv_path = args.output_dir / f"{prefix}.csv"
    json_path = args.output_dir / f"{prefix}.json"
    manifest_path = args.output_dir / f"{prefix}.manifest.json"
    existing = tuple(path for path in (csv_path, json_path, manifest_path) if path.exists())
    if existing and not args.overwrite:
        raise FileExistsError("refusing to overwrite existing artifacts: " + ", ".join(map(str, existing)))
    cells = run_finite_h1_polymorphism_eligibility_audit(
        spec,
        master_seeds=seeds,
        endpoint_padding_fraction=args.endpoint_padding_fraction,
        stage_generations=args.stage_generations,
        hold_generations=args.hold_generations,
        nested_barrier_points=points,
        interaction_separation_threshold=args.interaction_separation_threshold,
        maximum_normalized_bracket_width=args.maximum_normalized_bracket_width,
        polymorphism_epsilon=args.polymorphism_epsilon,
    )
    write_finite_h1_polymorphism_eligibility_artifacts(cells, csv_path=csv_path, json_path=json_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = {
        "runner": "causal_model.finite_h1_polymorphism_eligibility_audit_cli",
        "campaign": "finite_h1_polymorphism_eligibility_v1",
        "code_revision": args.code_revision,
        "spec": asdict(spec),
        "master_seeds": list(seeds),
        "h1_source_preparation": {
            "endpoint_padding_fraction": args.endpoint_padding_fraction,
            "stage_generations": args.stage_generations,
            "hold_generations": args.hold_generations,
            "nested_barrier_points": list(points),
            "interaction_separation_threshold": args.interaction_separation_threshold,
            "maximum_normalized_bracket_width": args.maximum_normalized_bracket_width,
        },
        "genetic_eligibility_rule": {
            "polymorphism_epsilon": args.polymorphism_epsilon,
            "h2_primary_landscape": "equal_isolated",
            "h2_dynamic_warning_eligible": "epsilon < p < 1-epsilon; H-alpha > threshold; H-gamma > threshold at generation zero",
            "h3_genetic_contrast_eligible": "epsilon < p < 1-epsilon; H-alpha > 0; H-gamma > 0 at generation zero",
            "preexisting_warning_policy": "A warning already present at generation zero is baseline-ineligible for H2 lead-time inference, not a warning lead.",
        },
        "selection_policy": "Every declared seed, pair, and replicate is retained. An unavailable H1 source or a baseline-ineligible genetic state remains explicit in the artifact.",
        "cell_count": len(cells),
        "outputs": {"csv": str(csv_path), "json": str(json_path), "manifest": str(manifest_path)},
    }
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
    return {"csv": csv_path, "json": json_path, "manifest": manifest_path}


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
