"""Command-line runner for multi-master-seed finite H1 boundary robustness."""
from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict, replace
from pathlib import Path
from typing import Sequence

from causal_model.finite_h1_boundary_seed_ensemble import (
    DEFAULT_ENSEMBLE_MASTER_SEEDS,
    run_finite_h1_boundary_seed_ensemble,
    write_finite_h1_boundary_seed_ensemble_artifacts,
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
            "Repeat the finite H1 nested-grid boundary audit across a declared ensemble of independent master seeds."
        )
    )
    parser.add_argument("--profile", choices=tuple(PROFILE_FACTORIES), default=PROFILE_STANDARD)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/h1_boundary_seed_ensemble"))
    parser.add_argument("--prefix", default=None)
    parser.add_argument("--replicates", type=int, default=None)
    parser.add_argument("--generations", type=int, default=None)
    parser.add_argument(
        "--ensemble-master-seed",
        action="append",
        type=int,
        dest="ensemble_master_seeds",
        help="One independent master seed. Repeat; default 20260630,20260631,20260632,20260633,20260634.",
    )
    parser.add_argument("--endpoint-padding-fraction", type=float, default=0.5)
    parser.add_argument("--stage-generations", type=int, default=30)
    parser.add_argument(
        "--barrier-points",
        action="append",
        type=int,
        dest="nested_barrier_points",
        help="Nested barrier-grid point count. Repeat; default 25,49,97.",
    )
    parser.add_argument("--interaction-separation-threshold", type=float, default=0.05)
    parser.add_argument("--low-state-threshold", type=float, default=0.25)
    parser.add_argument("--high-state-threshold", type=float, default=0.75)
    parser.add_argument("--maximum-normalized-bracket-width", type=float, default=0.03)
    parser.add_argument("--code-revision", default=os.environ.get("GITHUB_SHA"))
    parser.add_argument("--overwrite", action="store_true")
    return parser


def profile_spec(
    profile: str,
    *,
    replicates: int | None = None,
    generations: int | None = None,
) -> ExperimentSpec:
    if profile not in PROFILE_FACTORIES:
        raise ValueError(f"unknown profile: {profile}")
    updates = {
        name: value
        for name, value in (("replicates", replicates), ("generations", generations))
        if value is not None
    }
    return replace(PROFILE_FACTORIES[profile](), **updates)


def run_from_namespace(args: argparse.Namespace) -> dict[str, Path]:
    spec = profile_spec(args.profile, replicates=args.replicates, generations=args.generations)
    seeds = tuple(DEFAULT_ENSEMBLE_MASTER_SEEDS if args.ensemble_master_seeds is None else args.ensemble_master_seeds)
    grids = tuple((25, 49, 97) if args.nested_barrier_points is None else args.nested_barrier_points)
    seed_tag = "_".join(str(value) for value in seeds)
    prefix = args.prefix or f"h1_boundary_seed_ensemble_{spec.profile}_seeds_{seed_tag}"
    csv_path = args.output_dir / f"{prefix}.csv"
    json_path = args.output_dir / f"{prefix}.json"
    manifest_path = args.output_dir / f"{prefix}.manifest.json"
    existing = tuple(path for path in (csv_path, json_path, manifest_path) if path.exists())
    if existing and not args.overwrite:
        raise FileExistsError(
            "refusing to overwrite existing artifacts: " + ", ".join(str(path) for path in existing)
        )
    cells = run_finite_h1_boundary_seed_ensemble(
        spec,
        master_seeds=seeds,
        endpoint_padding_fraction=args.endpoint_padding_fraction,
        stage_generations=args.stage_generations,
        nested_barrier_points=grids,
        interaction_separation_threshold=args.interaction_separation_threshold,
        low_state_threshold=args.low_state_threshold,
        high_state_threshold=args.high_state_threshold,
        maximum_normalized_bracket_width=args.maximum_normalized_bracket_width,
    )
    write_finite_h1_boundary_seed_ensemble_artifacts(cells, csv_path=csv_path, json_path=json_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = {
        "runner": "causal_model.finite_h1_boundary_seed_ensemble_cli",
        "campaign": "finite_h1_boundary_seed_ensemble_v1",
        "code_revision": args.code_revision,
        "spec": asdict(spec),
        "ensemble_master_seeds": list(seeds),
        "endpoint_padding_fraction": args.endpoint_padding_fraction,
        "nested_barrier_points": list(grids),
        "audit_arguments": {
            "stage_generations": args.stage_generations,
            "interaction_separation_threshold": args.interaction_separation_threshold,
            "low_state_threshold": args.low_state_threshold,
            "high_state_threshold": args.high_state_threshold,
            "maximum_normalized_bracket_width": args.maximum_normalized_bracket_width,
        },
        "cell_count": len(cells),
        "selection_policy": (
            "Every declared (area_reference, interaction_feedback) pair is evaluated under every declared master seed on the same "
            "endpoint range and nested grids. No seed-specific boundary, loop, trait-switch, or support outcome controls inclusion."
        ),
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
