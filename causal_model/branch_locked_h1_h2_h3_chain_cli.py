"""Command-line runner for the branch-locked same-replicate H1--H2--H3 audit."""
from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict, replace
from pathlib import Path
from typing import Sequence

from causal_model.branch_locked_h1_h2_h3_chain import (
    DEFAULT_CHAIN_MASTER_SEEDS,
    run_branch_locked_h1_h2_h3_chain_audit,
    write_branch_locked_h1_h2_h3_chain_artifacts,
)
from causal_model.multipatch_criticality_experiments import (
    PROFILE_FULL,
    PROFILE_QUICK,
    PROFILE_STANDARD,
    ExperimentSpec,
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
            "Calibrate finite H1 anchors, then evaluate same-replicate H2 warning order and high-branch H3 contrasts."
        )
    )
    parser.add_argument("--profile", choices=tuple(PROFILE_FACTORIES), default=PROFILE_STANDARD)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/branch_locked_h1_h2_h3"))
    parser.add_argument("--prefix", default=None)
    parser.add_argument("--replicates", type=int, default=None)
    parser.add_argument("--generations", type=int, default=None)
    parser.add_argument(
        "--master-seed",
        action="append",
        type=int,
        dest="master_seeds",
        help="One independent master seed. Repeat; default 20260630..20260634.",
    )
    parser.add_argument("--h1-endpoint-padding-fraction", type=float, default=0.5)
    parser.add_argument("--h1-stage-generations", type=int, default=30)
    parser.add_argument(
        "--h1-barrier-points",
        action="append",
        type=int,
        dest="h1_nested_barrier_points",
        help="Nested H1 calibration grid. Repeat; default 25,49,97.",
    )
    parser.add_argument("--interaction-separation-threshold", type=float, default=0.05)
    parser.add_argument("--h1-maximum-normalized-bracket-width", type=float, default=0.03)
    parser.add_argument("--code-revision", default=os.environ.get("GITHUB_SHA"))
    parser.add_argument("--overwrite", action="store_true")
    return parser


def profile_spec(profile: str, *, replicates: int | None, generations: int | None) -> ExperimentSpec:
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
    seeds = tuple(DEFAULT_CHAIN_MASTER_SEEDS if args.master_seeds is None else args.master_seeds)
    grids = tuple((25, 49, 97) if args.h1_nested_barrier_points is None else args.h1_nested_barrier_points)
    seed_tag = "_".join(str(seed) for seed in seeds)
    prefix = args.prefix or f"branch_locked_h1_h2_h3_{spec.profile}_seeds_{seed_tag}"
    csv_path = args.output_dir / f"{prefix}.csv"
    json_path = args.output_dir / f"{prefix}.json"
    manifest_path = args.output_dir / f"{prefix}.manifest.json"
    existing = tuple(path for path in (csv_path, json_path, manifest_path) if path.exists())
    if existing and not args.overwrite:
        raise FileExistsError("refusing to overwrite existing artifacts: " + ", ".join(str(path) for path in existing))
    cells = run_branch_locked_h1_h2_h3_chain_audit(
        spec,
        master_seeds=seeds,
        h1_endpoint_padding_fraction=args.h1_endpoint_padding_fraction,
        h1_stage_generations=args.h1_stage_generations,
        h1_nested_barrier_points=grids,
        interaction_separation_threshold=args.interaction_separation_threshold,
        h1_maximum_normalized_bracket_width=args.h1_maximum_normalized_bracket_width,
    )
    write_branch_locked_h1_h2_h3_chain_artifacts(cells, csv_path=csv_path, json_path=json_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = {
        "runner": "causal_model.branch_locked_h1_h2_h3_chain_cli",
        "campaign": "branch_locked_h1_h2_h3_chain_v1",
        "code_revision": args.code_revision,
        "spec": asdict(spec),
        "master_seeds": list(seeds),
        "h1_calibration": {
            "endpoint_padding_fraction": args.h1_endpoint_padding_fraction,
            "stage_generations": args.h1_stage_generations,
            "nested_barrier_points": list(grids),
            "interaction_separation_threshold": args.interaction_separation_threshold,
            "maximum_normalized_bracket_width": args.h1_maximum_normalized_bracket_width,
        },
        "h2_primary_endpoint": {
            "landscape": "equal_isolated",
            "branch": "high_start",
            "warnings": ["H_alpha", "H_gamma"],
            "censoring_policy": "A warning order is valid only when both warning and realised trait-loss events occur; absent events remain censored.",
        },
        "h3_primary_endpoint": {
            "landscape_comparison": "equal_isolated_minus_one_large",
            "branch": "high_start",
            "required_signs": {
                "interaction": "isolated < one_large",
                "local_effective_size": "isolated < one_large",
                "realised_high_trait_mass": "isolated < one_large",
            },
            "migration_interpretation": "allele-frequency mixing only; no demographic rescue or recolonisation is modelled",
        },
        "cell_count": len(cells),
        "selection_policy": (
            "Every declared seed, parameter pair, and replicate first receives the same H1 calibration. "
            "A missing H1 anchor or censored H2 event remains in the artifact and denominator ledger; no H2/H3 outcome selects records for inclusion."
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
