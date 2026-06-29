"""Command-line runner for the finite H1 stage-duration robustness audit."""
from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict, replace
from pathlib import Path
from typing import Sequence

from causal_model.finite_h1_hysteresis_duration_audit import (
    run_finite_h1_hysteresis_duration_audit,
    write_finite_h1_hysteresis_duration_artifacts,
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
            "Repeat one-large finite H1 hysteresis continuation across stage durations to distinguish persistent route memory from short-stage effects."
        )
    )
    parser.add_argument("--profile", choices=tuple(PROFILE_FACTORIES), default=PROFILE_STANDARD)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/h1_hysteresis_duration"))
    parser.add_argument("--prefix", default=None)
    parser.add_argument("--replicates", type=int, default=None)
    parser.add_argument("--generations", type=int, default=None)
    parser.add_argument("--master-seed", type=int, default=None)
    parser.add_argument(
        "--stage-generations",
        action="append",
        type=int,
        dest="stage_duration_values",
        help="One continuation duration. Repeat in strictly increasing order; default 5,10,30,80.",
    )
    parser.add_argument("--barrier-points", type=int, default=7)
    parser.add_argument("--barrier-padding", type=float, default=0.1)
    parser.add_argument("--interaction-separation-threshold", type=float, default=0.05)
    parser.add_argument("--low-state-threshold", type=float, default=0.25)
    parser.add_argument("--high-state-threshold", type=float, default=0.75)
    parser.add_argument("--gap-stability-tolerance", type=float, default=0.05)
    parser.add_argument("--code-revision", default=os.environ.get("GITHUB_SHA"))
    parser.add_argument("--overwrite", action="store_true")
    return parser


def profile_spec(
    profile: str,
    *,
    replicates: int | None = None,
    generations: int | None = None,
    master_seed: int | None = None,
) -> ExperimentSpec:
    if profile not in PROFILE_FACTORIES:
        raise ValueError(f"unknown profile: {profile}")
    updates = {
        key: value
        for key, value in (
            ("replicates", replicates),
            ("generations", generations),
            ("master_seed", master_seed),
        )
        if value is not None
    }
    return replace(PROFILE_FACTORIES[profile](), **updates)


def run_from_namespace(args: argparse.Namespace) -> dict[str, Path]:
    spec = profile_spec(
        args.profile,
        replicates=args.replicates,
        generations=args.generations,
        master_seed=args.master_seed,
    )
    durations = tuple((5, 10, 30, 80) if args.stage_duration_values is None else args.stage_duration_values)
    prefix = args.prefix or f"h1_hysteresis_duration_{spec.profile}"
    csv_path = args.output_dir / f"{prefix}.csv"
    json_path = args.output_dir / f"{prefix}.json"
    manifest_path = args.output_dir / f"{prefix}.manifest.json"
    existing = tuple(path for path in (csv_path, json_path, manifest_path) if path.exists())
    if existing and not args.overwrite:
        raise FileExistsError(
            "refusing to overwrite existing artifacts: " + ", ".join(str(path) for path in existing)
        )
    cells = run_finite_h1_hysteresis_duration_audit(
        spec,
        stage_generations=durations,
        barrier_points=args.barrier_points,
        barrier_padding=args.barrier_padding,
        interaction_separation_threshold=args.interaction_separation_threshold,
        low_state_threshold=args.low_state_threshold,
        high_state_threshold=args.high_state_threshold,
        gap_stability_tolerance=args.gap_stability_tolerance,
    )
    write_finite_h1_hysteresis_duration_artifacts(cells, csv_path=csv_path, json_path=json_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = {
        "runner": "causal_model.finite_h1_hysteresis_duration_audit_cli",
        "campaign": "finite_h1_hysteresis_duration_v1",
        "code_revision": args.code_revision,
        "spec": asdict(spec),
        "stage_generations": list(durations),
        "audit_arguments": {
            "barrier_points": args.barrier_points,
            "barrier_padding": args.barrier_padding,
            "interaction_separation_threshold": args.interaction_separation_threshold,
            "low_state_threshold": args.low_state_threshold,
            "high_state_threshold": args.high_state_threshold,
            "gap_stability_tolerance": args.gap_stability_tolerance,
        },
        "cell_count": len(cells),
        "selection_policy": (
            "Every declared (area_reference, interaction_feedback) pair is retained. A pair without one-large strict canonical "
            "bistability is recorded as unavailable; no finite support, gap, or trait outcome controls inclusion."
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
