"""Command-line runner for paired trait/genetic/full baseline comparisons."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict, replace
from pathlib import Path
from typing import Sequence

from causal_model.multipatch_criticality_experiments import (
    ExperimentSpec,
    LandscapeScenario,
    PROFILE_FULL,
    PROFILE_QUICK,
    PROFILE_STANDARD,
    SCENARIO_EQUAL_ISOLATED,
    SCENARIO_EQUAL_MIGRATING,
    SCENARIO_ONE_LARGE,
    full_profile,
    scenario_equal_isolated,
    scenario_equal_migrating,
    scenario_one_large,
    standard_profile,
)
from causal_model.paired_baseline_comparisons import (
    BASELINE_IDS,
    baseline_definition,
    comparison_quick_profile,
    run_paired_baseline_comparisons,
    write_paired_baseline_artifacts,
)

PROFILE_FACTORIES = {
    PROFILE_QUICK: comparison_quick_profile,
    PROFILE_STANDARD: standard_profile,
    PROFILE_FULL: full_profile,
}
SCENARIO_CHOICES = (
    SCENARIO_ONE_LARGE,
    SCENARIO_EQUAL_ISOLATED,
    SCENARIO_EQUAL_MIGRATING,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run matched trait-only, genetic-only, and full eco-genetic ablations and write CSV, JSON, and a manifest."
        )
    )
    parser.add_argument("--profile", choices=tuple(PROFILE_FACTORIES), default=PROFILE_QUICK)
    parser.add_argument(
        "--scenario",
        choices=SCENARIO_CHOICES,
        action="append",
        dest="scenarios",
        help="Landscape scenario to run. Repeat to select several; default runs all three.",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/baseline_comparison"))
    parser.add_argument("--prefix", default=None, help="Output stem; defaults to baseline_comparison_<profile>.")
    parser.add_argument("--replicates", type=int, default=None, help="Override profile replicate count.")
    parser.add_argument("--generations", type=int, default=None, help="Override profile generation count.")
    parser.add_argument("--master-seed", type=int, default=None, help="Override profile master random seed.")
    parser.add_argument("--overwrite", action="store_true", help="Allow replacement of existing output files.")
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
    updates: dict[str, int] = {}
    for name, value in (("replicates", replicates), ("generations", generations), ("master_seed", master_seed)):
        if value is not None:
            updates[name] = value
    return replace(PROFILE_FACTORIES[profile](), **updates)


def select_scenarios(spec: ExperimentSpec, names: Sequence[str] | None) -> tuple[LandscapeScenario, ...]:
    factories = {
        SCENARIO_ONE_LARGE: scenario_one_large,
        SCENARIO_EQUAL_ISOLATED: scenario_equal_isolated,
        SCENARIO_EQUAL_MIGRATING: scenario_equal_migrating,
    }
    selected = SCENARIO_CHOICES if not names else tuple(dict.fromkeys(names))
    return tuple(factories[name](spec) for name in selected)


def output_paths(output_dir: Path, prefix: str) -> tuple[Path, Path, Path]:
    return (
        output_dir / f"{prefix}.csv",
        output_dir / f"{prefix}.json",
        output_dir / f"{prefix}.manifest.json",
    )


def run_from_namespace(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    spec = profile_spec(
        args.profile,
        replicates=args.replicates,
        generations=args.generations,
        master_seed=args.master_seed,
    )
    scenarios = select_scenarios(spec, args.scenarios)
    prefix = args.prefix or f"baseline_comparison_{spec.profile}"
    csv_path, json_path, manifest_path = output_paths(args.output_dir, prefix)
    destinations = (csv_path, json_path, manifest_path)
    existing = tuple(path for path in destinations if path.exists())
    if existing and not args.overwrite:
        raise FileExistsError(
            f"refusing to overwrite existing artifacts: {', '.join(str(path) for path in existing)}; pass --overwrite to replace"
        )

    cells = run_paired_baseline_comparisons(spec, scenarios=scenarios)
    write_paired_baseline_artifacts(cells, csv_path=csv_path, json_path=json_path)
    baseline_parameters = spec.base_parameters
    manifest = {
        "runner": "causal_model.paired_baseline_cli",
        "profile": spec.profile,
        "spec": asdict(spec),
        "scenario_ids": [scenario.scenario_id for scenario in scenarios],
        "baseline_ids": list(BASELINE_IDS),
        "baseline_definitions": [
            baseline_definition(baseline_parameters, baseline_id).as_dict() for baseline_id in BASELINE_IDS
        ],
        "shared_seed_rule": "Each baseline in the same cell and replicate uses the same derived seed.",
        "outputs": {"csv": str(csv_path), "json": str(json_path)},
        "cell_count": len(cells),
        "replicate_count_per_cell": spec.replicates,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
    return csv_path, json_path, manifest_path


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        csv_path, json_path, manifest_path = run_from_namespace(args)
    except (FileExistsError, ValueError) as error:
        build_parser().error(str(error))
    print(f"Wrote CSV: {csv_path}")
    print(f"Wrote JSON: {json_path}")
    print(f"Wrote manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
