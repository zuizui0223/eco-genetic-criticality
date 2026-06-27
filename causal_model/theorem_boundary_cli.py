"""Command-line runner for reproducible theorem-boundary phase diagrams."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict, replace
from pathlib import Path
from typing import Sequence

from causal_model.multipatch_criticality_experiments import (
    ExperimentSpec,
    PROFILE_FULL,
    PROFILE_QUICK,
    PROFILE_STANDARD,
    SCENARIO_EQUAL_ISOLATED,
    SCENARIO_EQUAL_MIGRATING,
    SCENARIO_ONE_LARGE,
    LandscapeScenario,
    full_profile,
    quick_profile,
    scenario_equal_isolated,
    scenario_equal_migrating,
    scenario_one_large,
    standard_profile,
)
from causal_model.theorem_boundary_phase_diagram import (
    run_theorem_boundary_phase_diagram,
    write_theorem_boundary_phase_artifacts,
)

PROFILE_FACTORIES = {
    PROFILE_QUICK: quick_profile,
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
            "Run a censoring-aware theorem-boundary phase diagram and write CSV, JSON, and a reproducibility manifest."
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
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/theorem_boundary"))
    parser.add_argument("--prefix", default=None, help="Output stem; defaults to theorem_boundary_<profile>.")
    parser.add_argument("--replicates", type=int, default=None, help="Override profile replicate count.")
    parser.add_argument("--generations", type=int, default=None, help="Override profile generation count.")
    parser.add_argument("--master-seed", type=int, default=None, help="Override profile master random seed.")
    parser.add_argument("--tolerance", type=float, default=1e-12, help="Canonical-update equality tolerance.")
    parser.add_argument("--overwrite", action="store_true", help="Allow replacement of existing output files.")
    return parser


def profile_spec(
    profile: str,
    *,
    replicates: int | None = None,
    generations: int | None = None,
    master_seed: int | None = None,
) -> ExperimentSpec:
    """Build a named profile with explicit, validated command-line overrides."""
    if profile not in PROFILE_FACTORIES:
        raise ValueError(f"unknown profile: {profile}")
    updates: dict[str, int] = {}
    for name, value in (
        ("replicates", replicates),
        ("generations", generations),
        ("master_seed", master_seed),
    ):
        if value is not None:
            updates[name] = value
    return replace(PROFILE_FACTORIES[profile](), **updates)


def select_scenarios(spec: ExperimentSpec, names: Sequence[str] | None) -> tuple[LandscapeScenario, ...]:
    """Resolve CLI scenario names in a stable, caller-visible order."""
    factories = {
        SCENARIO_ONE_LARGE: scenario_one_large,
        SCENARIO_EQUAL_ISOLATED: scenario_equal_isolated,
        SCENARIO_EQUAL_MIGRATING: scenario_equal_migrating,
    }
    selected = SCENARIO_CHOICES if not names else tuple(dict.fromkeys(names))
    return tuple(factories[name](spec) for name in selected)


def output_paths(output_dir: Path, prefix: str) -> tuple[Path, Path, Path]:
    """Return CSV, full JSON, and manifest destinations for a run."""
    return (
        output_dir / f"{prefix}.csv",
        output_dir / f"{prefix}.json",
        output_dir / f"{prefix}.manifest.json",
    )


def run_from_namespace(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    """Execute a parsed CLI request and return the written artifact paths."""
    spec = profile_spec(
        args.profile,
        replicates=args.replicates,
        generations=args.generations,
        master_seed=args.master_seed,
    )
    scenarios = select_scenarios(spec, args.scenarios)
    prefix = args.prefix or f"theorem_boundary_{spec.profile}"
    csv_path, json_path, manifest_path = output_paths(args.output_dir, prefix)
    destinations = (csv_path, json_path, manifest_path)
    existing = tuple(path for path in destinations if path.exists())
    if existing and not args.overwrite:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"refusing to overwrite existing artifacts: {names}; pass --overwrite to replace")

    cells = run_theorem_boundary_phase_diagram(spec, scenarios=scenarios, tolerance=args.tolerance)
    write_theorem_boundary_phase_artifacts(cells, csv_path=csv_path, json_path=json_path)
    manifest = {
        "runner": "causal_model.theorem_boundary_cli",
        "profile": spec.profile,
        "spec": asdict(spec),
        "scenario_ids": [scenario.scenario_id for scenario in scenarios],
        "tolerance": args.tolerance,
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
