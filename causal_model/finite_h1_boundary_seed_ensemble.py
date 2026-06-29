"""Multi-master-seed robustness audit for finite H1 boundary brackets.

A single independent master seed can show that a nested-grid collapse/recovery
bracket is not a replay of an earlier stochastic stream.  It is still only one
master-seed realization.  This module repeats the *same predeclared* boundary-
resolution design across a fixed ensemble of independent master seeds and keeps
the seed-specific raw continuations alongside pooled and seed-level summaries.

The audit is deliberately conservative:

* every declared master seed is evaluated for every declared H1 parameter pair;
* a pair is never dropped because it lacks a loop or trait-switch response;
* seed-level and pooled replicate probabilities are reported separately; and
* ``all_master_seed_runs_fully_supported`` requires 100% support within every
  declared master-seed run, not merely a high pooled proportion.

This is Type S evidence for the finite closure.  Multiple master seeds make a
numerical robustness claim stronger, but do not turn a finite simulation into a
bifurcation theorem.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from statistics import median
from typing import Iterable, Mapping, Sequence

from causal_model.finite_h1_boundary_resolution_audit import (
    DEFAULT_NESTED_BARRIER_POINTS,
    FiniteH1BoundaryResolutionCell,
    run_finite_h1_boundary_resolution_audit,
)
from causal_model.multipatch_criticality_experiments import ExperimentSpec, ParameterCell

DEFAULT_ENSEMBLE_MASTER_SEEDS = (20260630, 20260631, 20260632, 20260633, 20260634)


@dataclass(frozen=True)
class FiniteH1BoundarySeedRun:
    """One master-seed realization of one finite H1 parameter cell."""

    master_seed: int
    cell: FiniteH1BoundaryResolutionCell

    def as_dict(self) -> dict[str, object]:
        return {"master_seed": self.master_seed, "cell": self.cell.as_dict()}


@dataclass(frozen=True)
class FiniteH1BoundarySeedEnsembleCell:
    """Aligned multi-master-seed evidence for one (A_ref, kappa) pair."""

    experiment_id: str
    profile: str
    pair_index: int
    parameters: ParameterCell
    canonical_bistable_barrier_interval: tuple[float, float] | None
    endpoint_padding_fraction: float
    stage_generations: int
    nested_barrier_points: tuple[int, ...]
    maximum_normalized_bracket_width: float
    seed_runs: tuple[FiniteH1BoundarySeedRun, ...]
    summary: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "pair_index": self.pair_index,
            "parameters": asdict(self.parameters),
            "canonical_bistable_barrier_interval": None
            if self.canonical_bistable_barrier_interval is None
            else list(self.canonical_bistable_barrier_interval),
            "endpoint_padding_fraction": self.endpoint_padding_fraction,
            "stage_generations": self.stage_generations,
            "nested_barrier_points": list(self.nested_barrier_points),
            "maximum_normalized_bracket_width": self.maximum_normalized_bracket_width,
            "master_seed_count": len(self.seed_runs),
            "seed_runs": [run.as_dict() for run in self.seed_runs],
            "summary": dict(self.summary),
        }

    def to_csv_row(self) -> dict[str, object]:
        row: dict[str, object] = {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "pair_index": self.pair_index,
            "master_seed_count": len(self.seed_runs),
            "master_seeds": ",".join(str(run.master_seed) for run in self.seed_runs),
            "endpoint_padding_fraction": self.endpoint_padding_fraction,
            "stage_generations": self.stage_generations,
            "nested_barrier_points": ",".join(str(value) for value in self.nested_barrier_points),
            "maximum_normalized_bracket_width": self.maximum_normalized_bracket_width,
            "canonical_bistable_barrier_lower": None
            if self.canonical_bistable_barrier_interval is None
            else self.canonical_bistable_barrier_interval[0],
            "canonical_bistable_barrier_upper": None
            if self.canonical_bistable_barrier_interval is None
            else self.canonical_bistable_barrier_interval[1],
        }
        row.update(asdict(self.parameters))
        row.update(_flatten_mapping(self.summary))
        return row


def run_finite_h1_boundary_seed_ensemble(
    spec: ExperimentSpec,
    *,
    master_seeds: Sequence[int] = DEFAULT_ENSEMBLE_MASTER_SEEDS,
    endpoint_padding_fraction: float = 0.5,
    stage_generations: int = 30,
    nested_barrier_points: Sequence[int] = DEFAULT_NESTED_BARRIER_POINTS,
    interaction_separation_threshold: float = 0.05,
    low_state_threshold: float = 0.25,
    high_state_threshold: float = 0.75,
    maximum_normalized_bracket_width: float = 0.03,
) -> tuple[FiniteH1BoundarySeedEnsembleCell, ...]:
    """Run the same nested-grid boundary audit across independent master seeds."""
    seeds = _validate_master_seeds(master_seeds)
    campaigns = tuple(
        run_finite_h1_boundary_resolution_audit(
            replace(spec, master_seed=seed),
            endpoint_padding_fraction=endpoint_padding_fraction,
            stage_generations=stage_generations,
            nested_barrier_points=nested_barrier_points,
            interaction_separation_threshold=interaction_separation_threshold,
            low_state_threshold=low_state_threshold,
            high_state_threshold=high_state_threshold,
            maximum_normalized_bracket_width=maximum_normalized_bracket_width,
        )
        for seed in seeds
    )
    expected_count = len(campaigns[0])
    if any(len(campaign) != expected_count for campaign in campaigns):
        raise RuntimeError("master-seed campaigns produced unequal parameter-cell counts")

    ensemble: list[FiniteH1BoundarySeedEnsembleCell] = []
    for cell_index in range(expected_count):
        source = tuple(campaign[cell_index] for campaign in campaigns)
        reference = source[0]
        identity = _cell_identity(reference)
        if any(_cell_identity(cell) != identity for cell in source[1:]):
            raise RuntimeError("master-seed campaigns produced misaligned parameter cells")
        seed_runs = tuple(
            FiniteH1BoundarySeedRun(master_seed=seed, cell=cell)
            for seed, cell in zip(seeds, source, strict=True)
        )
        ensemble.append(
            FiniteH1BoundarySeedEnsembleCell(
                experiment_id=reference.experiment_id,
                profile=reference.profile,
                pair_index=reference.pair_index,
                parameters=reference.parameters,
                canonical_bistable_barrier_interval=reference.canonical_bistable_barrier_interval,
                endpoint_padding_fraction=reference.endpoint_padding_fraction,
                stage_generations=reference.stage_generations,
                nested_barrier_points=reference.nested_barrier_points,
                maximum_normalized_bracket_width=reference.maximum_normalized_bracket_width,
                seed_runs=seed_runs,
                summary=_summarise_seed_runs(seed_runs),
            )
        )
    return tuple(ensemble)


def write_finite_h1_boundary_seed_ensemble_artifacts(
    cells: Iterable[FiniteH1BoundarySeedEnsembleCell],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> None:
    """Write complete seed-specific records and an aligned ensemble ledger."""
    values = tuple(cells)
    if not values:
        raise ValueError("cells must be nonempty")
    csv_target = Path(csv_path)
    json_target = Path(json_path)
    csv_target.parent.mkdir(parents=True, exist_ok=True)
    json_target.parent.mkdir(parents=True, exist_ok=True)
    rows = [cell.to_csv_row() for cell in values]
    fieldnames = sorted({key for row in rows for key in row})
    with csv_target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    with json_target.open("w", encoding="utf-8") as handle:
        json.dump([cell.as_dict() for cell in values], handle, indent=2, sort_keys=True)


def _summarise_seed_runs(seed_runs: Sequence[FiniteH1BoundarySeedRun]) -> dict[str, object]:
    if not seed_runs:
        raise ValueError("seed_runs must be nonempty")
    per_seed: dict[str, object] = {}
    stable_probabilities: list[float] = []
    mechanism_probabilities: list[float] = []
    pooled_stable: list[bool] = []
    pooled_mechanism: list[bool] = []
    pooled_available: list[bool] = []

    for run in seed_runs:
        cell_summary = run.cell.summary
        stable = cell_summary.get("resolution_stable_loop_supported_probability")
        mechanism = cell_summary.get("resolution_stable_h1_loop_mechanism_supported_probability")
        available = cell_summary.get("available_nested_grid_probability")
        per_seed[str(run.master_seed)] = {
            "available_nested_grid_probability": available,
            "resolution_stable_loop_supported_probability": stable,
            "resolution_stable_h1_loop_mechanism_supported_probability": mechanism,
            "loop_on_two_finest_grids_probability": cell_summary.get("loop_on_two_finest_grids_probability"),
            "boundary_location_stable_probability": cell_summary.get("boundary_location_stable_probability"),
            "finest_grid_bracket_width_fraction": cell_summary.get("finest_grid_bracket_width_fraction"),
        }
        if stable is not None:
            stable_probabilities.append(float(stable))
        if mechanism is not None:
            mechanism_probabilities.append(float(mechanism))
        for replicate in run.cell.replicates:
            if replicate.resolution_stable_loop_supported is not None:
                pooled_available.append(True)
                pooled_stable.append(bool(replicate.resolution_stable_loop_supported))
                pooled_mechanism.append(bool(replicate.resolution_stable_h1_loop_mechanism_supported))
            else:
                pooled_available.append(False)

    summary: dict[str, object] = {
        "master_seed_count": len(seed_runs),
        "replicates_per_master_seed": len(seed_runs[0].cell.replicates),
        "total_replicates": sum(len(run.cell.replicates) for run in seed_runs),
        "by_master_seed": per_seed,
        "available_master_seed_probability": _probability(
            value is not None
            for value in (
                run.cell.summary.get("resolution_stable_loop_supported_probability")
                for run in seed_runs
            )
        ),
        "pooled_available_replicate_probability": _probability(pooled_available),
        "pooled_resolution_stable_loop_supported_probability": None
        if not pooled_stable
        else _probability(pooled_stable),
        "pooled_resolution_stable_h1_loop_mechanism_supported_probability": None
        if not pooled_mechanism
        else _probability(pooled_mechanism),
        "master_seed_resolution_stable_loop_supported_probability": _summary(stable_probabilities),
        "master_seed_resolution_stable_h1_loop_mechanism_supported_probability": _summary(mechanism_probabilities),
        "all_master_seed_runs_fully_resolution_stable": None
        if len(stable_probabilities) != len(seed_runs)
        else all(value == 1.0 for value in stable_probabilities),
        "all_master_seed_runs_fully_resolution_stable_h1_mechanism": None
        if len(mechanism_probabilities) != len(seed_runs)
        else all(value == 1.0 for value in mechanism_probabilities),
    }
    return summary


def _validate_master_seeds(values: Sequence[int]) -> tuple[int, ...]:
    seeds = tuple(int(value) for value in values)
    if len(seeds) < 2:
        raise ValueError("master_seeds must contain at least two independent seeds")
    if len(set(seeds)) != len(seeds):
        raise ValueError("master_seeds must be distinct")
    if any(seed < 0 for seed in seeds):
        raise ValueError("master_seeds must be non-negative")
    return seeds


def _cell_identity(cell: FiniteH1BoundaryResolutionCell) -> str:
    return json.dumps(
        {
            "pair_index": cell.pair_index,
            "parameters": asdict(cell.parameters),
            "interval": cell.canonical_bistable_barrier_interval,
            "endpoint_padding_fraction": cell.endpoint_padding_fraction,
            "stage_generations": cell.stage_generations,
            "nested_barrier_points": cell.nested_barrier_points,
            "maximum_normalized_bracket_width": cell.maximum_normalized_bracket_width,
        },
        sort_keys=True,
    )


def _probability(values: Iterable[bool]) -> float:
    observed = tuple(values)
    if not observed:
        raise ValueError("values must be nonempty")
    return sum(observed) / len(observed)


def _summary(values: Iterable[float]) -> dict[str, float | None]:
    observed = tuple(float(value) for value in values)
    if not observed:
        return _empty_summary()
    return {
        "mean": sum(observed) / len(observed),
        "median": median(observed),
        "minimum": min(observed),
        "maximum": max(observed),
    }


def _empty_summary() -> dict[str, None]:
    return {"mean": None, "median": None, "minimum": None, "maximum": None}


def _flatten_mapping(mapping: Mapping[str, object], prefix: str = "") -> dict[str, object]:
    flat: dict[str, object] = {}
    for key, value in mapping.items():
        name = f"{prefix}{key}" if not prefix else f"{prefix}.{key}"
        if isinstance(value, Mapping):
            flat.update(_flatten_mapping(value, name))
        else:
            flat[name] = value
    return flat
