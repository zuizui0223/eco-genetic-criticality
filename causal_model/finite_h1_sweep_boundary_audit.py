"""Endpoint-expansion audit for finite H1 continuation boundaries.

The finite H1 duration ladder can establish that rising/falling route memory
persists at long stage durations.  It may still fail to observe a finite
collapse/recovery pair when the continuation begins and ends only a small
absolute distance beyond the canonical bistable interval.  This audit expands
those endpoints in *fractions of the canonical interval width* and records the
first sweep that observes both finite route transitions.

A detected transition is a grid-limited finite estimate, not an exact critical
barrier.  Every estimate carries the barrier-grid step that bounds its
resolution.  The audit is Type S and does not prove a finite bifurcation.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from statistics import median
from typing import Iterable, Mapping, Sequence

from causal_model.canonical_h1_bifurcation import canonical_bistable_barrier_interval
from causal_model.finite_h1_hysteresis_audit import (
    FiniteH1HysteresisCell,
    FiniteH1HysteresisReplicate,
    run_finite_h1_hysteresis_audit,
)
from causal_model.multipatch_criticality_experiments import (
    ExperimentSpec,
    ParameterCell,
    scenario_one_large,
)

DEFAULT_ENDPOINT_PADDING_FRACTIONS = (0.1, 0.5, 1.0, 2.0)


@dataclass(frozen=True)
class FiniteH1SweepObservation:
    """One endpoint-expanded continuation sweep for one same-seed replicate."""

    endpoint_padding_fraction: float
    absolute_padding: float
    barrier_step: float
    finite_hysteresis_supported: bool | None
    finite_h1_hysteresis_mechanism_supported: bool | None
    rising_collapse_barrier: float | None
    falling_recovery_barrier: float | None
    jump_boundary_gap: float | None
    finite_loop_closed: bool | None
    finite_h1_loop_mechanism_supported: bool | None

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class FiniteH1SweepReplicate:
    """Same-seed endpoint-expansion ladder for one one-large H1 replicate."""

    replicate_index: int
    seed: int
    observations: tuple[FiniteH1SweepObservation, ...]
    first_closed_padding_fraction: float | None
    first_closed_barrier_step: float | None
    first_rising_collapse_barrier: float | None
    first_falling_recovery_barrier: float | None
    first_jump_boundary_gap: float | None
    finite_loop_closed_at_any_padding: bool | None
    finite_h1_loop_mechanism_supported_at_any_padding: bool | None

    def as_dict(self) -> dict[str, object]:
        return {
            "replicate_index": self.replicate_index,
            "seed": self.seed,
            "observations": [observation.as_dict() for observation in self.observations],
            "first_closed_padding_fraction": self.first_closed_padding_fraction,
            "first_closed_barrier_step": self.first_closed_barrier_step,
            "first_rising_collapse_barrier": self.first_rising_collapse_barrier,
            "first_falling_recovery_barrier": self.first_falling_recovery_barrier,
            "first_jump_boundary_gap": self.first_jump_boundary_gap,
            "finite_loop_closed_at_any_padding": self.finite_loop_closed_at_any_padding,
            "finite_h1_loop_mechanism_supported_at_any_padding": self.finite_h1_loop_mechanism_supported_at_any_padding,
        }


@dataclass(frozen=True)
class FiniteH1SweepBoundaryCell:
    """Endpoint-expanded finite boundary evidence for one (A_ref, kappa) pair."""

    experiment_id: str
    profile: str
    pair_index: int
    parameters: ParameterCell
    canonical_bistable_barrier_interval: tuple[float, float] | None
    endpoint_padding_fractions: tuple[float, ...]
    stage_generations: int
    barrier_points: int
    interaction_separation_threshold: float
    low_state_threshold: float
    high_state_threshold: float
    sweep_cells: tuple[FiniteH1HysteresisCell, ...]
    replicates: tuple[FiniteH1SweepReplicate, ...]
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
            "endpoint_padding_fractions": list(self.endpoint_padding_fractions),
            "stage_generations": self.stage_generations,
            "barrier_points": self.barrier_points,
            "interaction_separation_threshold": self.interaction_separation_threshold,
            "low_state_threshold": self.low_state_threshold,
            "high_state_threshold": self.high_state_threshold,
            "replicate_count": len(self.replicates),
            "sweep_cells": [cell.as_dict() for cell in self.sweep_cells],
            "replicates": [replicate.as_dict() for replicate in self.replicates],
            "summary": dict(self.summary),
        }

    def to_csv_row(self) -> dict[str, object]:
        row: dict[str, object] = {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "pair_index": self.pair_index,
            "replicate_count": len(self.replicates),
            "endpoint_padding_fractions": ",".join(str(value) for value in self.endpoint_padding_fractions),
            "stage_generations": self.stage_generations,
            "barrier_points": self.barrier_points,
            "interaction_separation_threshold": self.interaction_separation_threshold,
            "low_state_threshold": self.low_state_threshold,
            "high_state_threshold": self.high_state_threshold,
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


def run_finite_h1_sweep_boundary_audit(
    spec: ExperimentSpec,
    *,
    endpoint_padding_fractions: Sequence[float] = DEFAULT_ENDPOINT_PADDING_FRACTIONS,
    stage_generations: int = 30,
    barrier_points: int = 25,
    interaction_separation_threshold: float = 0.05,
    low_state_threshold: float = 0.25,
    high_state_threshold: float = 0.75,
) -> tuple[FiniteH1SweepBoundaryCell, ...]:
    """Expand finite continuation endpoints until route transitions are observed.

    One midpoint-labelled one-large continuation family is run for every
    ``(area_reference, interaction_feedback)`` pair.  The endpoint padding for
    one sweep is ``fraction * (canonical_upper - canonical_lower)``.  Results
    at all fractions are retained; no finite outcome chooses which fractions
    appear in the design.
    """
    fractions = _validate_padding_fractions(endpoint_padding_fractions)
    if stage_generations < 1:
        raise ValueError("stage_generations must be positive")
    if barrier_points < 3:
        raise ValueError("barrier_points must be at least three")
    if interaction_separation_threshold < 0.0:
        raise ValueError("interaction_separation_threshold must be non-negative")
    if not 0.0 <= low_state_threshold < high_state_threshold <= 1.0:
        raise ValueError("state thresholds must satisfy 0 <= low < high <= 1")

    cells: list[FiniteH1SweepBoundaryCell] = []
    pair_index = 0
    for area_reference in spec.area_reference_values:
        for interaction_feedback in spec.interaction_feedback_values:
            interval = canonical_bistable_barrier_interval(
                interaction_feedback,
                spec.total_area,
                area_reference,
            )
            midpoint = spec.interaction_barrier_values[0] if interval is None else sum(interval) / 2.0
            pair_spec = replace(
                spec,
                area_reference_values=(area_reference,),
                interaction_feedback_values=(interaction_feedback,),
                interaction_barrier_values=(midpoint,),
                master_seed=_pair_seed(spec.master_seed, pair_index),
            )
            scenario = scenario_one_large(pair_spec)
            sweep_cells = tuple(
                run_finite_h1_hysteresis_audit(
                    pair_spec,
                    scenarios=(scenario,),
                    barrier_points=barrier_points,
                    barrier_padding=_absolute_padding(interval, fraction),
                    stage_generations=stage_generations,
                    interaction_separation_threshold=interaction_separation_threshold,
                    low_state_threshold=low_state_threshold,
                    high_state_threshold=high_state_threshold,
                )[0]
                for fraction in fractions
            )
            replicates = _build_replicates(
                sweep_cells,
                interval=interval,
                fractions=fractions,
                barrier_points=barrier_points,
            )
            cells.append(
                FiniteH1SweepBoundaryCell(
                    experiment_id=spec.experiment_id,
                    profile=spec.profile,
                    pair_index=pair_index,
                    parameters=sweep_cells[0].parameters,
                    canonical_bistable_barrier_interval=interval,
                    endpoint_padding_fractions=fractions,
                    stage_generations=stage_generations,
                    barrier_points=barrier_points,
                    interaction_separation_threshold=interaction_separation_threshold,
                    low_state_threshold=low_state_threshold,
                    high_state_threshold=high_state_threshold,
                    sweep_cells=sweep_cells,
                    replicates=replicates,
                    summary=_summarise(replicates, fractions),
                )
            )
            pair_index += 1
    return tuple(cells)


def write_finite_h1_sweep_boundary_artifacts(
    cells: Iterable[FiniteH1SweepBoundaryCell],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> None:
    """Write complete endpoint-expansion sweeps and a flat cell summary."""
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


def _build_replicates(
    sweep_cells: Sequence[FiniteH1HysteresisCell],
    *,
    interval: tuple[float, float] | None,
    fractions: tuple[float, ...],
    barrier_points: int,
) -> tuple[FiniteH1SweepReplicate, ...]:
    by_fraction = [cell.replicates for cell in sweep_cells]
    count = len(by_fraction[0])
    if any(len(values) != count for values in by_fraction):
        raise RuntimeError("endpoint sweeps produced unequal replicate counts")
    output: list[FiniteH1SweepReplicate] = []
    for replicate_index in range(count):
        source = tuple(values[replicate_index] for values in by_fraction)
        if any(record.replicate_index != replicate_index for record in source):
            raise RuntimeError("endpoint sweep replicate order is inconsistent")
        if len({record.seed for record in source}) != 1:
            raise RuntimeError("same-seed endpoint pairing was not preserved")
        observations = tuple(
            _observation(
                record,
                endpoint_padding_fraction=fraction,
                absolute_padding=_absolute_padding(interval, fraction),
                barrier_step=_barrier_step(interval, fraction, barrier_points),
            )
            for fraction, record in zip(fractions, source, strict=True)
        )
        if any(observation.finite_loop_closed is None for observation in observations):
            output.append(
                FiniteH1SweepReplicate(
                    replicate_index=replicate_index,
                    seed=source[0].seed,
                    observations=observations,
                    first_closed_padding_fraction=None,
                    first_closed_barrier_step=None,
                    first_rising_collapse_barrier=None,
                    first_falling_recovery_barrier=None,
                    first_jump_boundary_gap=None,
                    finite_loop_closed_at_any_padding=None,
                    finite_h1_loop_mechanism_supported_at_any_padding=None,
                )
            )
            continue
        first = next((item for item in observations if item.finite_loop_closed), None)
        first_mechanism = next((item for item in observations if item.finite_h1_loop_mechanism_supported), None)
        output.append(
            FiniteH1SweepReplicate(
                replicate_index=replicate_index,
                seed=source[0].seed,
                observations=observations,
                first_closed_padding_fraction=None if first is None else first.endpoint_padding_fraction,
                first_closed_barrier_step=None if first is None else first.barrier_step,
                first_rising_collapse_barrier=None if first is None else first.rising_collapse_barrier,
                first_falling_recovery_barrier=None if first is None else first.falling_recovery_barrier,
                first_jump_boundary_gap=None if first is None else first.jump_boundary_gap,
                finite_loop_closed_at_any_padding=first is not None,
                finite_h1_loop_mechanism_supported_at_any_padding=first_mechanism is not None,
            )
        )
    return tuple(output)


def _observation(
    replicate: FiniteH1HysteresisReplicate,
    *,
    endpoint_padding_fraction: float,
    absolute_padding: float,
    barrier_step: float,
) -> FiniteH1SweepObservation:
    available = replicate.finite_hysteresis_supported is not None
    loop_closed = (
        None
        if not available
        else replicate.rising_collapse_barrier is not None
        and replicate.falling_recovery_barrier is not None
        and _required(replicate.jump_boundary_gap) > 0.0
    )
    mechanism_closed = (
        None
        if loop_closed is None
        else bool(loop_closed) and bool(replicate.finite_h1_hysteresis_mechanism_supported)
    )
    return FiniteH1SweepObservation(
        endpoint_padding_fraction=endpoint_padding_fraction,
        absolute_padding=absolute_padding,
        barrier_step=barrier_step,
        finite_hysteresis_supported=replicate.finite_hysteresis_supported,
        finite_h1_hysteresis_mechanism_supported=replicate.finite_h1_hysteresis_mechanism_supported,
        rising_collapse_barrier=replicate.rising_collapse_barrier,
        falling_recovery_barrier=replicate.falling_recovery_barrier,
        jump_boundary_gap=replicate.jump_boundary_gap,
        finite_loop_closed=loop_closed,
        finite_h1_loop_mechanism_supported=mechanism_closed,
    )


def _summarise(
    replicates: Sequence[FiniteH1SweepReplicate],
    fractions: Sequence[float],
) -> dict[str, object]:
    if not replicates:
        raise ValueError("replicates must be nonempty")
    available = tuple(rep for rep in replicates if rep.finite_loop_closed_at_any_padding is not None)
    summary: dict[str, object] = {
        "replicate_count": len(replicates),
        "available_endpoint_sweep_probability": len(available) / len(replicates),
        "by_endpoint_padding_fraction": {},
    }
    for index, fraction in enumerate(fractions):
        observations = tuple(rep.observations[index] for rep in available)
        loop_closed = tuple(obs for obs in observations if obs.finite_loop_closed)
        summary["by_endpoint_padding_fraction"][str(fraction)] = {
            "finite_hysteresis_supported_probability": None
            if not observations
            else _probability(bool(obs.finite_hysteresis_supported) for obs in observations),
            "finite_h1_hysteresis_mechanism_supported_probability": None
            if not observations
            else _probability(bool(obs.finite_h1_hysteresis_mechanism_supported) for obs in observations),
            "finite_loop_closed_probability": None
            if not observations
            else _probability(bool(obs.finite_loop_closed) for obs in observations),
            "finite_h1_loop_mechanism_supported_probability": None
            if not observations
            else _probability(bool(obs.finite_h1_loop_mechanism_supported) for obs in observations),
            "rising_collapse_barrier": _summary(
                obs.rising_collapse_barrier for obs in loop_closed if obs.rising_collapse_barrier is not None
            ),
            "falling_recovery_barrier": _summary(
                obs.falling_recovery_barrier for obs in loop_closed if obs.falling_recovery_barrier is not None
            ),
            "jump_boundary_gap": _summary(
                obs.jump_boundary_gap for obs in loop_closed if obs.jump_boundary_gap is not None
            ),
            "barrier_step": _summary(obs.barrier_step for obs in observations),
        }
    if not available:
        summary.update(
            {
                "finite_loop_closed_at_any_padding_probability": None,
                "finite_h1_loop_mechanism_supported_at_any_padding_probability": None,
                "first_closed_padding_fraction": _empty_summary(),
                "first_closed_barrier_step": _empty_summary(),
                "first_rising_collapse_barrier": _empty_summary(),
                "first_falling_recovery_barrier": _empty_summary(),
                "first_jump_boundary_gap": _empty_summary(),
            }
        )
        return summary
    closed = tuple(rep for rep in available if rep.finite_loop_closed_at_any_padding)
    summary.update(
        {
            "finite_loop_closed_at_any_padding_probability": _probability(
                bool(rep.finite_loop_closed_at_any_padding) for rep in available
            ),
            "finite_h1_loop_mechanism_supported_at_any_padding_probability": _probability(
                bool(rep.finite_h1_loop_mechanism_supported_at_any_padding) for rep in available
            ),
            "first_closed_padding_fraction": _summary(
                rep.first_closed_padding_fraction for rep in closed if rep.first_closed_padding_fraction is not None
            ),
            "first_closed_barrier_step": _summary(
                rep.first_closed_barrier_step for rep in closed if rep.first_closed_barrier_step is not None
            ),
            "first_rising_collapse_barrier": _summary(
                rep.first_rising_collapse_barrier for rep in closed if rep.first_rising_collapse_barrier is not None
            ),
            "first_falling_recovery_barrier": _summary(
                rep.first_falling_recovery_barrier for rep in closed if rep.first_falling_recovery_barrier is not None
            ),
            "first_jump_boundary_gap": _summary(
                rep.first_jump_boundary_gap for rep in closed if rep.first_jump_boundary_gap is not None
            ),
        }
    )
    return summary


def _validate_padding_fractions(values: Sequence[float]) -> tuple[float, ...]:
    fractions = tuple(float(value) for value in values)
    if not fractions:
        raise ValueError("endpoint_padding_fractions must be nonempty")
    if any(value <= 0.0 for value in fractions):
        raise ValueError("endpoint padding fractions must be positive")
    if tuple(sorted(fractions)) != fractions or len(set(fractions)) != len(fractions):
        raise ValueError("endpoint padding fractions must be strictly increasing")
    return fractions


def _absolute_padding(interval: tuple[float, float] | None, fraction: float) -> float:
    if interval is None:
        return 0.1
    return fraction * (interval[1] - interval[0])


def _barrier_step(interval: tuple[float, float] | None, fraction: float, points: int) -> float:
    if interval is None:
        return 0.0
    width = interval[1] - interval[0]
    return (width + 2.0 * _absolute_padding(interval, fraction)) / (points - 1)


def _pair_seed(master_seed: int, pair_index: int) -> int:
    return (master_seed * 1_000_003 + pair_index * 10_007 + 401) % (2**31 - 1)


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


def _required(value):
    if value is None:
        raise RuntimeError("unexpected missing finite sweep boundary value")
    return value


def _flatten_mapping(mapping: Mapping[str, object], prefix: str = "") -> dict[str, object]:
    flat: dict[str, object] = {}
    for key, value in mapping.items():
        name = f"{prefix}{key}" if not prefix else f"{prefix}.{key}"
        if isinstance(value, Mapping):
            flat.update(_flatten_mapping(value, name))
        else:
            flat[name] = value
    return flat
