"""Stage-duration robustness audit for finite H1 continuation results.

A finite rising/falling continuation can show route memory because the stage is
shorter than the relaxation time, not because the route-specific states persist
under a quasi-static barrier change.  The original finite H1 hysteresis audit
uses one declared number of generations per barrier stage.  This module repeats
that *same* one-large continuation protocol across an ordered stage-duration
ladder and preserves the full result at every duration.

The central finite result is deliberately narrow:

``convergence_robust_hysteresis_supported`` requires route-memory support at the
two longest declared durations and an absolute change in maximum internal route
gap no larger than a declared stability tolerance.  It does not prove a theorem
or equilibrium convergence.  It is a Type S check against the specific concern
that a short continuation stage manufactured apparent hysteresis.
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


@dataclass(frozen=True)
class H1HysteresisDurationObservation:
    """One replicate's finite continuation outcome at one stage duration."""

    stage_generations: int
    finite_hysteresis_supported: bool | None
    finite_h1_hysteresis_mechanism_supported: bool | None
    maximum_internal_interaction_gap: float | None
    jump_boundary_gap: float | None
    internal_potential_switch_barrier_count: int | None

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class H1HysteresisDurationReplicate:
    """Same-seed stage-duration ladder for one one-large parameter replicate."""

    replicate_index: int
    seed: int
    observations: tuple[H1HysteresisDurationObservation, ...]
    longest_pair_gap_change: float | None
    longest_pair_gap_stable: bool | None
    support_at_all_durations: bool | None
    mechanism_support_at_all_durations: bool | None
    convergence_robust_hysteresis_supported: bool | None
    convergence_robust_h1_mechanism_supported: bool | None

    def as_dict(self) -> dict[str, object]:
        return {
            "replicate_index": self.replicate_index,
            "seed": self.seed,
            "observations": [observation.as_dict() for observation in self.observations],
            "longest_pair_gap_change": self.longest_pair_gap_change,
            "longest_pair_gap_stable": self.longest_pair_gap_stable,
            "support_at_all_durations": self.support_at_all_durations,
            "mechanism_support_at_all_durations": self.mechanism_support_at_all_durations,
            "convergence_robust_hysteresis_supported": self.convergence_robust_hysteresis_supported,
            "convergence_robust_h1_mechanism_supported": self.convergence_robust_h1_mechanism_supported,
        }


@dataclass(frozen=True)
class H1HysteresisDurationCell:
    """One (A_ref, kappa) one-large continuation robustness family."""

    experiment_id: str
    profile: str
    pair_index: int
    parameters: ParameterCell
    canonical_bistable_barrier_interval: tuple[float, float] | None
    stage_generations: tuple[int, ...]
    barrier_points: int
    barrier_padding: float
    interaction_separation_threshold: float
    gap_stability_tolerance: float
    duration_cells: tuple[FiniteH1HysteresisCell, ...]
    replicates: tuple[H1HysteresisDurationReplicate, ...]
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
            "stage_generations": list(self.stage_generations),
            "barrier_points": self.barrier_points,
            "barrier_padding": self.barrier_padding,
            "interaction_separation_threshold": self.interaction_separation_threshold,
            "gap_stability_tolerance": self.gap_stability_tolerance,
            "replicate_count": len(self.replicates),
            "duration_cells": [cell.as_dict() for cell in self.duration_cells],
            "replicates": [replicate.as_dict() for replicate in self.replicates],
            "summary": dict(self.summary),
        }

    def to_csv_row(self) -> dict[str, object]:
        row: dict[str, object] = {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "pair_index": self.pair_index,
            "replicate_count": len(self.replicates),
            "stage_generations": ",".join(str(value) for value in self.stage_generations),
            "barrier_points": self.barrier_points,
            "barrier_padding": self.barrier_padding,
            "interaction_separation_threshold": self.interaction_separation_threshold,
            "gap_stability_tolerance": self.gap_stability_tolerance,
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


def run_finite_h1_hysteresis_duration_audit(
    spec: ExperimentSpec,
    *,
    stage_generations: Sequence[int] = (5, 10, 30, 80),
    barrier_points: int = 7,
    barrier_padding: float = 0.1,
    interaction_separation_threshold: float = 0.05,
    low_state_threshold: float = 0.25,
    high_state_threshold: float = 0.75,
    gap_stability_tolerance: float = 0.05,
) -> tuple[H1HysteresisDurationCell, ...]:
    """Repeat one-large continuation over an ordered stage-duration ladder.

    Exactly one barrier is used per ``(area_reference, interaction_feedback)``
    pair: the midpoint of its one-large strict canonical bistable interval.
    The finite continuation itself spans the full interval, so the arbitrary
    raw-barrier grids used by broad phase diagrams do not duplicate the same
    hysteresis path. Pairs without a strict canonical interval are retained as
    unavailable rather than being silently omitted.
    """
    durations = _validate_stage_generations(stage_generations)
    if barrier_points < 3:
        raise ValueError("barrier_points must be at least three")
    if barrier_padding <= 0.0:
        raise ValueError("barrier_padding must be positive")
    if interaction_separation_threshold < 0.0:
        raise ValueError("interaction_separation_threshold must be non-negative")
    if gap_stability_tolerance < 0.0:
        raise ValueError("gap_stability_tolerance must be non-negative")

    cells: list[H1HysteresisDurationCell] = []
    pair_index = 0
    for area_reference in spec.area_reference_values:
        for interaction_feedback in spec.interaction_feedback_values:
            interval = canonical_bistable_barrier_interval(
                interaction_feedback,
                spec.total_area,
                area_reference,
            )
            midpoint = (
                spec.interaction_barrier_values[0]
                if interval is None
                else (interval[0] + interval[1]) / 2.0
            )
            pair_seed = _pair_seed(spec.master_seed, pair_index)
            pair_spec = replace(
                spec,
                area_reference_values=(area_reference,),
                interaction_feedback_values=(interaction_feedback,),
                interaction_barrier_values=(midpoint,),
                master_seed=pair_seed,
            )
            scenario = scenario_one_large(pair_spec)
            duration_cells = tuple(
                run_finite_h1_hysteresis_audit(
                    pair_spec,
                    scenarios=(scenario,),
                    barrier_points=barrier_points,
                    barrier_padding=barrier_padding,
                    stage_generations=duration,
                    interaction_separation_threshold=interaction_separation_threshold,
                    low_state_threshold=low_state_threshold,
                    high_state_threshold=high_state_threshold,
                )[0]
                for duration in durations
            )
            parameter_cell = duration_cells[0].parameters
            replicates = _build_replicates(
                duration_cells,
                durations=durations,
                gap_stability_tolerance=gap_stability_tolerance,
            )
            cells.append(
                H1HysteresisDurationCell(
                    experiment_id=spec.experiment_id,
                    profile=spec.profile,
                    pair_index=pair_index,
                    parameters=parameter_cell,
                    canonical_bistable_barrier_interval=interval,
                    stage_generations=durations,
                    barrier_points=barrier_points,
                    barrier_padding=barrier_padding,
                    interaction_separation_threshold=interaction_separation_threshold,
                    gap_stability_tolerance=gap_stability_tolerance,
                    duration_cells=duration_cells,
                    replicates=replicates,
                    summary=_summarise_replicates(replicates, durations),
                )
            )
            pair_index += 1
    return tuple(cells)


def write_finite_h1_hysteresis_duration_artifacts(
    cells: Iterable[H1HysteresisDurationCell],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> None:
    """Write complete duration-ladder trajectories and flat summaries."""
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
    duration_cells: Sequence[FiniteH1HysteresisCell],
    *,
    durations: tuple[int, ...],
    gap_stability_tolerance: float,
) -> tuple[H1HysteresisDurationReplicate, ...]:
    records_by_duration = [cell.replicates for cell in duration_cells]
    count = len(records_by_duration[0])
    if any(len(records) != count for records in records_by_duration):
        raise RuntimeError("stage-duration runs produced unequal replicate counts")
    output: list[H1HysteresisDurationReplicate] = []
    for replicate_index in range(count):
        source = tuple(records[replicate_index] for records in records_by_duration)
        if any(record.replicate_index != replicate_index for record in source):
            raise RuntimeError("stage-duration replicate order is inconsistent")
        if len({record.seed for record in source}) != 1:
            raise RuntimeError("same-seed duration pairing was not preserved")
        observations = tuple(
            _observation(duration, record) for duration, record in zip(durations, source, strict=True)
        )
        available = all(observation.finite_hysteresis_supported is not None for observation in observations)
        if not available:
            output.append(
                H1HysteresisDurationReplicate(
                    replicate_index=replicate_index,
                    seed=source[0].seed,
                    observations=observations,
                    longest_pair_gap_change=None,
                    longest_pair_gap_stable=None,
                    support_at_all_durations=None,
                    mechanism_support_at_all_durations=None,
                    convergence_robust_hysteresis_supported=None,
                    convergence_robust_h1_mechanism_supported=None,
                )
            )
            continue
        largest_pair = observations[-2:]
        previous_gap = _required(largest_pair[0].maximum_internal_interaction_gap)
        largest_gap = _required(largest_pair[1].maximum_internal_interaction_gap)
        gap_change = abs(largest_gap - previous_gap)
        gap_stable = gap_change <= gap_stability_tolerance
        support_all = all(bool(observation.finite_hysteresis_supported) for observation in observations)
        mechanism_all = all(bool(observation.finite_h1_hysteresis_mechanism_supported) for observation in observations)
        output.append(
            H1HysteresisDurationReplicate(
                replicate_index=replicate_index,
                seed=source[0].seed,
                observations=observations,
                longest_pair_gap_change=gap_change,
                longest_pair_gap_stable=gap_stable,
                support_at_all_durations=support_all,
                mechanism_support_at_all_durations=mechanism_all,
                convergence_robust_hysteresis_supported=(
                    bool(largest_pair[0].finite_hysteresis_supported)
                    and bool(largest_pair[1].finite_hysteresis_supported)
                    and gap_stable
                ),
                convergence_robust_h1_mechanism_supported=(
                    bool(largest_pair[0].finite_h1_hysteresis_mechanism_supported)
                    and bool(largest_pair[1].finite_h1_hysteresis_mechanism_supported)
                    and gap_stable
                ),
            )
        )
    return tuple(output)


def _observation(
    duration: int,
    replicate: FiniteH1HysteresisReplicate,
) -> H1HysteresisDurationObservation:
    return H1HysteresisDurationObservation(
        stage_generations=duration,
        finite_hysteresis_supported=replicate.finite_hysteresis_supported,
        finite_h1_hysteresis_mechanism_supported=replicate.finite_h1_hysteresis_mechanism_supported,
        maximum_internal_interaction_gap=replicate.maximum_internal_interaction_gap,
        jump_boundary_gap=replicate.jump_boundary_gap,
        internal_potential_switch_barrier_count=None
        if replicate.internal_potential_switch_barriers is None
        else len(replicate.internal_potential_switch_barriers),
    )


def _summarise_replicates(
    replicates: Sequence[H1HysteresisDurationReplicate],
    durations: Sequence[int],
) -> dict[str, object]:
    if not replicates:
        raise ValueError("replicates must be nonempty")
    available = tuple(
        replicate for replicate in replicates if replicate.convergence_robust_hysteresis_supported is not None
    )
    summary: dict[str, object] = {
        "replicate_count": len(replicates),
        "available_duration_ladder_probability": len(available) / len(replicates),
        "by_stage_generations": {},
    }
    for index, duration in enumerate(durations):
        observations = tuple(replicate.observations[index] for replicate in available)
        valid_gap = tuple(
            observation.maximum_internal_interaction_gap
            for observation in observations
            if observation.maximum_internal_interaction_gap is not None
        )
        summary["by_stage_generations"][str(duration)] = {
            "finite_hysteresis_supported_probability": None
            if not observations
            else _probability(bool(observation.finite_hysteresis_supported) for observation in observations),
            "finite_h1_hysteresis_mechanism_supported_probability": None
            if not observations
            else _probability(bool(observation.finite_h1_hysteresis_mechanism_supported) for observation in observations),
            "maximum_internal_interaction_gap": _summary(valid_gap),
            "jump_boundary_available_probability": None
            if not observations
            else _probability(observation.jump_boundary_gap is not None for observation in observations),
        }
    if not available:
        summary.update(
            {
                "support_at_all_durations_probability": None,
                "mechanism_support_at_all_durations_probability": None,
                "longest_pair_gap_change": _empty_summary(),
                "longest_pair_gap_stable_probability": None,
                "convergence_robust_hysteresis_supported_probability": None,
                "convergence_robust_h1_mechanism_supported_probability": None,
            }
        )
        return summary
    summary.update(
        {
            "support_at_all_durations_probability": _probability(
                bool(replicate.support_at_all_durations) for replicate in available
            ),
            "mechanism_support_at_all_durations_probability": _probability(
                bool(replicate.mechanism_support_at_all_durations) for replicate in available
            ),
            "longest_pair_gap_change": _summary(
                _required(replicate.longest_pair_gap_change) for replicate in available
            ),
            "longest_pair_gap_stable_probability": _probability(
                bool(replicate.longest_pair_gap_stable) for replicate in available
            ),
            "convergence_robust_hysteresis_supported_probability": _probability(
                bool(replicate.convergence_robust_hysteresis_supported) for replicate in available
            ),
            "convergence_robust_h1_mechanism_supported_probability": _probability(
                bool(replicate.convergence_robust_h1_mechanism_supported) for replicate in available
            ),
        }
    )
    return summary


def _validate_stage_generations(values: Sequence[int]) -> tuple[int, ...]:
    durations = tuple(int(value) for value in values)
    if len(durations) < 2:
        raise ValueError("stage_generations must contain at least two durations")
    if any(value < 1 for value in durations):
        raise ValueError("stage generations must be positive")
    if tuple(sorted(durations)) != durations or len(set(durations)) != len(durations):
        raise ValueError("stage generations must be strictly increasing")
    return durations


def _pair_seed(master_seed: int, pair_index: int) -> int:
    return (master_seed * 1_000_003 + pair_index * 10_007 + 271) % (2**31 - 1)


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
        raise RuntimeError("unexpected missing duration-ladder value")
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
