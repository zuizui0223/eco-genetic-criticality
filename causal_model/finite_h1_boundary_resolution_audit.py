"""Nested-grid resolution audit for finite H1 collapse/recovery boundaries.

The endpoint-expansion audit establishes a common continuation range wide enough
to observe finite collapse and recovery.  Its boundary locations are nevertheless
coarse-grid estimates.  This module fixes one common endpoint padding fraction
and repeats the same one-large continuation on *nested* barrier grids.

The default grids (25, 49, 97) have 24, 48, and 96 intervals respectively, so
every coarser grid is a subset of the next finer grid.  For every detected
transition, the audit stores the two adjacent barriers that bracket the first
threshold crossing, rather than treating the first observed crossing grid point
as an exact finite critical value.

``resolution_stable_loop_supported`` is an operational Type S predicate.  It
requires threshold-defined loops on the two finest grids, a finest-grid bracket
no wider than a predeclared fraction of canonical-interval width, and boundary
locations that move by no more than one penultimate-grid step.  It is evidence
of grid-resolution robustness for this closure, not a bifurcation theorem.
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
    FiniteH1ContinuationStage,
    FiniteH1HysteresisCell,
    FiniteH1HysteresisReplicate,
    run_finite_h1_hysteresis_audit,
)
from causal_model.multipatch_criticality_experiments import (
    ExperimentSpec,
    ParameterCell,
    scenario_one_large,
)

DEFAULT_NESTED_BARRIER_POINTS = (25, 49, 97)


@dataclass(frozen=True)
class FiniteH1BoundaryBracket:
    """Adjacent grid barriers bracketing a first finite threshold crossing."""

    lower_barrier: float
    upper_barrier: float
    crossing_barrier: float
    prior_barrier: float

    @property
    def width(self) -> float:
        return self.upper_barrier - self.lower_barrier

    @property
    def midpoint(self) -> float:
        return (self.lower_barrier + self.upper_barrier) / 2.0

    def as_dict(self) -> dict[str, float]:
        return {
            "lower_barrier": self.lower_barrier,
            "upper_barrier": self.upper_barrier,
            "crossing_barrier": self.crossing_barrier,
            "prior_barrier": self.prior_barrier,
            "width": self.width,
            "midpoint": self.midpoint,
        }


@dataclass(frozen=True)
class FiniteH1BoundaryResolutionObservation:
    """One same-seed continuation at one nested barrier-grid resolution."""

    barrier_points: int
    barrier_step: float | None
    finite_hysteresis_supported: bool | None
    finite_h1_hysteresis_mechanism_supported: bool | None
    rising_collapse_bracket: FiniteH1BoundaryBracket | None
    falling_recovery_bracket: FiniteH1BoundaryBracket | None
    finite_loop_bracket_supported: bool | None
    finite_h1_loop_bracket_mechanism_supported: bool | None
    finite_gap_lower_bound: float | None
    finite_gap_upper_bound: float | None

    def as_dict(self) -> dict[str, object]:
        return {
            "barrier_points": self.barrier_points,
            "barrier_step": self.barrier_step,
            "finite_hysteresis_supported": self.finite_hysteresis_supported,
            "finite_h1_hysteresis_mechanism_supported": self.finite_h1_hysteresis_mechanism_supported,
            "rising_collapse_bracket": None
            if self.rising_collapse_bracket is None
            else self.rising_collapse_bracket.as_dict(),
            "falling_recovery_bracket": None
            if self.falling_recovery_bracket is None
            else self.falling_recovery_bracket.as_dict(),
            "finite_loop_bracket_supported": self.finite_loop_bracket_supported,
            "finite_h1_loop_bracket_mechanism_supported": self.finite_h1_loop_bracket_mechanism_supported,
            "finite_gap_lower_bound": self.finite_gap_lower_bound,
            "finite_gap_upper_bound": self.finite_gap_upper_bound,
        }


@dataclass(frozen=True)
class FiniteH1BoundaryResolutionReplicate:
    """One same-seed nested-grid family for one one-large parameter pair."""

    replicate_index: int
    seed: int
    observations: tuple[FiniteH1BoundaryResolutionObservation, ...]
    finest_rising_midpoint_shift: float | None
    finest_falling_midpoint_shift: float | None
    finest_grid_bracket_width_fraction: float | None
    loop_on_two_finest_grids: bool | None
    boundary_location_stable: bool | None
    resolution_stable_loop_supported: bool | None
    resolution_stable_h1_loop_mechanism_supported: bool | None

    def as_dict(self) -> dict[str, object]:
        return {
            "replicate_index": self.replicate_index,
            "seed": self.seed,
            "observations": [item.as_dict() for item in self.observations],
            "finest_rising_midpoint_shift": self.finest_rising_midpoint_shift,
            "finest_falling_midpoint_shift": self.finest_falling_midpoint_shift,
            "finest_grid_bracket_width_fraction": self.finest_grid_bracket_width_fraction,
            "loop_on_two_finest_grids": self.loop_on_two_finest_grids,
            "boundary_location_stable": self.boundary_location_stable,
            "resolution_stable_loop_supported": self.resolution_stable_loop_supported,
            "resolution_stable_h1_loop_mechanism_supported": self.resolution_stable_h1_loop_mechanism_supported,
        }


@dataclass(frozen=True)
class FiniteH1BoundaryResolutionCell:
    """Nested-grid finite boundary evidence for one (A_ref, kappa) pair."""

    experiment_id: str
    profile: str
    pair_index: int
    parameters: ParameterCell
    canonical_bistable_barrier_interval: tuple[float, float] | None
    endpoint_padding_fraction: float
    absolute_padding: float
    stage_generations: int
    nested_barrier_points: tuple[int, ...]
    interaction_separation_threshold: float
    low_state_threshold: float
    high_state_threshold: float
    maximum_normalized_bracket_width: float
    resolution_cells: tuple[FiniteH1HysteresisCell, ...]
    replicates: tuple[FiniteH1BoundaryResolutionReplicate, ...]
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
            "absolute_padding": self.absolute_padding,
            "stage_generations": self.stage_generations,
            "nested_barrier_points": list(self.nested_barrier_points),
            "interaction_separation_threshold": self.interaction_separation_threshold,
            "low_state_threshold": self.low_state_threshold,
            "high_state_threshold": self.high_state_threshold,
            "maximum_normalized_bracket_width": self.maximum_normalized_bracket_width,
            "replicate_count": len(self.replicates),
            "resolution_cells": [cell.as_dict() for cell in self.resolution_cells],
            "replicates": [replicate.as_dict() for replicate in self.replicates],
            "summary": dict(self.summary),
        }

    def to_csv_row(self) -> dict[str, object]:
        row: dict[str, object] = {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "pair_index": self.pair_index,
            "replicate_count": len(self.replicates),
            "endpoint_padding_fraction": self.endpoint_padding_fraction,
            "absolute_padding": self.absolute_padding,
            "stage_generations": self.stage_generations,
            "nested_barrier_points": ",".join(str(value) for value in self.nested_barrier_points),
            "interaction_separation_threshold": self.interaction_separation_threshold,
            "low_state_threshold": self.low_state_threshold,
            "high_state_threshold": self.high_state_threshold,
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


def run_finite_h1_boundary_resolution_audit(
    spec: ExperimentSpec,
    *,
    endpoint_padding_fraction: float = 0.5,
    stage_generations: int = 30,
    nested_barrier_points: Sequence[int] = DEFAULT_NESTED_BARRIER_POINTS,
    interaction_separation_threshold: float = 0.05,
    low_state_threshold: float = 0.25,
    high_state_threshold: float = 0.75,
    maximum_normalized_bracket_width: float = 0.03,
) -> tuple[FiniteH1BoundaryResolutionCell, ...]:
    """Refine finite collapse/recovery brackets on nested grids.

    A shared endpoint padding fraction is used for all canonical-H1-applicable
    pairs.  This is a confirmation design: it does not choose pair-specific
    endpoint ranges from their observed first-loop padding.
    """
    grids = _validate_nested_points(nested_barrier_points)
    if endpoint_padding_fraction <= 0.0:
        raise ValueError("endpoint_padding_fraction must be positive")
    if stage_generations < 1:
        raise ValueError("stage_generations must be positive")
    if interaction_separation_threshold < 0.0:
        raise ValueError("interaction_separation_threshold must be non-negative")
    if not 0.0 <= low_state_threshold < high_state_threshold <= 1.0:
        raise ValueError("state thresholds must satisfy 0 <= low < high <= 1")
    if maximum_normalized_bracket_width <= 0.0:
        raise ValueError("maximum_normalized_bracket_width must be positive")

    cells: list[FiniteH1BoundaryResolutionCell] = []
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
            absolute_padding = _absolute_padding(interval, endpoint_padding_fraction)
            resolution_cells = tuple(
                run_finite_h1_hysteresis_audit(
                    pair_spec,
                    scenarios=(scenario,),
                    barrier_points=points,
                    barrier_padding=absolute_padding,
                    stage_generations=stage_generations,
                    interaction_separation_threshold=interaction_separation_threshold,
                    low_state_threshold=low_state_threshold,
                    high_state_threshold=high_state_threshold,
                )[0]
                for points in grids
            )
            replicates = _build_replicates(
                resolution_cells,
                interval=interval,
                points=grids,
                low_state_threshold=low_state_threshold,
                high_state_threshold=high_state_threshold,
                maximum_normalized_bracket_width=maximum_normalized_bracket_width,
            )
            cells.append(
                FiniteH1BoundaryResolutionCell(
                    experiment_id=spec.experiment_id,
                    profile=spec.profile,
                    pair_index=pair_index,
                    parameters=resolution_cells[0].parameters,
                    canonical_bistable_barrier_interval=interval,
                    endpoint_padding_fraction=endpoint_padding_fraction,
                    absolute_padding=absolute_padding,
                    stage_generations=stage_generations,
                    nested_barrier_points=grids,
                    interaction_separation_threshold=interaction_separation_threshold,
                    low_state_threshold=low_state_threshold,
                    high_state_threshold=high_state_threshold,
                    maximum_normalized_bracket_width=maximum_normalized_bracket_width,
                    resolution_cells=resolution_cells,
                    replicates=replicates,
                    summary=_summarise(replicates, grids),
                )
            )
            pair_index += 1
    return tuple(cells)


def write_finite_h1_boundary_resolution_artifacts(
    cells: Iterable[FiniteH1BoundaryResolutionCell],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> None:
    """Write complete nested-grid audit records and flat cell summaries."""
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
    resolution_cells: Sequence[FiniteH1HysteresisCell],
    *,
    interval: tuple[float, float] | None,
    points: tuple[int, ...],
    low_state_threshold: float,
    high_state_threshold: float,
    maximum_normalized_bracket_width: float,
) -> tuple[FiniteH1BoundaryResolutionReplicate, ...]:
    by_resolution = [cell.replicates for cell in resolution_cells]
    count = len(by_resolution[0])
    if any(len(values) != count for values in by_resolution):
        raise RuntimeError("nested grids produced unequal replicate counts")
    output: list[FiniteH1BoundaryResolutionReplicate] = []
    for replicate_index in range(count):
        source = tuple(values[replicate_index] for values in by_resolution)
        if any(record.replicate_index != replicate_index for record in source):
            raise RuntimeError("nested-grid replicate order is inconsistent")
        if len({record.seed for record in source}) != 1:
            raise RuntimeError("same-seed nested-grid pairing was not preserved")
        observations = tuple(
            _observation(
                record,
                barrier_points=point,
                interval=interval,
                low_state_threshold=low_state_threshold,
                high_state_threshold=high_state_threshold,
            )
            for point, record in zip(points, source, strict=True)
        )
        if any(observation.finite_loop_bracket_supported is None for observation in observations):
            output.append(
                FiniteH1BoundaryResolutionReplicate(
                    replicate_index=replicate_index,
                    seed=source[0].seed,
                    observations=observations,
                    finest_rising_midpoint_shift=None,
                    finest_falling_midpoint_shift=None,
                    finest_grid_bracket_width_fraction=None,
                    loop_on_two_finest_grids=None,
                    boundary_location_stable=None,
                    resolution_stable_loop_supported=None,
                    resolution_stable_h1_loop_mechanism_supported=None,
                )
            )
            continue
        penultimate, finest = observations[-2:]
        loops_on_two_finest = bool(penultimate.finite_loop_bracket_supported) and bool(finest.finite_loop_bracket_supported)
        if not loops_on_two_finest:
            output.append(
                FiniteH1BoundaryResolutionReplicate(
                    replicate_index=replicate_index,
                    seed=source[0].seed,
                    observations=observations,
                    finest_rising_midpoint_shift=None,
                    finest_falling_midpoint_shift=None,
                    finest_grid_bracket_width_fraction=None,
                    loop_on_two_finest_grids=False,
                    boundary_location_stable=False,
                    resolution_stable_loop_supported=False,
                    resolution_stable_h1_loop_mechanism_supported=False,
                )
            )
            continue
        assert penultimate.rising_collapse_bracket is not None
        assert penultimate.falling_recovery_bracket is not None
        assert finest.rising_collapse_bracket is not None
        assert finest.falling_recovery_bracket is not None
        rising_shift = abs(finest.rising_collapse_bracket.midpoint - penultimate.rising_collapse_bracket.midpoint)
        falling_shift = abs(finest.falling_recovery_bracket.midpoint - penultimate.falling_recovery_bracket.midpoint)
        assert penultimate.barrier_step is not None
        interval_width = _interval_width(interval)
        finest_width_fraction = max(
            finest.rising_collapse_bracket.width,
            finest.falling_recovery_bracket.width,
        ) / interval_width
        location_stable = (
            rising_shift <= penultimate.barrier_step
            and falling_shift <= penultimate.barrier_step
            and finest_width_fraction <= maximum_normalized_bracket_width
        )
        output.append(
            FiniteH1BoundaryResolutionReplicate(
                replicate_index=replicate_index,
                seed=source[0].seed,
                observations=observations,
                finest_rising_midpoint_shift=rising_shift,
                finest_falling_midpoint_shift=falling_shift,
                finest_grid_bracket_width_fraction=finest_width_fraction,
                loop_on_two_finest_grids=True,
                boundary_location_stable=location_stable,
                resolution_stable_loop_supported=location_stable,
                resolution_stable_h1_loop_mechanism_supported=(
                    location_stable
                    and bool(penultimate.finite_h1_loop_bracket_mechanism_supported)
                    and bool(finest.finite_h1_loop_bracket_mechanism_supported)
                ),
            )
        )
    return tuple(output)


def _observation(
    replicate: FiniteH1HysteresisReplicate,
    *,
    barrier_points: int,
    interval: tuple[float, float] | None,
    low_state_threshold: float,
    high_state_threshold: float,
) -> FiniteH1BoundaryResolutionObservation:
    if replicate.finite_hysteresis_supported is None:
        return FiniteH1BoundaryResolutionObservation(
            barrier_points=barrier_points,
            barrier_step=None,
            finite_hysteresis_supported=None,
            finite_h1_hysteresis_mechanism_supported=None,
            rising_collapse_bracket=None,
            falling_recovery_bracket=None,
            finite_loop_bracket_supported=None,
            finite_h1_loop_bracket_mechanism_supported=None,
            finite_gap_lower_bound=None,
            finite_gap_upper_bound=None,
        )
    rising = _required(replicate.rising)
    falling = _required(replicate.falling)
    rising_bracket = _rising_bracket(rising, low_state_threshold)
    falling_bracket = _falling_bracket(falling, high_state_threshold)
    if interval is None:
        barrier_step = None
    else:
        lower = min(stage.barrier for stage in rising)
        upper = max(stage.barrier for stage in rising)
        barrier_step = (upper - lower) / (barrier_points - 1)
    if rising_bracket is None or falling_bracket is None:
        return FiniteH1BoundaryResolutionObservation(
            barrier_points=barrier_points,
            barrier_step=barrier_step,
            finite_hysteresis_supported=replicate.finite_hysteresis_supported,
            finite_h1_hysteresis_mechanism_supported=replicate.finite_h1_hysteresis_mechanism_supported,
            rising_collapse_bracket=rising_bracket,
            falling_recovery_bracket=falling_bracket,
            finite_loop_bracket_supported=False,
            finite_h1_loop_bracket_mechanism_supported=False,
            finite_gap_lower_bound=None,
            finite_gap_upper_bound=None,
        )
    gap_lower = rising_bracket.lower_barrier - falling_bracket.upper_barrier
    gap_upper = rising_bracket.upper_barrier - falling_bracket.lower_barrier
    loop_supported = gap_lower > 0.0
    return FiniteH1BoundaryResolutionObservation(
        barrier_points=barrier_points,
        barrier_step=barrier_step,
        finite_hysteresis_supported=replicate.finite_hysteresis_supported,
        finite_h1_hysteresis_mechanism_supported=replicate.finite_h1_hysteresis_mechanism_supported,
        rising_collapse_bracket=rising_bracket,
        falling_recovery_bracket=falling_bracket,
        finite_loop_bracket_supported=loop_supported,
        finite_h1_loop_bracket_mechanism_supported=(
            loop_supported and bool(replicate.finite_h1_hysteresis_mechanism_supported)
        ),
        finite_gap_lower_bound=gap_lower,
        finite_gap_upper_bound=gap_upper,
    )


def _rising_bracket(
    stages: Sequence[FiniteH1ContinuationStage],
    low_state_threshold: float,
) -> FiniteH1BoundaryBracket | None:
    for index, stage in enumerate(stages):
        if stage.terminal_interaction_mean <= low_state_threshold:
            if index == 0:
                return None
            prior = stages[index - 1]
            return FiniteH1BoundaryBracket(
                lower_barrier=prior.barrier,
                upper_barrier=stage.barrier,
                crossing_barrier=stage.barrier,
                prior_barrier=prior.barrier,
            )
    return None


def _falling_bracket(
    stages: Sequence[FiniteH1ContinuationStage],
    high_state_threshold: float,
) -> FiniteH1BoundaryBracket | None:
    for index, stage in enumerate(stages):
        if stage.terminal_interaction_mean >= high_state_threshold:
            if index == 0:
                return None
            prior = stages[index - 1]
            return FiniteH1BoundaryBracket(
                lower_barrier=stage.barrier,
                upper_barrier=prior.barrier,
                crossing_barrier=stage.barrier,
                prior_barrier=prior.barrier,
            )
    return None


def _summarise(
    replicates: Sequence[FiniteH1BoundaryResolutionReplicate],
    points: Sequence[int],
) -> dict[str, object]:
    if not replicates:
        raise ValueError("replicates must be nonempty")
    available = tuple(rep for rep in replicates if rep.resolution_stable_loop_supported is not None)
    summary: dict[str, object] = {
        "replicate_count": len(replicates),
        "available_nested_grid_probability": len(available) / len(replicates),
        "by_barrier_points": {},
    }
    for index, point in enumerate(points):
        observations = tuple(rep.observations[index] for rep in available)
        closed = tuple(obs for obs in observations if obs.finite_loop_bracket_supported)
        summary["by_barrier_points"][str(point)] = {
            "finite_hysteresis_supported_probability": None
            if not observations
            else _probability(bool(obs.finite_hysteresis_supported) for obs in observations),
            "finite_h1_hysteresis_mechanism_supported_probability": None
            if not observations
            else _probability(bool(obs.finite_h1_hysteresis_mechanism_supported) for obs in observations),
            "finite_loop_bracket_supported_probability": None
            if not observations
            else _probability(bool(obs.finite_loop_bracket_supported) for obs in observations),
            "finite_h1_loop_bracket_mechanism_supported_probability": None
            if not observations
            else _probability(bool(obs.finite_h1_loop_bracket_mechanism_supported) for obs in observations),
            "barrier_step": _summary(obs.barrier_step for obs in observations if obs.barrier_step is not None),
            "rising_collapse_midpoint": _summary(
                obs.rising_collapse_bracket.midpoint
                for obs in closed
                if obs.rising_collapse_bracket is not None
            ),
            "falling_recovery_midpoint": _summary(
                obs.falling_recovery_bracket.midpoint
                for obs in closed
                if obs.falling_recovery_bracket is not None
            ),
            "finite_gap_lower_bound": _summary(
                obs.finite_gap_lower_bound for obs in closed if obs.finite_gap_lower_bound is not None
            ),
            "finite_gap_upper_bound": _summary(
                obs.finite_gap_upper_bound for obs in closed if obs.finite_gap_upper_bound is not None
            ),
        }
    if not available:
        summary.update(
            {
                "loop_on_two_finest_grids_probability": None,
                "boundary_location_stable_probability": None,
                "resolution_stable_loop_supported_probability": None,
                "resolution_stable_h1_loop_mechanism_supported_probability": None,
                "finest_rising_midpoint_shift": _empty_summary(),
                "finest_falling_midpoint_shift": _empty_summary(),
                "finest_grid_bracket_width_fraction": _empty_summary(),
            }
        )
        return summary
    summary.update(
        {
            "loop_on_two_finest_grids_probability": _probability(
                bool(rep.loop_on_two_finest_grids) for rep in available
            ),
            "boundary_location_stable_probability": _probability(
                bool(rep.boundary_location_stable) for rep in available
            ),
            "resolution_stable_loop_supported_probability": _probability(
                bool(rep.resolution_stable_loop_supported) for rep in available
            ),
            "resolution_stable_h1_loop_mechanism_supported_probability": _probability(
                bool(rep.resolution_stable_h1_loop_mechanism_supported) for rep in available
            ),
            "finest_rising_midpoint_shift": _summary(
                _required(rep.finest_rising_midpoint_shift) for rep in available if rep.finest_rising_midpoint_shift is not None
            ),
            "finest_falling_midpoint_shift": _summary(
                _required(rep.finest_falling_midpoint_shift) for rep in available if rep.finest_falling_midpoint_shift is not None
            ),
            "finest_grid_bracket_width_fraction": _summary(
                _required(rep.finest_grid_bracket_width_fraction)
                for rep in available
                if rep.finest_grid_bracket_width_fraction is not None
            ),
        }
    )
    return summary


def _validate_nested_points(values: Sequence[int]) -> tuple[int, ...]:
    points = tuple(int(value) for value in values)
    if len(points) < 2:
        raise ValueError("nested_barrier_points must contain at least two grids")
    if any(value < 3 for value in points):
        raise ValueError("each barrier grid must have at least three points")
    if tuple(sorted(points)) != points or len(set(points)) != len(points):
        raise ValueError("nested barrier points must be strictly increasing")
    subdivisions = tuple(value - 1 for value in points)
    if any(later % earlier != 0 for earlier, later in zip(subdivisions, subdivisions[1:])):
        raise ValueError("nested barrier grids require each finer subdivision count to be divisible by the prior count")
    return points


def _absolute_padding(interval: tuple[float, float] | None, fraction: float) -> float:
    if interval is None:
        return 0.1
    return fraction * _interval_width(interval)


def _interval_width(interval: tuple[float, float] | None) -> float:
    if interval is None:
        raise RuntimeError("canonical interval is required for a finite boundary width")
    return interval[1] - interval[0]


def _pair_seed(master_seed: int, pair_index: int) -> int:
    return (master_seed * 1_000_003 + pair_index * 10_007 + 613) % (2**31 - 1)


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
        raise RuntimeError("unexpected missing finite H1 boundary-resolution value")
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
