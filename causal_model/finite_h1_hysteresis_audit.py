"""Finite continuation audit for H1 hysteresis.

Canonical H1 has an analytic hysteresis construction for the one-state map.
The finite coupled simulator has density variation, realised-trait and allele
feedback, finite recruitment, and optional migration.  This module performs a
separate Type S continuation experiment: it carries the terminal finite state
through rising and falling barrier paths and asks whether route histories remain
separated inside the canonical bistable interval.

The audit is deliberately not a theorem certificate.  It reports both support
and counterexamples, plus theorem-boundary audits for every continuation stage.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from statistics import median
from typing import Iterable, Mapping, Sequence

from causal_model.canonical_h1_bifurcation import (
    CanonicalH1Certificate,
    canonical_bistable_barrier_interval,
    canonical_h1_certificate,
)
from causal_model.h1_theorem_boundary_audit import H1TheoremBoundaryAudit, audit_h1_theorem_boundary
from causal_model.multipatch_criticality_dynamics import DynamicsParameters, SimulationResult, SimulationSnapshot, simulate
from causal_model.multipatch_criticality_experiments import (
    ExperimentSpec,
    LandscapeScenario,
    ParameterCell,
    default_scenarios,
    derived_seed,
    parameter_grid,
    parameters_for_cell,
)


@dataclass(frozen=True)
class FiniteH1ContinuationStage:
    """Terminal state from one carried-over finite continuation stage."""

    stage_index: int
    barrier: float
    stage_seed: int
    terminal_interaction_mean: float
    terminal_potential_high_trait_viable: bool
    terminal_high_trait_mass_mean: float
    terminal_local_effective_size_mean: float
    h1_scope: H1TheoremBoundaryAudit

    def as_dict(self) -> dict[str, object]:
        return {
            "stage_index": self.stage_index,
            "barrier": self.barrier,
            "stage_seed": self.stage_seed,
            "terminal_interaction_mean": self.terminal_interaction_mean,
            "terminal_potential_high_trait_viable": self.terminal_potential_high_trait_viable,
            "terminal_high_trait_mass_mean": self.terminal_high_trait_mass_mean,
            "terminal_local_effective_size_mean": self.terminal_local_effective_size_mean,
            "h1_scope": self.h1_scope.as_dict(),
        }


@dataclass(frozen=True)
class FiniteH1HysteresisReplicate:
    """Rising and falling finite continuation paths for one parameter replicate."""

    replicate_index: int
    seed: int
    rising: tuple[FiniteH1ContinuationStage, ...] | None
    falling: tuple[FiniteH1ContinuationStage, ...] | None
    internal_barrier_gaps: tuple[tuple[float, float], ...] | None
    internal_potential_switch_barriers: tuple[float, ...] | None
    maximum_internal_interaction_gap: float | None
    rising_collapse_barrier: float | None
    falling_recovery_barrier: float | None
    jump_boundary_gap: float | None
    finite_hysteresis_supported: bool | None
    finite_h1_hysteresis_mechanism_supported: bool | None

    def as_dict(self) -> dict[str, object]:
        return {
            "replicate_index": self.replicate_index,
            "seed": self.seed,
            "rising": None if self.rising is None else [stage.as_dict() for stage in self.rising],
            "falling": None if self.falling is None else [stage.as_dict() for stage in self.falling],
            "internal_barrier_gaps": None
            if self.internal_barrier_gaps is None
            else [{"barrier": barrier, "rising_minus_falling_interaction": gap} for barrier, gap in self.internal_barrier_gaps],
            "internal_potential_switch_barriers": None
            if self.internal_potential_switch_barriers is None
            else list(self.internal_potential_switch_barriers),
            "maximum_internal_interaction_gap": self.maximum_internal_interaction_gap,
            "rising_collapse_barrier": self.rising_collapse_barrier,
            "falling_recovery_barrier": self.falling_recovery_barrier,
            "jump_boundary_gap": self.jump_boundary_gap,
            "finite_hysteresis_supported": self.finite_hysteresis_supported,
            "finite_h1_hysteresis_mechanism_supported": self.finite_h1_hysteresis_mechanism_supported,
        }


@dataclass(frozen=True)
class FiniteH1HysteresisCell:
    """Finite hysteresis evidence for one landscape and canonical parameter pair."""

    experiment_id: str
    profile: str
    scenario_id: str
    parameters: ParameterCell
    canonical_midpoint_h1: CanonicalH1Certificate
    canonical_bistable_barrier_interval: tuple[float, float] | None
    barriers: tuple[float, ...] | None
    stage_generations: int
    interaction_separation_threshold: float
    low_state_threshold: float
    high_state_threshold: float
    replicates: tuple[FiniteH1HysteresisReplicate, ...]
    summary: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "scenario_id": self.scenario_id,
            "parameters": asdict(self.parameters),
            "canonical_midpoint_h1": asdict(self.canonical_midpoint_h1),
            "canonical_bistable_barrier_interval": None
            if self.canonical_bistable_barrier_interval is None
            else list(self.canonical_bistable_barrier_interval),
            "barriers": None if self.barriers is None else list(self.barriers),
            "stage_generations": self.stage_generations,
            "interaction_separation_threshold": self.interaction_separation_threshold,
            "low_state_threshold": self.low_state_threshold,
            "high_state_threshold": self.high_state_threshold,
            "replicate_count": len(self.replicates),
            "replicates": [replicate.as_dict() for replicate in self.replicates],
            "summary": dict(self.summary),
        }

    def to_csv_row(self) -> dict[str, object]:
        row: dict[str, object] = {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "scenario_id": self.scenario_id,
            "replicate_count": len(self.replicates),
            "stage_generations": self.stage_generations,
            "interaction_separation_threshold": self.interaction_separation_threshold,
            "low_state_threshold": self.low_state_threshold,
            "high_state_threshold": self.high_state_threshold,
        }
        row.update(asdict(self.parameters))
        row.update(_flatten_mapping(self.summary))
        return row


def run_finite_h1_hysteresis_audit(
    spec: ExperimentSpec,
    scenarios: Sequence[LandscapeScenario] | None = None,
    *,
    barrier_points: int = 7,
    barrier_padding: float = 0.1,
    stage_generations: int = 10,
    interaction_separation_threshold: float = 0.1,
    low_state_threshold: float = 0.25,
    high_state_threshold: float = 0.75,
) -> tuple[FiniteH1HysteresisCell, ...]:
    """Carry finite states through rising and falling barriers.

    The audit is applicable only where the local canonical map has a strict
    bistable interval and its midpoint has a branch-dependent high-trait mode.
    For applicable cells, rising continuation starts above the canonical high
    branch at a low barrier; falling continuation starts below the canonical low
    branch at a high barrier.  At each stage, the complete terminal finite state
    becomes the next stage's initial state.

    `finite_hysteresis_supported` is a predeclared route-memory predicate: at
    least one shared barrier strictly inside the canonical interval has rising
    interaction minus falling interaction larger than the declared threshold.
    The stronger mechanism predicate additionally requires a potential
    high-trait switch at one such internally shared barrier.
    """
    if barrier_points < 3:
        raise ValueError("barrier_points must be at least three")
    if barrier_padding <= 0.0:
        raise ValueError("barrier_padding must be positive")
    if stage_generations < 1:
        raise ValueError("stage_generations must be positive")
    if interaction_separation_threshold < 0.0:
        raise ValueError("interaction_separation_threshold must be non-negative")
    if not 0.0 <= low_state_threshold < high_state_threshold <= 1.0:
        raise ValueError("state thresholds must satisfy 0 <= low < high <= 1")

    selected = tuple(default_scenarios(spec) if scenarios is None else scenarios)
    if not selected:
        raise ValueError("scenarios must be nonempty")
    if len({scenario.scenario_id for scenario in selected}) != len(selected):
        raise ValueError("scenario identifiers must be unique")
    if any(len(set(scenario.patch_areas)) > 1 for scenario in selected):
        raise ValueError("finite H1 hysteresis audit requires uniform local patch areas")

    cells: list[FiniteH1HysteresisCell] = []
    for scenario in selected:
        local_area = scenario.patch_areas[0]
        for cell in parameter_grid(spec):
            base = parameters_for_cell(spec, scenario, cell, seed=spec.master_seed)
            interval = canonical_bistable_barrier_interval(
                cell.interaction_feedback,
                local_area,
                cell.area_reference,
            )
            midpoint = cell.interaction_barrier if interval is None else (interval[0] + interval[1]) / 2.0
            canonical = canonical_h1_certificate(
                feedback_strength=cell.interaction_feedback,
                area=local_area,
                area_reference=cell.area_reference,
                barrier=midpoint,
                trait_parameters=base,
            )
            barriers = None
            if interval is not None and canonical.branch_dependent_high_trait_mode:
                barriers = _barrier_grid(interval, barrier_points=barrier_points, padding=barrier_padding)
            pairs = tuple(
                _run_hysteresis_replicate(
                    spec,
                    scenario,
                    cell,
                    interval=interval,
                    canonical=canonical,
                    barriers=barriers,
                    replicate_index=index,
                    stage_generations=stage_generations,
                    interaction_separation_threshold=interaction_separation_threshold,
                    low_state_threshold=low_state_threshold,
                    high_state_threshold=high_state_threshold,
                )
                for index in range(spec.replicates)
            )
            cells.append(
                FiniteH1HysteresisCell(
                    experiment_id=spec.experiment_id,
                    profile=spec.profile,
                    scenario_id=scenario.scenario_id,
                    parameters=cell,
                    canonical_midpoint_h1=canonical,
                    canonical_bistable_barrier_interval=interval,
                    barriers=barriers,
                    stage_generations=stage_generations,
                    interaction_separation_threshold=interaction_separation_threshold,
                    low_state_threshold=low_state_threshold,
                    high_state_threshold=high_state_threshold,
                    replicates=pairs,
                    summary=_summarise_hysteresis_pairs(pairs, canonical, interval),
                )
            )
    return tuple(cells)


def write_finite_h1_hysteresis_artifacts(
    cells: Iterable[FiniteH1HysteresisCell],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> None:
    """Write flat summaries and complete finite continuation trajectories."""
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


def _barrier_grid(interval: tuple[float, float], *, barrier_points: int, padding: float) -> tuple[float, ...]:
    lower, upper = interval
    start = lower - padding
    stop = upper + padding
    step = (stop - start) / (barrier_points - 1)
    return tuple(start + index * step for index in range(barrier_points))


def _run_hysteresis_replicate(
    spec: ExperimentSpec,
    scenario: LandscapeScenario,
    cell: ParameterCell,
    *,
    interval: tuple[float, float] | None,
    canonical: CanonicalH1Certificate,
    barriers: tuple[float, ...] | None,
    replicate_index: int,
    stage_generations: int,
    interaction_separation_threshold: float,
    low_state_threshold: float,
    high_state_threshold: float,
) -> FiniteH1HysteresisReplicate:
    seed = derived_seed(spec.master_seed, cell.cell_index, replicate_index)
    if interval is None or barriers is None or not canonical.branch_dependent_high_trait_mode:
        return FiniteH1HysteresisReplicate(
            replicate_index=replicate_index,
            seed=seed,
            rising=None,
            falling=None,
            internal_barrier_gaps=None,
            internal_potential_switch_barriers=None,
            maximum_internal_interaction_gap=None,
            rising_collapse_barrier=None,
            falling_recovery_barrier=None,
            jump_boundary_gap=None,
            finite_hysteresis_supported=None,
            finite_h1_hysteresis_mechanism_supported=None,
        )

    assert canonical.low_stable_branch is not None
    assert canonical.high_stable_branch is not None
    base = parameters_for_cell(spec, scenario, cell, seed=seed)
    rising_initial = _with_uniform_initial_interaction(base, canonical.high_stable_branch.interaction)
    falling_initial = _with_uniform_initial_interaction(base, canonical.low_stable_branch.interaction)
    rising = _run_route(
        rising_initial,
        barriers,
        route_code=1,
        replicate_seed=seed,
        stage_generations=stage_generations,
    )
    falling = _run_route(
        falling_initial,
        tuple(reversed(barriers)),
        route_code=2,
        replicate_seed=seed,
        stage_generations=stage_generations,
    )
    falling_by_barrier = {stage.barrier: stage for stage in falling}
    internal_gaps: list[tuple[float, float]] = []
    potential_switches: list[float] = []
    lower, upper = interval
    for stage in rising:
        if lower < stage.barrier < upper:
            reverse_stage = falling_by_barrier[stage.barrier]
            gap = stage.terminal_interaction_mean - reverse_stage.terminal_interaction_mean
            internal_gaps.append((stage.barrier, gap))
            if stage.terminal_potential_high_trait_viable and not reverse_stage.terminal_potential_high_trait_viable:
                potential_switches.append(stage.barrier)
    maximum_gap = max((gap for _, gap in internal_gaps), default=0.0)
    rising_collapse = _first_barrier(rising, lambda stage: stage.terminal_interaction_mean <= low_state_threshold)
    falling_recovery = _first_barrier(falling, lambda stage: stage.terminal_interaction_mean >= high_state_threshold)
    jump_gap = None if rising_collapse is None or falling_recovery is None else rising_collapse - falling_recovery
    hysteresis = maximum_gap > interaction_separation_threshold
    mechanism = hysteresis and bool(potential_switches)
    return FiniteH1HysteresisReplicate(
        replicate_index=replicate_index,
        seed=seed,
        rising=rising,
        falling=falling,
        internal_barrier_gaps=tuple(internal_gaps),
        internal_potential_switch_barriers=tuple(potential_switches),
        maximum_internal_interaction_gap=maximum_gap,
        rising_collapse_barrier=rising_collapse,
        falling_recovery_barrier=falling_recovery,
        jump_boundary_gap=jump_gap,
        finite_hysteresis_supported=hysteresis,
        finite_h1_hysteresis_mechanism_supported=mechanism,
    )


def _run_route(
    initial: DynamicsParameters,
    barriers: Sequence[float],
    *,
    route_code: int,
    replicate_seed: int,
    stage_generations: int,
) -> tuple[FiniteH1ContinuationStage, ...]:
    parameters = initial
    stages: list[FiniteH1ContinuationStage] = []
    for index, barrier in enumerate(barriers):
        stage_seed = _stage_seed(replicate_seed, route_code, index)
        result = simulate(
            replace(
                parameters,
                interaction_barrier=barrier,
                generations=stage_generations,
                random_seed=stage_seed,
            )
        )
        terminal = result.snapshots[-1]
        stages.append(_stage_from_result(index, barrier, stage_seed, result))
        parameters = _parameters_from_terminal(parameters, terminal)
    return tuple(stages)


def _stage_seed(replicate_seed: int, route_code: int, stage_index: int) -> int:
    return (replicate_seed * 1_000_003 + route_code * 10_007 + stage_index * 101 + 29) % (2**31 - 1)


def _with_uniform_initial_interaction(parameters: DynamicsParameters, interaction: float) -> DynamicsParameters:
    return replace(parameters, initial_interaction=tuple(interaction for _ in parameters.patch_areas))


def _parameters_from_terminal(parameters: DynamicsParameters, terminal: SimulationSnapshot) -> DynamicsParameters:
    return replace(
        parameters,
        initial_population=terminal.population,
        initial_interaction=terminal.interaction,
        initial_high_allele_frequency=terminal.high_allele_frequency,
        initial_trait_distribution=(),
        initial_trait_abundance=tuple(summary.abundance for summary in terminal.trait_occupancy),
    )


def _stage_from_result(
    index: int,
    barrier: float,
    stage_seed: int,
    result: SimulationResult,
) -> FiniteH1ContinuationStage:
    terminal = result.snapshots[-1]
    occupancy = terminal.trait_occupancy
    return FiniteH1ContinuationStage(
        stage_index=index,
        barrier=barrier,
        stage_seed=stage_seed,
        terminal_interaction_mean=_mean(terminal.interaction),
        terminal_potential_high_trait_viable=any(
            summary.high_trait_component_present for summary in terminal.trait_space
        ),
        terminal_high_trait_mass_mean=_mean(tuple(summary.high_trait_mass for summary in occupancy)),
        terminal_local_effective_size_mean=_mean(terminal.effective_size),
        h1_scope=audit_h1_theorem_boundary(result),
    )


def _first_barrier(
    stages: Sequence[FiniteH1ContinuationStage],
    predicate,
) -> float | None:
    return next((stage.barrier for stage in stages if predicate(stage)), None)


def _summarise_hysteresis_pairs(
    pairs: Sequence[FiniteH1HysteresisReplicate],
    canonical: CanonicalH1Certificate,
    interval: tuple[float, float] | None,
) -> dict[str, object]:
    available = tuple(pair for pair in pairs if pair.finite_hysteresis_supported is not None)
    summary: dict[str, object] = {
        "canonical_context": _canonical_context(canonical, interval),
        "finite_continuation_available_probability": len(available) / len(pairs),
    }
    if not available:
        summary.update(
            {
                "finite_hysteresis_supported_probability": None,
                "finite_h1_hysteresis_mechanism_supported_probability": None,
                "maximum_internal_interaction_gap": _empty_summary(),
                "jump_boundary_gap": _empty_summary(),
                "internal_potential_switch_probability": None,
                "h1_theorem_scope": {"rising": None, "falling": None},
            }
        )
        return summary
    all_rising = tuple(stage for pair in available for stage in _required(pair.rising))
    all_falling = tuple(stage for pair in available for stage in _required(pair.falling))
    summary.update(
        {
            "finite_hysteresis_supported_probability": _probability(
                bool(pair.finite_hysteresis_supported) for pair in available
            ),
            "finite_h1_hysteresis_mechanism_supported_probability": _probability(
                bool(pair.finite_h1_hysteresis_mechanism_supported) for pair in available
            ),
            "maximum_internal_interaction_gap": _summary(
                _required(pair.maximum_internal_interaction_gap) for pair in available
            ),
            "jump_boundary_gap": _summary(
                pair.jump_boundary_gap for pair in available if pair.jump_boundary_gap is not None
            ),
            "internal_potential_switch_probability": _probability(
                bool(_required(pair.internal_potential_switch_barriers)) for pair in available
            ),
            "h1_theorem_scope": {
                "rising": _summarise_scope(stage.h1_scope for stage in all_rising),
                "falling": _summarise_scope(stage.h1_scope for stage in all_falling),
            },
        }
    )
    return summary


def _canonical_context(
    canonical: CanonicalH1Certificate,
    interval: tuple[float, float] | None,
) -> dict[str, object]:
    return {
        "gain": canonical.bifurcation.gain,
        "strict_bistability_certified_at_midpoint": canonical.bifurcation.strict_bistability_certified,
        "branch_dependent_high_trait_mode_at_midpoint": canonical.branch_dependent_high_trait_mode,
        "bistable_interval_lower": None if interval is None else interval[0],
        "bistable_interval_upper": None if interval is None else interval[1],
        "low_branch_interaction": None if canonical.low_stable_branch is None else canonical.low_stable_branch.interaction,
        "high_branch_interaction": None if canonical.high_stable_branch is None else canonical.high_stable_branch.interaction,
    }


def _summarise_scope(audits: Iterable[H1TheoremBoundaryAudit]) -> dict[str, object]:
    values = tuple(audits)
    return {
        "patchwise_canonical_update_probability": _probability(
            audit.patchwise_canonical_update_certified for audit in values
        ),
        "single_patch_canonical_theorem_limit_probability": _probability(
            audit.single_patch_canonical_theorem_limit_certified for audit in values
        ),
        "maximum_canonical_update_residual": _summary(
            audit.maximum_canonical_update_residual for audit in values
        ),
    }


def _mean(values: Sequence[float]) -> float:
    if not values:
        raise ValueError("values must be nonempty")
    return sum(float(value) for value in values) / len(values)


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
        raise RuntimeError("unexpected missing finite continuation value")
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
