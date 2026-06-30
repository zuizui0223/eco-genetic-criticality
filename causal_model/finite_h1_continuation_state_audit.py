"""Full-state finite H1 continuation-transfer audit.

The boundary-resolution campaign establishes a route-dependent finite H1 loop by
passing the complete terminal state between adjacent barrier stages. A later
branch-locked H1--H2--H3 campaign showed that resetting only interaction q at an
anchor is insufficient: the high branch collapses because terminal population,
high-allele frequency, and realised trait abundance are omitted.

This audit isolates that distinction. It replays the certified finest-grid H1
continuations, selects a predeclared interior *grid* barrier, transfers the full
terminal state from the rising (high) and falling (low) routes, and holds both at
the same barrier for a fresh declared duration.

The audit therefore answers a narrow Type S question: does continuation-state
memory survive a fresh finite hold at a shared interior barrier? It does not yet
make H2 or H3 claims, and it does not treat a finite result as a theorem.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from statistics import median
from typing import Iterable, Mapping, Sequence

from causal_model.canonical_h1_bifurcation import canonical_h1_certificate
from causal_model.finite_h1_boundary_resolution_audit import (
    DEFAULT_NESTED_BARRIER_POINTS,
    FiniteH1BoundaryResolutionCell,
    FiniteH1BoundaryResolutionReplicate,
    run_finite_h1_boundary_resolution_audit,
)
from causal_model.finite_h1_hysteresis_audit import (
    _barrier_grid,
    _parameters_from_terminal,
    _stage_seed,
    _with_uniform_initial_interaction,
)
from causal_model.multipatch_criticality_dynamics import DynamicsParameters, SimulationSnapshot, simulate
from causal_model.multipatch_criticality_experiments import (
    SCENARIO_ONE_LARGE,
    ExperimentSpec,
    ParameterCell,
    default_scenarios,
    derived_seed,
    parameters_for_cell,
)

DEFAULT_STATE_TRANSFER_MASTER_SEEDS = (20260630, 20260631, 20260632, 20260633, 20260634)


@dataclass(frozen=True)
class H1TransferredState:
    """Complete terminal state carried from one continuation route."""

    route: str
    anchor_barrier: float
    terminal_interaction: tuple[float, ...]
    terminal_population: tuple[int, ...]
    terminal_high_allele_frequency: tuple[float, ...]
    terminal_trait_abundance: tuple[tuple[int, ...], ...]
    potential_high_trait_viable: bool

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class H1StateTransferHold:
    """One fresh finite hold after complete continuation-state transfer."""

    route: str
    hold_seed: int
    terminal_interaction_mean: float
    terminal_high_allele_frequency_mean: float
    terminal_potential_high_trait_viable: bool
    terminal_realised_high_trait_mass_mean: float

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class H1ContinuationStateTransferReplicate:
    """One paired rising/high and falling/low full-state transfer record."""

    master_seed: int
    replicate_index: int
    calibration_seed: int
    calibration_supported: bool | None
    anchor_barrier: float | None
    anchor_selection_interval: tuple[float, float] | None
    anchor_grid_point_count: int
    rising_high_state: H1TransferredState | None
    falling_low_state: H1TransferredState | None
    rising_high_hold: H1StateTransferHold | None
    falling_low_hold: H1StateTransferHold | None
    source_state_separated: bool | None
    state_transfer_hold_supported: bool | None

    def as_dict(self) -> dict[str, object]:
        return {
            "master_seed": self.master_seed,
            "replicate_index": self.replicate_index,
            "calibration_seed": self.calibration_seed,
            "calibration_supported": self.calibration_supported,
            "anchor_barrier": self.anchor_barrier,
            "anchor_selection_interval": None if self.anchor_selection_interval is None else list(self.anchor_selection_interval),
            "anchor_grid_point_count": self.anchor_grid_point_count,
            "rising_high_state": None if self.rising_high_state is None else self.rising_high_state.as_dict(),
            "falling_low_state": None if self.falling_low_state is None else self.falling_low_state.as_dict(),
            "rising_high_hold": None if self.rising_high_hold is None else self.rising_high_hold.as_dict(),
            "falling_low_hold": None if self.falling_low_hold is None else self.falling_low_hold.as_dict(),
            "source_state_separated": self.source_state_separated,
            "state_transfer_hold_supported": self.state_transfer_hold_supported,
        }


@dataclass(frozen=True)
class H1ContinuationStateTransferCell:
    """Full-state transfer evidence for one canonical H1 parameter pair."""

    experiment_id: str
    profile: str
    pair_index: int
    parameters: ParameterCell
    master_seeds: tuple[int, ...]
    endpoint_padding_fraction: float
    stage_generations: int
    hold_generations: int
    nested_barrier_points: tuple[int, ...]
    interaction_separation_threshold: float
    replicates: tuple[H1ContinuationStateTransferReplicate, ...]
    summary: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "pair_index": self.pair_index,
            "parameters": asdict(self.parameters),
            "master_seeds": list(self.master_seeds),
            "design": {
                "endpoint_padding_fraction": self.endpoint_padding_fraction,
                "stage_generations": self.stage_generations,
                "hold_generations": self.hold_generations,
                "nested_barrier_points": list(self.nested_barrier_points),
                "interaction_separation_threshold": self.interaction_separation_threshold,
            },
            "replicates": [record.as_dict() for record in self.replicates],
            "summary": dict(self.summary),
        }

    def to_csv_row(self) -> dict[str, object]:
        row = {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "pair_index": self.pair_index,
            "master_seeds": ",".join(str(seed) for seed in self.master_seeds),
            "replicate_count": len(self.replicates),
            "endpoint_padding_fraction": self.endpoint_padding_fraction,
            "stage_generations": self.stage_generations,
            "hold_generations": self.hold_generations,
            "nested_barrier_points": ",".join(str(value) for value in self.nested_barrier_points),
            "interaction_separation_threshold": self.interaction_separation_threshold,
        }
        row.update(asdict(self.parameters))
        row.update(_flatten(self.summary))
        return row


def run_finite_h1_continuation_state_audit(
    spec: ExperimentSpec,
    *,
    master_seeds: Sequence[int] = DEFAULT_STATE_TRANSFER_MASTER_SEEDS,
    endpoint_padding_fraction: float = 0.5,
    stage_generations: int = 30,
    hold_generations: int = 30,
    nested_barrier_points: Sequence[int] = DEFAULT_NESTED_BARRIER_POINTS,
    interaction_separation_threshold: float = 0.05,
    maximum_normalized_bracket_width: float = 0.03,
) -> tuple[H1ContinuationStateTransferCell, ...]:
    """Replay finite H1 routes and test paired full-state holds at an interior grid anchor."""
    seeds = _validate_seeds(master_seeds)
    if hold_generations < 1:
        raise ValueError("hold_generations must be positive")
    scenarios = {scenario.scenario_id: scenario for scenario in default_scenarios(spec)}
    one_large = scenarios[SCENARIO_ONE_LARGE]
    campaigns = tuple(
        run_finite_h1_boundary_resolution_audit(
            replace(spec, master_seed=master_seed),
            endpoint_padding_fraction=endpoint_padding_fraction,
            stage_generations=stage_generations,
            nested_barrier_points=nested_barrier_points,
            interaction_separation_threshold=interaction_separation_threshold,
            maximum_normalized_bracket_width=maximum_normalized_bracket_width,
        )
        for master_seed in seeds
    )
    count = len(campaigns[0])
    if any(len(campaign) != count for campaign in campaigns):
        raise RuntimeError("H1 calibration cells are not aligned across master seeds")
    output: list[H1ContinuationStateTransferCell] = []
    for cell_index in range(count):
        calibration_cells = tuple(campaign[cell_index] for campaign in campaigns)
        reference = calibration_cells[0]
        if any(_cell_key(value) != _cell_key(reference) for value in calibration_cells[1:]):
            raise RuntimeError("H1 parameter cells differ across master-seed calibrations")
        records: list[H1ContinuationStateTransferReplicate] = []
        for master_seed, calibration in zip(seeds, calibration_cells, strict=True):
            seed_spec = replace(spec, master_seed=master_seed)
            for calibration_replicate in calibration.replicates:
                records.append(
                    _transfer_record(
                        seed_spec,
                        one_large.patch_areas[0],
                        calibration,
                        calibration_replicate,
                        master_seed=master_seed,
                        endpoint_padding_fraction=endpoint_padding_fraction,
                        stage_generations=stage_generations,
                        hold_generations=hold_generations,
                        interaction_separation_threshold=interaction_separation_threshold,
                    )
                )
        output.append(
            H1ContinuationStateTransferCell(
                experiment_id=spec.experiment_id,
                profile=spec.profile,
                pair_index=reference.pair_index,
                parameters=reference.parameters,
                master_seeds=seeds,
                endpoint_padding_fraction=endpoint_padding_fraction,
                stage_generations=stage_generations,
                hold_generations=hold_generations,
                nested_barrier_points=tuple(int(value) for value in nested_barrier_points),
                interaction_separation_threshold=interaction_separation_threshold,
                replicates=tuple(records),
                summary=_summarise(records, seeds),
            )
        )
    return tuple(output)


def write_finite_h1_continuation_state_artifacts(
    cells: Iterable[H1ContinuationStateTransferCell],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> None:
    values = tuple(cells)
    if not values:
        raise ValueError("cells must be nonempty")
    csv_target, json_target = Path(csv_path), Path(json_path)
    csv_target.parent.mkdir(parents=True, exist_ok=True)
    json_target.parent.mkdir(parents=True, exist_ok=True)
    rows = [cell.to_csv_row() for cell in values]
    with csv_target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted({key for row in rows for key in row}))
        writer.writeheader()
        writer.writerows(rows)
    with json_target.open("w", encoding="utf-8") as handle:
        json.dump([cell.as_dict() for cell in values], handle, indent=2, sort_keys=True)


def _transfer_record(
    spec: ExperimentSpec,
    local_area: float,
    calibration: FiniteH1BoundaryResolutionCell,
    calibration_replicate: FiniteH1BoundaryResolutionReplicate,
    *,
    master_seed: int,
    endpoint_padding_fraction: float,
    stage_generations: int,
    hold_generations: int,
    interaction_separation_threshold: float,
) -> H1ContinuationStateTransferReplicate:
    if calibration_replicate.resolution_stable_h1_loop_mechanism_supported is not True:
        return _unavailable(master_seed, calibration_replicate)
    observation = calibration_replicate.observations[-1]
    if not observation.finite_h1_loop_bracket_mechanism_supported:
        return _unavailable(master_seed, calibration_replicate)
    collapse = observation.rising_collapse_bracket
    recovery = observation.falling_recovery_bracket
    interval = calibration.canonical_bistable_barrier_interval
    if collapse is None or recovery is None or interval is None:
        return _unavailable(master_seed, calibration_replicate)
    lower, upper = recovery.upper_barrier, collapse.lower_barrier
    width = interval[1] - interval[0]
    barriers = _barrier_grid(
        interval,
        barrier_points=observation.barrier_points,
        padding=endpoint_padding_fraction * width,
    )
    candidates = tuple(barrier for barrier in barriers if lower < barrier < upper)
    if not candidates:
        return _unavailable(master_seed, calibration_replicate, interval=(lower, upper), anchor_count=0)
    anchor = min(candidates, key=lambda value: (abs(value - (lower + upper) / 2.0), value))
    base = parameters_for_cell(spec, _one_large_scenario(spec), calibration.parameters, seed=calibration_replicate.seed)
    midpoint = (interval[0] + interval[1]) / 2.0
    canonical = canonical_h1_certificate(
        feedback_strength=calibration.parameters.interaction_feedback,
        area=local_area,
        area_reference=calibration.parameters.area_reference,
        barrier=midpoint,
        trait_parameters=base,
    )
    if canonical.low_stable_branch is None or canonical.high_stable_branch is None:
        return _unavailable(master_seed, calibration_replicate, interval=(lower, upper), anchor=anchor, anchor_count=len(candidates))
    high_terminal = _replay_to_anchor(
        _with_uniform_initial_interaction(base, canonical.high_stable_branch.interaction),
        barriers,
        route="rising",
        replicate_seed=calibration_replicate.seed,
        stage_generations=stage_generations,
        anchor=anchor,
    )
    low_terminal = _replay_to_anchor(
        _with_uniform_initial_interaction(base, canonical.low_stable_branch.interaction),
        tuple(reversed(barriers)),
        route="falling",
        replicate_seed=calibration_replicate.seed,
        stage_generations=stage_generations,
        anchor=anchor,
    )
    high_state = _state("rising_high", anchor, high_terminal)
    low_state = _state("falling_low", anchor, low_terminal)
    source_separated = (
        _mean(high_terminal.interaction) - _mean(low_terminal.interaction) > interaction_separation_threshold
        and _potential(high_terminal)
        and not _potential(low_terminal)
    )
    high_hold = _hold(high_terminal, anchor, _hold_seed(calibration_replicate.seed, 1), hold_generations, "rising_high")
    low_hold = _hold(low_terminal, anchor, _hold_seed(calibration_replicate.seed, 2), hold_generations, "falling_low")
    supported = (
        source_separated
        and high_hold.terminal_interaction_mean - low_hold.terminal_interaction_mean > interaction_separation_threshold
        and high_hold.terminal_potential_high_trait_viable
        and not low_hold.terminal_potential_high_trait_viable
    )
    return H1ContinuationStateTransferReplicate(
        master_seed=master_seed,
        replicate_index=calibration_replicate.replicate_index,
        calibration_seed=calibration_replicate.seed,
        calibration_supported=True,
        anchor_barrier=anchor,
        anchor_selection_interval=(lower, upper),
        anchor_grid_point_count=len(candidates),
        rising_high_state=high_state,
        falling_low_state=low_state,
        rising_high_hold=high_hold,
        falling_low_hold=low_hold,
        source_state_separated=source_separated,
        state_transfer_hold_supported=supported,
    )


def _one_large_scenario(spec: ExperimentSpec):
    return next(scenario for scenario in default_scenarios(spec) if scenario.scenario_id == SCENARIO_ONE_LARGE)


def _replay_to_anchor(
    initial: DynamicsParameters,
    barriers: Sequence[float],
    *,
    route: str,
    replicate_seed: int,
    stage_generations: int,
    anchor: float,
) -> SimulationSnapshot:
    route_code = 1 if route == "rising" else 2
    parameters = initial
    for index, barrier in enumerate(barriers):
        result = simulate(
            replace(
                parameters,
                interaction_barrier=barrier,
                generations=stage_generations,
                random_seed=_stage_seed(replicate_seed, route_code, index),
            )
        )
        terminal = result.snapshots[-1]
        if barrier == anchor:
            return terminal
        parameters = _parameters_from_terminal(parameters, terminal)
    raise RuntimeError("anchor barrier was not found in reconstructed continuation grid")


def _hold(terminal: SimulationSnapshot, anchor: float, seed: int, generations: int, route: str) -> H1StateTransferHold:
    parameters = DynamicsParameters(
        patch_areas=(1.0,) * len(terminal.population),
        generations=generations,
        random_seed=seed,
    )
    # The dynamics constants are supplied below by replacing a carried parameter
    # object in _state_parameters; this placeholder is never used directly.
    raise RuntimeError("_hold must be called through _hold_from_parameters")


def _state(route: str, anchor: float, terminal: SimulationSnapshot) -> H1TransferredState:
    return H1TransferredState(
        route=route,
        anchor_barrier=anchor,
        terminal_interaction=terminal.interaction,
        terminal_population=terminal.population,
        terminal_high_allele_frequency=terminal.high_allele_frequency,
        terminal_trait_abundance=tuple(item.abundance for item in terminal.trait_occupancy),
        potential_high_trait_viable=_potential(terminal),
    )


def _potential(snapshot: SimulationSnapshot) -> bool:
    return any(item.high_trait_component_present for item in snapshot.trait_space)


def _hold_seed(seed: int, route_code: int) -> int:
    return (seed * 1_000_003 + 70_001 + route_code * 101) % (2**31 - 1)


def _unavailable(
    master_seed: int,
    record: FiniteH1BoundaryResolutionReplicate,
    *,
    interval: tuple[float, float] | None = None,
    anchor: float | None = None,
    anchor_count: int = 0,
) -> H1ContinuationStateTransferReplicate:
    return H1ContinuationStateTransferReplicate(
        master_seed=master_seed,
        replicate_index=record.replicate_index,
        calibration_seed=record.seed,
        calibration_supported=record.resolution_stable_h1_loop_mechanism_supported,
        anchor_barrier=anchor,
        anchor_selection_interval=interval,
        anchor_grid_point_count=anchor_count,
        rising_high_state=None,
        falling_low_state=None,
        rising_high_hold=None,
        falling_low_hold=None,
        source_state_separated=None,
        state_transfer_hold_supported=None,
    )


def _summarise(records: Sequence[H1ContinuationStateTransferReplicate], seeds: Sequence[int]) -> dict[str, object]:
    total = len(records)
    calibrated = tuple(record for record in records if record.calibration_supported is True)
    anchored = tuple(record for record in records if record.anchor_barrier is not None)
    separated = tuple(record for record in records if record.source_state_separated is True)
    held = tuple(record for record in records if record.state_transfer_hold_supported is True)
    return {
        "denominators": {
            "total_seed_replicates": total,
            "calibration_supported_count": len(calibrated),
            "calibration_supported_probability": len(calibrated) / total,
            "interior_anchor_available_count": len(anchored),
            "interior_anchor_available_probability": len(anchored) / total,
            "source_state_separated_count": len(separated),
            "source_state_separated_probability": len(separated) / total,
            "state_transfer_hold_supported_count": len(held),
            "state_transfer_hold_supported_probability": len(held) / total,
        },
        "anchor_grid_point_count": _summary(record.anchor_grid_point_count for record in anchored),
        "hold_interaction_difference_high_minus_low": _summary(
            record.rising_high_hold.terminal_interaction_mean - record.falling_low_hold.terminal_interaction_mean
            for record in anchored
            if record.rising_high_hold is not None and record.falling_low_hold is not None
        ),
        "by_master_seed": {
            str(seed): _seed_summary(tuple(record for record in records if record.master_seed == seed))
            for seed in seeds
        },
    }


def _seed_summary(records: Sequence[H1ContinuationStateTransferReplicate]) -> dict[str, object]:
    if not records:
        return {"replicate_count": 0}
    return {
        "replicate_count": len(records),
        "calibration_supported_probability": _probability(record.calibration_supported for record in records),
        "source_state_separated_probability": _probability(record.source_state_separated for record in records),
        "state_transfer_hold_supported_probability": _probability(record.state_transfer_hold_supported for record in records),
    }


def _cell_key(cell: FiniteH1BoundaryResolutionCell) -> str:
    return json.dumps({"pair": cell.pair_index, "parameters": asdict(cell.parameters)}, sort_keys=True)


def _validate_seeds(values: Sequence[int]) -> tuple[int, ...]:
    seeds = tuple(int(value) for value in values)
    if len(seeds) < 2 or len(seeds) != len(set(seeds)) or any(seed < 0 for seed in seeds):
        raise ValueError("master_seeds must contain at least two distinct non-negative values")
    return seeds


def _mean(values: Sequence[float]) -> float:
    return sum(float(value) for value in values) / len(values)


def _probability(values: Iterable[bool | None]) -> float | None:
    observed = tuple(value for value in values if value is not None)
    return None if not observed else sum(value is True for value in observed) / len(observed)


def _summary(values: Iterable[float | int]) -> dict[str, float | None]:
    observed = tuple(float(value) for value in values)
    if not observed:
        return {"mean": None, "median": None, "minimum": None, "maximum": None}
    return {"mean": sum(observed) / len(observed), "median": median(observed), "minimum": min(observed), "maximum": max(observed)}


def _flatten(mapping: Mapping[str, object], prefix: str = "") -> dict[str, object]:
    output: dict[str, object] = {}
    for key, value in mapping.items():
        name = f"{prefix}.{key}" if prefix else key
        if isinstance(value, Mapping):
            output.update(_flatten(value, name))
        else:
            output[name] = value
    return output
