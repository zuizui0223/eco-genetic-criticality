"""Conservation-preserving projection of finite H1 full states into fragments.

A one-large H1 state cannot be copied into equal fragments by repeating every
absolute state variable.  That would create individuals and realised trait-bin
abundance from nothing.  This module fixes a declared projection rule before a
fragmentation H2/H3 campaign uses it.

For a one-large full state projected to target patch areas, it preserves:

* total population exactly;
* the global total of every realised trait bin exactly;
* population-weighted high-allele frequency exactly; and
* population-weighted interaction q by copying its intensive mean to each
  target patch.

Population and trait-bin counts are allocated by deterministic largest-remainder
rounding subject to both row (target patch population) and column (trait-bin)
margins.  The audit replays the validated H1 high branch, performs its fresh
full-state hold, projects that source into one-large, equal-isolated, and
equal-migrating scenarios, and checks the actual simulator initial snapshot.

This is a Type S construction audit.  It validates the state-transfer bridge;
it does not yet make H2 warning or H3 outcome claims.
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
from causal_model.multipatch_criticality_dynamics import (
    DynamicsParameters,
    SimulationSnapshot,
    simulate,
)
from causal_model.multipatch_criticality_experiments import (
    SCENARIO_EQUAL_ISOLATED,
    SCENARIO_EQUAL_MIGRATING,
    SCENARIO_ONE_LARGE,
    ExperimentSpec,
    LandscapeScenario,
    ParameterCell,
    default_scenarios,
    parameters_for_cell,
)

DEFAULT_MASTER_SEEDS = (20260630, 20260631, 20260632, 20260633, 20260634)
SCENARIO_IDS = (SCENARIO_ONE_LARGE, SCENARIO_EQUAL_ISOLATED, SCENARIO_EQUAL_MIGRATING)


@dataclass(frozen=True)
class FullState:
    """Initialisable finite state, with absolute and intensive quantities explicit."""

    patch_areas: tuple[float, ...]
    population: tuple[int, ...]
    interaction: tuple[float, ...]
    high_allele_frequency: tuple[float, ...]
    trait_abundance: tuple[tuple[int, ...], ...]

    @classmethod
    def from_snapshot(cls, snapshot: SimulationSnapshot, patch_areas: Sequence[float]) -> "FullState":
        return cls(
            patch_areas=tuple(float(area) for area in patch_areas),
            population=tuple(int(value) for value in snapshot.population),
            interaction=tuple(float(value) for value in snapshot.interaction),
            high_allele_frequency=tuple(float(value) for value in snapshot.high_allele_frequency),
            trait_abundance=tuple(tuple(int(value) for value in item.abundance) for item in snapshot.trait_occupancy),
        )

    @property
    def total_population(self) -> int:
        return sum(self.population)

    @property
    def trait_bin_totals(self) -> tuple[int, ...]:
        return tuple(sum(row[index] for row in self.trait_abundance) for index in range(len(self.trait_abundance[0])))

    def population_weighted(self, values: Sequence[float]) -> float:
        return sum(population * value for population, value in zip(self.population, values)) / self.total_population


@dataclass(frozen=True)
class ProjectionInvariants:
    """Source-to-target conservation checks measured from simulator initial state."""

    source_total_population: int
    target_total_population: int
    population_conserved: bool
    source_weighted_interaction: float
    target_weighted_interaction: float
    interaction_replicated: bool
    source_weighted_high_allele_frequency: float
    target_weighted_high_allele_frequency: float
    allele_frequency_conserved: bool
    source_trait_bin_totals: tuple[int, ...]
    target_trait_bin_totals: tuple[int, ...]
    trait_bins_conserved: bool
    source_total_area: float
    target_total_area: float
    area_conserved: bool
    projection_supported: bool

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ProjectionScenarioRecord:
    """A declared scenario-specific projection of the same H1 high state."""

    scenario_id: str
    target_patch_areas: tuple[float, ...]
    target_initial_population: tuple[int, ...]
    target_initial_interaction: tuple[float, ...]
    target_initial_high_allele_frequency: tuple[float, ...]
    invariants: ProjectionInvariants

    def as_dict(self) -> dict[str, object]:
        return {
            "scenario_id": self.scenario_id,
            "target_patch_areas": list(self.target_patch_areas),
            "target_initial_population": list(self.target_initial_population),
            "target_initial_interaction": list(self.target_initial_interaction),
            "target_initial_high_allele_frequency": list(self.target_initial_high_allele_frequency),
            "invariants": self.invariants.as_dict(),
        }


@dataclass(frozen=True)
class FragmentProjectionReplicate:
    """One H1-conditioned high full state and all target projections."""

    master_seed: int
    replicate_index: int
    calibration_seed: int
    h1_full_state_hold_supported: bool | None
    anchor_barrier: float | None
    source_state: FullState | None
    scenario_projections: Mapping[str, ProjectionScenarioRecord] | None

    def as_dict(self) -> dict[str, object]:
        return {
            "master_seed": self.master_seed,
            "replicate_index": self.replicate_index,
            "calibration_seed": self.calibration_seed,
            "h1_full_state_hold_supported": self.h1_full_state_hold_supported,
            "anchor_barrier": self.anchor_barrier,
            "source_state": None if self.source_state is None else asdict(self.source_state),
            "scenario_projections": None
            if self.scenario_projections is None
            else {name: record.as_dict() for name, record in self.scenario_projections.items()},
        }


@dataclass(frozen=True)
class FragmentProjectionCell:
    """Projection evidence across seeds for a canonical H1 parameter pair."""

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
    replicates: tuple[FragmentProjectionReplicate, ...]
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
        row: dict[str, object] = {
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


def project_full_state(
    source: FullState,
    target_template: DynamicsParameters,
) -> tuple[DynamicsParameters, ProjectionInvariants]:
    """Project a source state to target patches under the declared conservation rule."""
    target_areas = tuple(float(value) for value in target_template.patch_areas)
    if len(source.trait_abundance) != len(source.population):
        raise ValueError("source trait_abundance must match source patch count")
    if not source.trait_abundance or not source.trait_abundance[0]:
        raise ValueError("source trait_abundance must be nonempty")
    if any(len(row) != len(source.trait_abundance[0]) for row in source.trait_abundance):
        raise ValueError("source trait-bin rows must have equal length")
    if sum(source.trait_bin_totals) != source.total_population:
        raise ValueError("source trait-bin totals must equal source total population")
    if abs(sum(source.patch_areas) - sum(target_areas)) > 1e-12:
        raise ValueError("source and target total area must match for fragmentation projection")

    target_population = _positive_largest_remainder(source.total_population, target_areas)
    target_trait_abundance = _contingency_allocate(target_population, source.trait_bin_totals)
    interaction_mean = source.population_weighted(source.interaction)
    allele_frequency_mean = source.population_weighted(source.high_allele_frequency)
    projected = replace(
        target_template,
        initial_population=target_population,
        initial_interaction=tuple(interaction_mean for _ in target_population),
        initial_high_allele_frequency=tuple(allele_frequency_mean for _ in target_population),
        initial_trait_distribution=(),
        initial_trait_abundance=target_trait_abundance,
    )
    initial = simulate(replace(projected, generations=1)).snapshots[0]
    target_bins = tuple(
        sum(item.abundance[index] for item in initial.trait_occupancy)
        for index in range(len(source.trait_bin_totals))
    )
    target_q = _population_weighted(initial.population, initial.interaction)
    target_p = _population_weighted(initial.population, initial.high_allele_frequency)
    invariants = ProjectionInvariants(
        source_total_population=source.total_population,
        target_total_population=sum(initial.population),
        population_conserved=sum(initial.population) == source.total_population,
        source_weighted_interaction=interaction_mean,
        target_weighted_interaction=target_q,
        interaction_replicated=(
            all(abs(value - interaction_mean) <= 1e-12 for value in initial.interaction)
            and abs(target_q - interaction_mean) <= 1e-12
        ),
        source_weighted_high_allele_frequency=allele_frequency_mean,
        target_weighted_high_allele_frequency=target_p,
        allele_frequency_conserved=abs(target_p - allele_frequency_mean) <= 1e-12,
        source_trait_bin_totals=source.trait_bin_totals,
        target_trait_bin_totals=target_bins,
        trait_bins_conserved=target_bins == source.trait_bin_totals,
        source_total_area=sum(source.patch_areas),
        target_total_area=sum(target_areas),
        area_conserved=abs(sum(source.patch_areas) - sum(target_areas)) <= 1e-12,
        projection_supported=False,
    )
    invariants = replace(
        invariants,
        projection_supported=(
            invariants.population_conserved
            and invariants.interaction_replicated
            and invariants.allele_frequency_conserved
            and invariants.trait_bins_conserved
            and invariants.area_conserved
        ),
    )
    return projected, invariants


def run_finite_h1_fragment_projection_audit(
    spec: ExperimentSpec,
    *,
    master_seeds: Sequence[int] = DEFAULT_MASTER_SEEDS,
    endpoint_padding_fraction: float = 0.5,
    stage_generations: int = 30,
    hold_generations: int = 30,
    nested_barrier_points: Sequence[int] = DEFAULT_NESTED_BARRIER_POINTS,
    interaction_separation_threshold: float = 0.05,
    maximum_normalized_bracket_width: float = 0.03,
) -> tuple[FragmentProjectionCell, ...]:
    """Prepare H1 high states then validate projections into all H3 landscapes."""
    seeds = _validate_master_seeds(master_seeds)
    scenarios = _scenario_map(spec)
    calibrations = tuple(
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
    cell_count = len(calibrations[0])
    if any(len(campaign) != cell_count for campaign in calibrations):
        raise RuntimeError("H1 calibration cell counts differ across master seeds")
    output: list[FragmentProjectionCell] = []
    for cell_index in range(cell_count):
        aligned = tuple(campaign[cell_index] for campaign in calibrations)
        reference = aligned[0]
        if any(_cell_key(value) != _cell_key(reference) for value in aligned[1:]):
            raise RuntimeError("H1 calibration cells are misaligned across master seeds")
        records: list[FragmentProjectionReplicate] = []
        for master_seed, calibration in zip(seeds, aligned, strict=True):
            seed_spec = replace(spec, master_seed=master_seed)
            for calibration_replicate in calibration.replicates:
                records.append(
                    _projection_record(
                        seed_spec,
                        calibration,
                        calibration_replicate,
                        scenarios,
                        master_seed=master_seed,
                        endpoint_padding_fraction=endpoint_padding_fraction,
                        stage_generations=stage_generations,
                        hold_generations=hold_generations,
                        interaction_separation_threshold=interaction_separation_threshold,
                    )
                )
        output.append(
            FragmentProjectionCell(
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


def write_finite_h1_fragment_projection_artifacts(
    cells: Iterable[FragmentProjectionCell],
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


def _projection_record(
    spec: ExperimentSpec,
    calibration: FiniteH1BoundaryResolutionCell,
    calibration_replicate: FiniteH1BoundaryResolutionReplicate,
    scenarios: Mapping[str, LandscapeScenario],
    *,
    master_seed: int,
    endpoint_padding_fraction: float,
    stage_generations: int,
    hold_generations: int,
    interaction_separation_threshold: float,
) -> FragmentProjectionReplicate:
    prepared = _prepare_high_source(
        spec,
        calibration,
        calibration_replicate,
        endpoint_padding_fraction=endpoint_padding_fraction,
        stage_generations=stage_generations,
        hold_generations=hold_generations,
        interaction_separation_threshold=interaction_separation_threshold,
    )
    if prepared is None:
        return FragmentProjectionReplicate(
            master_seed=master_seed,
            replicate_index=calibration_replicate.replicate_index,
            calibration_seed=calibration_replicate.seed,
            h1_full_state_hold_supported=False if calibration_replicate.resolution_stable_h1_loop_mechanism_supported else None,
            anchor_barrier=None,
            source_state=None,
            scenario_projections=None,
        )
    source, anchor = prepared
    projected: dict[str, ProjectionScenarioRecord] = {}
    for scenario_id in SCENARIO_IDS:
        scenario = scenarios[scenario_id]
        template = parameters_for_cell(spec, scenario, calibration.parameters, seed=calibration_replicate.seed)
        parameters, invariants = project_full_state(source, template)
        projected[scenario_id] = ProjectionScenarioRecord(
            scenario_id=scenario_id,
            target_patch_areas=parameters.patch_areas,
            target_initial_population=parameters.initial_population,
            target_initial_interaction=parameters.initial_interaction,
            target_initial_high_allele_frequency=parameters.initial_high_allele_frequency,
            invariants=invariants,
        )
    return FragmentProjectionReplicate(
        master_seed=master_seed,
        replicate_index=calibration_replicate.replicate_index,
        calibration_seed=calibration_replicate.seed,
        h1_full_state_hold_supported=True,
        anchor_barrier=anchor,
        source_state=source,
        scenario_projections=projected,
    )


def _prepare_high_source(
    spec: ExperimentSpec,
    calibration: FiniteH1BoundaryResolutionCell,
    calibration_replicate: FiniteH1BoundaryResolutionReplicate,
    *,
    endpoint_padding_fraction: float,
    stage_generations: int,
    hold_generations: int,
    interaction_separation_threshold: float,
) -> tuple[FullState, float] | None:
    if calibration_replicate.resolution_stable_h1_loop_mechanism_supported is not True:
        return None
    observation = calibration_replicate.observations[-1]
    collapse, recovery = observation.rising_collapse_bracket, observation.falling_recovery_bracket
    interval = calibration.canonical_bistable_barrier_interval
    if not observation.finite_h1_loop_bracket_mechanism_supported or collapse is None or recovery is None or interval is None:
        return None
    lower, upper = recovery.upper_barrier, collapse.lower_barrier
    barriers = _barrier_grid(
        interval,
        barrier_points=observation.barrier_points,
        padding=endpoint_padding_fraction * (interval[1] - interval[0]),
    )
    candidates = tuple(barrier for barrier in barriers if lower < barrier < upper)
    if not candidates:
        return None
    anchor = min(candidates, key=lambda value: (abs(value - (lower + upper) / 2.0), value))
    one_large = _scenario_map(spec)[SCENARIO_ONE_LARGE]
    base = parameters_for_cell(spec, one_large, calibration.parameters, seed=calibration_replicate.seed)
    canonical = canonical_h1_certificate(
        feedback_strength=calibration.parameters.interaction_feedback,
        area=one_large.patch_areas[0],
        area_reference=calibration.parameters.area_reference,
        barrier=(interval[0] + interval[1]) / 2.0,
        trait_parameters=base,
    )
    if canonical.high_stable_branch is None or canonical.low_stable_branch is None:
        return None
    high_terminal, high_carried = _replay_to_anchor(
        _with_uniform_initial_interaction(base, canonical.high_stable_branch.interaction),
        barriers,
        route_code=1,
        replicate_seed=calibration_replicate.seed,
        stage_generations=stage_generations,
        anchor=anchor,
    )
    low_terminal, low_carried = _replay_to_anchor(
        _with_uniform_initial_interaction(base, canonical.low_stable_branch.interaction),
        tuple(reversed(barriers)),
        route_code=2,
        replicate_seed=calibration_replicate.seed,
        stage_generations=stage_generations,
        anchor=anchor,
    )
    high_hold = simulate(replace(high_carried, interaction_barrier=anchor, generations=hold_generations, random_seed=_hold_seed(calibration_replicate.seed, 1))).snapshots[-1]
    low_hold = simulate(replace(low_carried, interaction_barrier=anchor, generations=hold_generations, random_seed=_hold_seed(calibration_replicate.seed, 2))).snapshots[-1]
    supported = (
        _mean(high_terminal.interaction) - _mean(low_terminal.interaction) > interaction_separation_threshold
        and _potential(high_terminal)
        and not _potential(low_terminal)
        and _mean(high_hold.interaction) - _mean(low_hold.interaction) > interaction_separation_threshold
        and _potential(high_hold)
        and not _potential(low_hold)
    )
    if not supported:
        return None
    return FullState.from_snapshot(high_hold, one_large.patch_areas), anchor


def _replay_to_anchor(
    initial: DynamicsParameters,
    barriers: Sequence[float],
    *,
    route_code: int,
    replicate_seed: int,
    stage_generations: int,
    anchor: float,
) -> tuple[SimulationSnapshot, DynamicsParameters]:
    parameters = initial
    for index, barrier in enumerate(barriers):
        result = simulate(replace(parameters, interaction_barrier=barrier, generations=stage_generations, random_seed=_stage_seed(replicate_seed, route_code, index)))
        terminal = result.snapshots[-1]
        carried = _parameters_from_terminal(parameters, terminal)
        if barrier == anchor:
            return terminal, carried
        parameters = carried
    raise RuntimeError("anchor is absent from reconstructed continuation grid")


def _positive_largest_remainder(total: int, weights: Sequence[float]) -> tuple[int, ...]:
    if total < len(weights):
        raise ValueError("total population must be at least the number of target patches")
    if not weights or any(weight <= 0.0 for weight in weights):
        raise ValueError("weights must be positive")
    residual = total - len(weights)
    weight_sum = sum(weights)
    raw = tuple(residual * weight / weight_sum for weight in weights)
    base = [int(value) for value in raw]
    remainder = residual - sum(base)
    for index in sorted(range(len(weights)), key=lambda item: (-(raw[item] - base[item]), item))[:remainder]:
        base[index] += 1
    return tuple(value + 1 for value in base)


def _contingency_allocate(row_totals: Sequence[int], column_totals: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    """Allocate integer counts with exact row and column margins."""
    if not row_totals or not column_totals or any(value < 0 for value in (*row_totals, *column_totals)):
        raise ValueError("row and column totals must be nonempty and non-negative")
    total = sum(row_totals)
    if total != sum(column_totals) or total <= 0:
        raise ValueError("row and column totals must have the same positive total")
    raw = [[row * column / total for column in column_totals] for row in row_totals]
    matrix = [[int(value) for value in values] for values in raw]
    row_deficit = [row - sum(values) for row, values in zip(row_totals, matrix)]
    column_deficit = [column - sum(matrix[row][index] for row in range(len(matrix))) for index, column in enumerate(column_totals)]
    while sum(row_deficit) > 0:
        candidates = [
            (raw[row][column] - matrix[row][column], -row, -column)
            for row, row_remaining in enumerate(row_deficit)
            if row_remaining > 0
            for column, column_remaining in enumerate(column_deficit)
            if column_remaining > 0
        ]
        if not candidates:
            raise RuntimeError("integer contingency allocation lost a feasible margin")
        _fraction, negative_row, negative_column = max(candidates)
        row, column = -negative_row, -negative_column
        matrix[row][column] += 1
        row_deficit[row] -= 1
        column_deficit[column] -= 1
    if any(column_deficit):
        raise RuntimeError("integer contingency allocation failed to fill column margins")
    return tuple(tuple(row) for row in matrix)


def _scenario_map(spec: ExperimentSpec) -> Mapping[str, LandscapeScenario]:
    values = {scenario.scenario_id: scenario for scenario in default_scenarios(spec)}
    if set(values) != set(SCENARIO_IDS):
        raise ValueError("default scenarios must contain one_large, equal_isolated, and equal_migrating")
    return values


def _population_weighted(population: Sequence[int], values: Sequence[float]) -> float:
    return sum(count * value for count, value in zip(population, values)) / sum(population)


def _potential(snapshot: SimulationSnapshot) -> bool:
    return any(item.high_trait_component_present for item in snapshot.trait_space)


def _hold_seed(seed: int, route_code: int) -> int:
    return (seed * 1_000_003 + 70_001 + route_code * 101) % (2**31 - 1)


def _cell_key(cell: FiniteH1BoundaryResolutionCell) -> str:
    return json.dumps({"pair_index": cell.pair_index, "parameters": asdict(cell.parameters)}, sort_keys=True)


def _validate_master_seeds(values: Sequence[int]) -> tuple[int, ...]:
    seeds = tuple(int(value) for value in values)
    if len(seeds) < 2 or len(seeds) != len(set(seeds)) or any(seed < 0 for seed in seeds):
        raise ValueError("master_seeds must contain at least two distinct non-negative values")
    return seeds


def _mean(values: Sequence[float]) -> float:
    return sum(float(value) for value in values) / len(values)


def _probability(values: Iterable[bool | None]) -> float | None:
    observed = tuple(value for value in values if value is not None)
    return None if not observed else sum(value is True for value in observed) / len(observed)


def _summarise(records: Sequence[FragmentProjectionReplicate], seeds: Sequence[int]) -> dict[str, object]:
    total = len(records)
    prepared = tuple(record for record in records if record.h1_full_state_hold_supported is True)
    by_scenario: dict[str, object] = {}
    for scenario_id in SCENARIO_IDS:
        values = tuple(
            record.scenario_projections[scenario_id].invariants
            for record in prepared
            if record.scenario_projections is not None
        )
        by_scenario[scenario_id] = {
            "h1_full_state_prepared_count": len(values),
            "projection_supported_count": sum(value.projection_supported for value in values),
            "projection_supported_probability_across_all_seed_replicates": sum(value.projection_supported for value in values) / total,
            "projection_supported_probability_conditional_on_h1_full_state": None if not values else sum(value.projection_supported for value in values) / len(values),
            "population_conserved_probability": None if not values else sum(value.population_conserved for value in values) / len(values),
            "trait_bins_conserved_probability": None if not values else sum(value.trait_bins_conserved for value in values) / len(values),
            "allele_frequency_conserved_probability": None if not values else sum(value.allele_frequency_conserved for value in values) / len(values),
            "interaction_replicated_probability": None if not values else sum(value.interaction_replicated for value in values) / len(values),
        }
    return {
        "denominators": {
            "total_seed_replicates": total,
            "h1_full_state_prepared_count": len(prepared),
            "h1_full_state_prepared_probability": len(prepared) / total,
        },
        "by_scenario": by_scenario,
        "by_master_seed": {
            str(seed): {
                "replicate_count": len(tuple(record for record in records if record.master_seed == seed)),
                "h1_full_state_prepared_probability": _probability(record.h1_full_state_hold_supported for record in records if record.master_seed == seed),
            }
            for seed in seeds
        },
    }


def _flatten(mapping: Mapping[str, object], prefix: str = "") -> dict[str, object]:
    output: dict[str, object] = {}
    for key, value in mapping.items():
        name = f"{prefix}.{key}" if prefix else key
        if isinstance(value, Mapping):
            output.update(_flatten(value, name))
        else:
            output[name] = value
    return output
