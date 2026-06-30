"""Branch-conditioned H2/H3 chain in the frozen mutation-H1 primary domain.

The chain begins with a mutation-conditioned finite H1 high full state.  It is
not restarted from a canonical q value alone.  For each declared primary cell,
the code replays the high continuation route to an interior anchor, applies the
fresh full-state hold, projects that complete state into each landscape under the
already validated conservation rule, and only then simulates the H2/H3 horizon.

Primary H3 comparison: one-large versus equal-isolated.
Migration comparison: allele-frequency mixing and F_ST modulation only; it is
not called demographic or ecological rescue.

All results are Type S for the symmetric-mutation finite closure.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from causal_model.canonical_h1_bifurcation import canonical_h1_certificate
from causal_model.finite_h1_boundary_resolution_audit import FiniteH1BoundaryResolutionCell, FiniteH1BoundaryResolutionReplicate, run_finite_h1_boundary_resolution_audit
from causal_model.finite_h1_fragment_projection_audit import FullState, ProjectionInvariants, _barrier_grid, _hold_seed, _mean, _potential, _scenario_map, project_full_state
from causal_model.finite_h1_hysteresis_audit import _parameters_from_terminal, _stage_seed, _with_uniform_initial_interaction
from causal_model.multipatch_criticality_dynamics import DynamicsParameters, SimulationSnapshot
from causal_model.multipatch_criticality_experiments import (
    SCENARIO_EQUAL_ISOLATED, SCENARIO_EQUAL_MIGRATING, SCENARIO_ONE_LARGE,
    ExperimentSpec, LandscapeScenario, ParameterCell, ReplicateSummary,
    parameters_for_cell, summarise_replicate,
)
from causal_model.mutation_h1_primary_domain import MutationH1DomainCell, domain_manifest, primary_analysis_cells
from causal_model.symmetric_allele_mutation_closure import patched_h1_mutation_runner, simulate_with_symmetric_allele_mutation

DEFAULT_PRIMARY_CHAIN_MASTER_SEEDS = (20260810, 20260811, 20260812, 20260813, 20260814)
_SCENARIOS = (SCENARIO_ONE_LARGE, SCENARIO_EQUAL_ISOLATED, SCENARIO_EQUAL_MIGRATING)


@dataclass(frozen=True)
class H2Comparison:
    warning_id: str
    baseline_eligible: bool
    warning_time: int | None
    trait_loss_time: int | None
    valid_pair: bool
    censored: bool
    warning_leads: bool | None

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ScenarioOutcome:
    scenario_id: str
    projection_supported: bool
    baseline_h_alpha: float
    baseline_h_gamma: float
    baseline_genetic_eligible: bool
    summary: ReplicateSummary
    h2: tuple[H2Comparison, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "scenario_id": self.scenario_id,
            "projection_supported": self.projection_supported,
            "baseline_h_alpha": self.baseline_h_alpha,
            "baseline_h_gamma": self.baseline_h_gamma,
            "baseline_genetic_eligible": self.baseline_genetic_eligible,
            "summary": self.summary.as_dict(),
            "h2": [value.as_dict() for value in self.h2],
        }


@dataclass(frozen=True)
class PrimaryChainReplicate:
    mutation_rate: float
    area_reference: float
    interaction_feedback: float
    master_seed: int
    replicate_index: int
    calibration_seed: int
    h1_full_state_hold_supported: bool | None
    anchor_barrier: float | None
    projections: Mapping[str, ProjectionInvariants] | None
    outcomes: Mapping[str, ScenarioOutcome] | None
    h3_fragmentation_pattern_supported: bool | None
    migration_fst_lower_than_isolation: bool | None
    same_replicate_halpha_chain_supported: bool | None
    same_replicate_hgamma_chain_supported: bool | None

    def as_dict(self) -> dict[str, object]:
        return {
            "mutation_rate": self.mutation_rate,
            "area_reference": self.area_reference,
            "interaction_feedback": self.interaction_feedback,
            "master_seed": self.master_seed,
            "replicate_index": self.replicate_index,
            "calibration_seed": self.calibration_seed,
            "h1_full_state_hold_supported": self.h1_full_state_hold_supported,
            "anchor_barrier": self.anchor_barrier,
            "projections": None if self.projections is None else {key: value.as_dict() for key, value in self.projections.items()},
            "outcomes": None if self.outcomes is None else {key: value.as_dict() for key, value in self.outcomes.items()},
            "h3_fragmentation_pattern_supported": self.h3_fragmentation_pattern_supported,
            "migration_fst_lower_than_isolation": self.migration_fst_lower_than_isolation,
            "same_replicate_halpha_chain_supported": self.same_replicate_halpha_chain_supported,
            "same_replicate_hgamma_chain_supported": self.same_replicate_hgamma_chain_supported,
        }


@dataclass(frozen=True)
class PrimaryChainCell:
    domain_cell: MutationH1DomainCell
    master_seeds: tuple[int, ...]
    records: tuple[PrimaryChainReplicate, ...]
    summary: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "domain_cell": self.domain_cell.as_dict(),
            "master_seeds": list(self.master_seeds),
            "records": [record.as_dict() for record in self.records],
            "summary": dict(self.summary),
        }

    def to_csv_row(self) -> dict[str, object]:
        row = self.domain_cell.as_dict()
        row.update({"master_seeds": ",".join(str(seed) for seed in self.master_seeds), "record_count": len(self.records)})
        row.update(_flatten(self.summary))
        return row


def run_mutation_primary_h1_h2_h3_chain(
    spec: ExperimentSpec,
    *,
    master_seeds: Sequence[int] = DEFAULT_PRIMARY_CHAIN_MASTER_SEEDS,
    endpoint_padding_fraction: float = 0.5,
    stage_generations: int = 30,
    hold_generations: int = 30,
    nested_barrier_points: Sequence[int] = (25, 49, 97),
    interaction_separation_threshold: float = 0.05,
    maximum_normalized_bracket_width: float = 0.03,
) -> tuple[PrimaryChainCell, ...]:
    """Run H2/H3 only for the frozen independent-validation primary domain."""
    seeds = _validate_seeds(master_seeds)
    if spec.replicates < 1 or spec.generations < 1:
        raise ValueError("spec requires positive generations and replicates")
    cells: list[PrimaryChainCell] = []
    for domain in primary_analysis_cells():
        records: list[PrimaryChainReplicate] = []
        for master_seed in seeds:
            target_spec = replace(
                spec,
                master_seed=master_seed,
                area_reference_values=(domain.area_reference,),
                interaction_feedback_values=(domain.interaction_feedback,),
                interaction_barrier_values=(0.5,),
            )
            with patched_h1_mutation_runner(domain.mutation_rate):
                calibration = run_finite_h1_boundary_resolution_audit(
                    target_spec,
                    endpoint_padding_fraction=endpoint_padding_fraction,
                    stage_generations=stage_generations,
                    nested_barrier_points=nested_barrier_points,
                    interaction_separation_threshold=interaction_separation_threshold,
                    maximum_normalized_bracket_width=maximum_normalized_bracket_width,
                )
            if len(calibration) != 1:
                raise RuntimeError("targeted mutation H1 calibration must yield exactly one parameter cell")
            source_cell = calibration[0]
            scenarios = _scenario_map(target_spec)
            records.extend(
                _run_record(
                    domain, target_spec, source_cell, calibration_record, scenarios,
                    master_seed=master_seed, endpoint_padding_fraction=endpoint_padding_fraction,
                    stage_generations=stage_generations, hold_generations=hold_generations,
                    interaction_separation_threshold=interaction_separation_threshold,
                )
                for calibration_record in source_cell.replicates
            )
        cells.append(PrimaryChainCell(domain, seeds, tuple(records), _summarise(records)))
    return tuple(cells)


def _run_record(
    domain: MutationH1DomainCell,
    spec: ExperimentSpec,
    calibration: FiniteH1BoundaryResolutionCell,
    record: FiniteH1BoundaryResolutionReplicate,
    scenarios: Mapping[str, LandscapeScenario],
    *, master_seed: int, endpoint_padding_fraction: float, stage_generations: int,
    hold_generations: int, interaction_separation_threshold: float,
) -> PrimaryChainReplicate:
    prepared = _prepare_mutation_high_state(
        domain.mutation_rate, spec, calibration, record,
        endpoint_padding_fraction=endpoint_padding_fraction, stage_generations=stage_generations,
        hold_generations=hold_generations, interaction_separation_threshold=interaction_separation_threshold,
    )
    base = dict(mutation_rate=domain.mutation_rate, area_reference=domain.area_reference,
                interaction_feedback=domain.interaction_feedback, master_seed=master_seed,
                replicate_index=record.replicate_index, calibration_seed=record.seed)
    if prepared is None:
        return PrimaryChainReplicate(**base, h1_full_state_hold_supported=False if record.resolution_stable_h1_loop_mechanism_supported else None,
            anchor_barrier=None, projections=None, outcomes=None, h3_fragmentation_pattern_supported=None,
            migration_fst_lower_than_isolation=None, same_replicate_halpha_chain_supported=None,
            same_replicate_hgamma_chain_supported=None)
    source, anchor = prepared
    projections: dict[str, ProjectionInvariants] = {}
    outcomes: dict[str, ScenarioOutcome] = {}
    anchor_cell = replace(calibration.parameters, interaction_barrier=anchor)
    for scenario_id in _SCENARIOS:
        template = parameters_for_cell(spec, scenarios[scenario_id], anchor_cell, seed=_outcome_seed(record.seed, scenario_id))
        projected, invariant = project_full_state(source, template)
        projections[scenario_id] = invariant
        result = simulate_with_symmetric_allele_mutation(
            replace(projected, generations=spec.generations, random_seed=_outcome_seed(record.seed, scenario_id)),
            mutation_rate=domain.mutation_rate,
        )
        initial = result.snapshots[0]
        baseline_eligible = initial.h_alpha > spec.h_alpha_warning_threshold and initial.h_gamma > spec.h_gamma_warning_threshold
        summary = summarise_replicate(result, replicate_index=record.replicate_index, seed=_outcome_seed(record.seed, scenario_id),
            h_alpha_warning_threshold=spec.h_alpha_warning_threshold, h_gamma_warning_threshold=spec.h_gamma_warning_threshold,
            fst_warning_threshold=spec.fst_warning_threshold, allele_loss_threshold=spec.allele_loss_threshold)
        outcomes[scenario_id] = ScenarioOutcome(
            scenario_id, invariant.projection_supported, initial.h_alpha, initial.h_gamma, baseline_eligible, summary,
            (_h2("H_alpha", summary.tau_H_alpha, summary.tau_trait_realised, baseline_eligible),
             _h2("H_gamma", summary.tau_H_gamma, summary.tau_trait_realised, baseline_eligible)),
        )
    large, isolated, migrating = (outcomes[SCENARIO_ONE_LARGE], outcomes[SCENARIO_EQUAL_ISOLATED], outcomes[SCENARIO_EQUAL_MIGRATING])
    h3 = (
        isolated.summary.final_q_by_patch and
        _mean(isolated.summary.final_q_by_patch) < _mean(large.summary.final_q_by_patch) and
        _mean(isolated.summary.final_effective_size_by_patch) < _mean(large.summary.final_effective_size_by_patch) and
        isolated.summary.realised_high_trait_mass_mean < large.summary.realised_high_trait_mass_mean
    )
    iso_h2 = {item.warning_id: item for item in isolated.h2}
    migration_fst = None if migrating.summary.fst is None or isolated.summary.fst is None else migrating.summary.fst < isolated.summary.fst
    return PrimaryChainReplicate(**base, h1_full_state_hold_supported=True, anchor_barrier=anchor,
        projections=projections, outcomes=outcomes, h3_fragmentation_pattern_supported=h3,
        migration_fst_lower_than_isolation=migration_fst,
        same_replicate_halpha_chain_supported=h3 and iso_h2["H_alpha"].warning_leads is True,
        same_replicate_hgamma_chain_supported=h3 and iso_h2["H_gamma"].warning_leads is True)


def _prepare_mutation_high_state(
    rate: float, spec: ExperimentSpec, calibration: FiniteH1BoundaryResolutionCell,
    record: FiniteH1BoundaryResolutionReplicate, *, endpoint_padding_fraction: float,
    stage_generations: int, hold_generations: int, interaction_separation_threshold: float,
) -> tuple[FullState, float] | None:
    if record.resolution_stable_h1_loop_mechanism_supported is not True:
        return None
    observation, interval = record.observations[-1], calibration.canonical_bistable_barrier_interval
    collapse, recovery = observation.rising_collapse_bracket, observation.falling_recovery_bracket
    if interval is None or collapse is None or recovery is None or not observation.finite_h1_loop_bracket_mechanism_supported:
        return None
    barriers = _barrier_grid(interval, barrier_points=observation.barrier_points, padding=endpoint_padding_fraction * (interval[1] - interval[0]))
    candidates = tuple(value for value in barriers if recovery.upper_barrier < value < collapse.lower_barrier)
    if not candidates:
        return None
    anchor = min(candidates, key=lambda value: (abs(value - (recovery.upper_barrier + collapse.lower_barrier) / 2.0), value))
    one_large = _scenario_map(spec)[SCENARIO_ONE_LARGE]
    base = parameters_for_cell(spec, one_large, calibration.parameters, seed=record.seed)
    canonical = canonical_h1_certificate(calibration.parameters.interaction_feedback, one_large.patch_areas[0], calibration.parameters.area_reference,
        (interval[0] + interval[1]) / 2.0, base)
    if canonical.high_stable_branch is None or canonical.low_stable_branch is None:
        return None
    high_terminal, high_carried = _replay(rate, _with_uniform_initial_interaction(base, canonical.high_stable_branch.interaction), barriers, 1, record.seed, stage_generations, anchor)
    low_terminal, low_carried = _replay(rate, _with_uniform_initial_interaction(base, canonical.low_stable_branch.interaction), tuple(reversed(barriers)), 2, record.seed, stage_generations, anchor)
    high_hold = simulate_with_symmetric_allele_mutation(replace(high_carried, interaction_barrier=anchor, generations=hold_generations, random_seed=_hold_seed(record.seed, 1)), mutation_rate=rate).snapshots[-1]
    low_hold = simulate_with_symmetric_allele_mutation(replace(low_carried, interaction_barrier=anchor, generations=hold_generations, random_seed=_hold_seed(record.seed, 2)), mutation_rate=rate).snapshots[-1]
    if not (_mean(high_terminal.interaction) - _mean(low_terminal.interaction) > interaction_separation_threshold and _potential(high_terminal) and not _potential(low_terminal)
            and _mean(high_hold.interaction) - _mean(low_hold.interaction) > interaction_separation_threshold and _potential(high_hold) and not _potential(low_hold)):
        return None
    return FullState.from_snapshot(high_hold, one_large.patch_areas), anchor


def _replay(rate: float, initial: DynamicsParameters, barriers: Sequence[float], route: int, seed: int, generations: int, anchor: float) -> tuple[SimulationSnapshot, DynamicsParameters]:
    parameters = initial
    for index, barrier in enumerate(barriers):
        terminal = simulate_with_symmetric_allele_mutation(replace(parameters, interaction_barrier=barrier, generations=generations, random_seed=_stage_seed(seed, route, index)), mutation_rate=rate).snapshots[-1]
        carried = _parameters_from_terminal(parameters, terminal)
        if barrier == anchor:
            return terminal, carried
        parameters = carried
    raise RuntimeError("anchor absent from reconstructed barrier grid")


def _h2(name: str, warning: int | None, trait: int | None, eligible: bool) -> H2Comparison:
    valid = eligible and warning is not None and trait is not None
    return H2Comparison(name, eligible, warning, trait, valid, not valid, None if not valid else warning < trait)


def _outcome_seed(seed: int, scenario_id: str) -> int:
    code = {SCENARIO_ONE_LARGE: 11, SCENARIO_EQUAL_ISOLATED: 23, SCENARIO_EQUAL_MIGRATING: 37}[scenario_id]
    return (seed * 1_000_003 + 90_001 + code) % (2**31 - 1)


def _validate_seeds(values: Sequence[int]) -> tuple[int, ...]:
    result = tuple(int(value) for value in values)
    if len(result) < 2 or len(result) != len(set(result)) or any(value < 0 for value in result):
        raise ValueError("master_seeds must contain at least two distinct non-negative values")
    return result


def _summarise(records: Sequence[PrimaryChainReplicate]) -> dict[str, object]:
    total = len(records)
    available = tuple(record for record in records if record.h1_full_state_hold_supported is True)
    return {
        "denominators": {"total_seed_replicates": total, "h1_full_state_hold_supported_count": len(available), "h1_full_state_hold_supported_probability": len(available) / total},
        "h3_fragmentation_pattern_supported_probability": None if not available else sum(record.h3_fragmentation_pattern_supported is True for record in available) / len(available),
        "same_replicate_halpha_chain_supported_probability": None if not available else sum(record.same_replicate_halpha_chain_supported is True for record in available) / len(available),
        "same_replicate_hgamma_chain_supported_probability": None if not available else sum(record.same_replicate_hgamma_chain_supported is True for record in available) / len(available),
        "migration_fst_lower_than_isolation_probability": None if not available else sum(record.migration_fst_lower_than_isolation is True for record in available if record.migration_fst_lower_than_isolation is not None) / max(1, sum(record.migration_fst_lower_than_isolation is not None for record in available)),
    }


def write_mutation_primary_h1_h2_h3_artifacts(cells: Iterable[PrimaryChainCell], *, csv_path: str | Path, json_path: str | Path, manifest_path: str | Path) -> None:
    values = tuple(cells)
    if not values:
        raise ValueError("cells must be nonempty")
    csv_target, json_target, manifest_target = Path(csv_path), Path(json_path), Path(manifest_path)
    for target in (csv_target, json_target, manifest_target):
        target.parent.mkdir(parents=True, exist_ok=True)
    rows = [cell.to_csv_row() for cell in values]
    with csv_target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted({key for row in rows for key in row}))
        writer.writeheader(); writer.writerows(rows)
    with json_target.open("w", encoding="utf-8") as handle:
        json.dump([cell.as_dict() for cell in values], handle, indent=2, sort_keys=True)
    with manifest_target.open("w", encoding="utf-8") as handle:
        json.dump({"campaign": "mutation_primary_full_state_h1_h2_h3_v1", "primary_domain": domain_manifest(), "primary_h3_comparison": "one_large vs equal_isolated", "migration_interpretation": "allele-frequency mixing only, not rescue"}, handle, indent=2, sort_keys=True)


def _flatten(mapping: Mapping[str, object], prefix: str = "") -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in mapping.items():
        name = f"{prefix}.{key}" if prefix else key
        if isinstance(value, Mapping): result.update(_flatten(value, name))
        else: result[name] = value
    return result
