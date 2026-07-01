"""Trait-loss-only calibration for the conditional H2-R proposition.

This module does not evaluate H-alpha, H-gamma, relative warnings, warning lead
counts, or lead-time magnitudes.  Its sole role is to choose, before H2-R
validation, an externally imposed monotone interaction-support deterioration
schedule that yields observable post-baseline realised high-trait loss.

For each frozen primary mutation-H1 cell, a new-seed H1 high full state is
reconstructed, held, and projected with the existing conservation rule into the
equal-isolated landscape.  The only changed quantity thereafter is the
interaction barrier in the q update.  It increases linearly by a declared
multiple of that cell's canonical bistable-interval width.

Selection uses only post-baseline trait-loss availability per seed block.  It
cannot inspect genetic-warning outcomes because they are neither calculated nor
written by this module.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

# Importing this adapter first patches the chain module's local canonical
# certificate reference to its keyword-only API adapter.  The H1 algorithm is
# unchanged; this merely preserves the already merged runtime compatibility fix.
import causal_model.mutation_primary_h1_h2_h3_runtime  # noqa: F401
from causal_model.finite_h1_boundary_resolution_audit import (
    FiniteH1BoundaryResolutionCell,
    FiniteH1BoundaryResolutionReplicate,
    run_finite_h1_boundary_resolution_audit,
)
from causal_model.finite_h1_fragment_projection_audit import (
    ProjectionInvariants,
    _scenario_map,
    project_full_state,
)
from causal_model.h2_relative_warning_contract import (
    DEFAULT_CALIBRATION_HORIZONS,
    DEFAULT_TOTAL_NORMALIZED_BARRIER_INCREASES,
    DeteriorationCalibrationCandidate,
    TARGET_TRAIT_LOSS_PROBABILITY_INTERVAL,
    h2r_protocol_manifest,
    select_trait_loss_only_calibration,
)
from causal_model.multipatch_criticality_dynamics import tau_trait_realised
from causal_model.multipatch_criticality_experiments import (
    SCENARIO_EQUAL_ISOLATED,
    ExperimentSpec,
    LandscapeScenario,
    ParameterCell,
    parameters_for_cell,
)
from causal_model.mutation_h1_primary_domain import (
    MutationH1DomainCell,
    domain_manifest,
    primary_analysis_cells,
)
from causal_model.mutation_primary_h1_h2_h3_chain import _prepare_mutation_high_state
from causal_model.symmetric_allele_mutation_closure import (
    simulate_with_symmetric_allele_mutation,
)

DEFAULT_H2R_CALIBRATION_MASTER_SEEDS = (20260910, 20260911, 20260912, 20260913, 20260914)


@dataclass(frozen=True)
class H2RTraitLossCalibrationRecord:
    """One attempted source-plus-deterioration trajectory; no warning fields exist."""

    mutation_rate: float
    area_reference: float
    interaction_feedback: float
    master_seed: int
    replicate_index: int
    calibration_seed: int
    horizon: int
    total_normalized_barrier_increase: float
    h1_resolution_supported: bool | None
    h1_full_state_source_prepared: bool
    anchor_barrier: float | None
    canonical_interval_width: float | None
    projection_supported: bool | None
    baseline_realised_high_trait_present: bool | None
    trajectory_seed: int | None
    barrier_first_generation: float | None
    barrier_final_generation: float | None
    trait_loss_time_post_baseline: int | None

    @property
    def eligible_for_trait_loss_denominator(self) -> bool:
        return (
            self.h1_full_state_source_prepared
            and self.projection_supported is True
            and self.baseline_realised_high_trait_present is True
        )

    @property
    def trait_loss_observed_post_baseline(self) -> bool | None:
        if not self.eligible_for_trait_loss_denominator:
            return None
        return self.trait_loss_time_post_baseline is not None

    def as_dict(self) -> dict[str, object]:
        result = asdict(self)
        result["eligible_for_trait_loss_denominator"] = self.eligible_for_trait_loss_denominator
        result["trait_loss_observed_post_baseline"] = self.trait_loss_observed_post_baseline
        return result


@dataclass(frozen=True)
class H2RTraitLossScheduleSummary:
    """Trait-loss-only evidence for one predeclared deterioration schedule."""

    horizon: int
    total_normalized_barrier_increase: float
    total_attempted_records: int
    h1_full_state_source_prepared_count: int
    projection_supported_count: int
    baseline_eligible_count: int
    trait_loss_count: int
    seed_block_baseline_eligible_counts: tuple[int, ...]
    seed_block_trait_loss_counts: tuple[int, ...]
    seed_block_trait_loss_probabilities: tuple[float | None, ...]

    @property
    def pooled_trait_loss_probability(self) -> float | None:
        if self.baseline_eligible_count == 0:
            return None
        return self.trait_loss_count / self.baseline_eligible_count

    @property
    def selectable_candidate(self) -> DeteriorationCalibrationCandidate | None:
        if any(value is None for value in self.seed_block_trait_loss_probabilities):
            return None
        return DeteriorationCalibrationCandidate(
            horizon=self.horizon,
            total_normalized_barrier_increase=self.total_normalized_barrier_increase,
            trait_loss_probability_by_seed_block=tuple(
                float(value) for value in self.seed_block_trait_loss_probabilities
            ),
        )

    def as_dict(self) -> dict[str, object]:
        result = asdict(self)
        result["pooled_trait_loss_probability"] = self.pooled_trait_loss_probability
        result["selection_candidate_available"] = self.selectable_candidate is not None
        result["selection_uses_warning_outcomes"] = False
        return result


@dataclass(frozen=True)
class H2RTraitLossCalibrationCell:
    """Complete calibration ledger and trait-loss-only schedule decision for one cell."""

    domain_cell: MutationH1DomainCell
    master_seeds: tuple[int, ...]
    records: tuple[H2RTraitLossCalibrationRecord, ...]
    schedules: tuple[H2RTraitLossScheduleSummary, ...]
    selection: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "domain_cell": self.domain_cell.as_dict(),
            "master_seeds": list(self.master_seeds),
            "records": [record.as_dict() for record in self.records],
            "schedules": [schedule.as_dict() for schedule in self.schedules],
            "selection": dict(self.selection),
        }

    def csv_rows(self) -> tuple[dict[str, object], ...]:
        base = self.domain_cell.as_dict()
        base["master_seeds"] = ",".join(str(seed) for seed in self.master_seeds)
        result = []
        for schedule in self.schedules:
            row = dict(base)
            row.update(schedule.as_dict())
            row["selection"] = json.dumps(dict(self.selection), sort_keys=True)
            result.append(row)
        return tuple(result)


def linear_normalized_barrier_schedule(
    *,
    anchor_barrier: float,
    canonical_interval_width: float,
    total_normalized_barrier_increase: float,
    horizon: int,
) -> tuple[float, ...]:
    """Build a strictly increasing post-baseline barrier schedule.

    Snapshot zero remains at the reconstructed interior anchor.  Entry zero of
    the returned sequence is used for generation one, after one equal ramp step.
    Thus the only exogenous change after projection is a monotone worsening of
    interaction support.
    """
    if horizon < 1:
        raise ValueError("horizon must be positive")
    if canonical_interval_width <= 0.0:
        raise ValueError("canonical_interval_width must be positive")
    if total_normalized_barrier_increase <= 0.0:
        raise ValueError("total_normalized_barrier_increase must be positive")
    total_increase = canonical_interval_width * total_normalized_barrier_increase
    return tuple(anchor_barrier + total_increase * generation / horizon for generation in range(1, horizon + 1))


def run_h2r_trait_loss_calibration(
    spec: ExperimentSpec,
    *,
    master_seeds: Sequence[int] = DEFAULT_H2R_CALIBRATION_MASTER_SEEDS,
    horizons: Sequence[int] = DEFAULT_CALIBRATION_HORIZONS,
    total_normalized_barrier_increases: Sequence[float] = DEFAULT_TOTAL_NORMALIZED_BARRIER_INCREASES,
    endpoint_padding_fraction: float = 0.5,
    stage_generations: int = 30,
    hold_generations: int = 30,
    nested_barrier_points: Sequence[int] = (25, 49, 97),
    interaction_separation_threshold: float = 0.05,
    maximum_normalized_bracket_width: float = 0.03,
) -> tuple[H2RTraitLossCalibrationCell, ...]:
    """Calibrate only post-baseline trait-loss availability in the primary domain."""
    seeds = _validate_master_seeds(master_seeds)
    candidate_horizons = _validate_horizons(horizons)
    candidate_increases = _validate_increases(total_normalized_barrier_increases)
    if spec.replicates < 1:
        raise ValueError("spec.replicates must be positive")

    output: list[H2RTraitLossCalibrationCell] = []
    for domain in primary_analysis_cells():
        records: list[H2RTraitLossCalibrationRecord] = []
        for master_seed in seeds:
            target_spec = _target_spec(spec, domain, master_seed)
            with _mutation_h1_runner(domain.mutation_rate):
                calibration = run_finite_h1_boundary_resolution_audit(
                    target_spec,
                    endpoint_padding_fraction=endpoint_padding_fraction,
                    stage_generations=stage_generations,
                    nested_barrier_points=nested_barrier_points,
                    interaction_separation_threshold=interaction_separation_threshold,
                    maximum_normalized_bracket_width=maximum_normalized_bracket_width,
                )
            if len(calibration) != 1:
                raise RuntimeError("targeted H2-R calibration must yield one H1 parameter cell")
            source_cell = calibration[0]
            isolated = _scenario_map(target_spec)[SCENARIO_EQUAL_ISOLATED]
            for source_record in source_cell.replicates:
                records.extend(
                    _records_from_source(
                        domain=domain,
                        spec=target_spec,
                        source_cell=source_cell,
                        source_record=source_record,
                        isolated=isolated,
                        master_seed=master_seed,
                        horizons=candidate_horizons,
                        increases=candidate_increases,
                        endpoint_padding_fraction=endpoint_padding_fraction,
                        stage_generations=stage_generations,
                        hold_generations=hold_generations,
                        interaction_separation_threshold=interaction_separation_threshold,
                    )
                )
        schedules = _summarise_schedules(records, seeds, candidate_horizons, candidate_increases)
        selection = _select_schedule(schedules)
        output.append(H2RTraitLossCalibrationCell(domain, seeds, tuple(records), schedules, selection))
    return tuple(output)


def _records_from_source(
    *,
    domain: MutationH1DomainCell,
    spec: ExperimentSpec,
    source_cell: FiniteH1BoundaryResolutionCell,
    source_record: FiniteH1BoundaryResolutionReplicate,
    isolated: LandscapeScenario,
    master_seed: int,
    horizons: Sequence[int],
    increases: Sequence[float],
    endpoint_padding_fraction: float,
    stage_generations: int,
    hold_generations: int,
    interaction_separation_threshold: float,
) -> tuple[H2RTraitLossCalibrationRecord, ...]:
    prepared = _prepare_mutation_high_state(
        domain.mutation_rate,
        spec,
        source_cell,
        source_record,
        endpoint_padding_fraction=endpoint_padding_fraction,
        stage_generations=stage_generations,
        hold_generations=hold_generations,
        interaction_separation_threshold=interaction_separation_threshold,
    )
    schedule_pairs = tuple((horizon, increase) for horizon in horizons for increase in increases)
    base = dict(
        mutation_rate=domain.mutation_rate,
        area_reference=domain.area_reference,
        interaction_feedback=domain.interaction_feedback,
        master_seed=master_seed,
        replicate_index=source_record.replicate_index,
        calibration_seed=source_record.seed,
        h1_resolution_supported=source_record.resolution_stable_h1_loop_mechanism_supported,
    )
    if prepared is None:
        return tuple(
            H2RTraitLossCalibrationRecord(
                **base,
                horizon=horizon,
                total_normalized_barrier_increase=increase,
                h1_full_state_source_prepared=False,
                anchor_barrier=None,
                canonical_interval_width=None,
                projection_supported=None,
                baseline_realised_high_trait_present=None,
                trajectory_seed=None,
                barrier_first_generation=None,
                barrier_final_generation=None,
                trait_loss_time_post_baseline=None,
            )
            for horizon, increase in schedule_pairs
        )

    source, anchor = prepared
    interval = source_cell.canonical_bistable_barrier_interval
    if interval is None or interval[1] <= interval[0]:
        raise RuntimeError("prepared H1 source requires a positive canonical bistable-interval width")
    width = interval[1] - interval[0]
    anchor_cell = _anchor_cell(source_cell.parameters, anchor)
    trajectory_seed = _trajectory_seed(source_record.seed)
    template = parameters_for_cell(spec, isolated, anchor_cell, seed=trajectory_seed)
    projected, invariants = project_full_state(source, template)
    if not invariants.projection_supported:
        return tuple(
            H2RTraitLossCalibrationRecord(
                **base,
                horizon=horizon,
                total_normalized_barrier_increase=increase,
                h1_full_state_source_prepared=True,
                anchor_barrier=anchor,
                canonical_interval_width=width,
                projection_supported=False,
                baseline_realised_high_trait_present=None,
                trajectory_seed=trajectory_seed,
                barrier_first_generation=None,
                barrier_final_generation=None,
                trait_loss_time_post_baseline=None,
            )
            for horizon, increase in schedule_pairs
        )

    records: list[H2RTraitLossCalibrationRecord] = []
    for horizon, increase in schedule_pairs:
        schedule = linear_normalized_barrier_schedule(
            anchor_barrier=anchor,
            canonical_interval_width=width,
            total_normalized_barrier_increase=increase,
            horizon=horizon,
        )
        result = simulate_with_symmetric_allele_mutation(
            _with_horizon(projected, horizon, trajectory_seed),
            mutation_rate=domain.mutation_rate,
            interaction_barrier_schedule=schedule,
        )
        baseline_present = any(
            item.realised_high_trait_occupied for item in result.snapshots[0].trait_occupancy
        )
        raw_loss_time = tau_trait_realised(result)
        post_baseline_loss = None if raw_loss_time in {None, 0} else raw_loss_time
        records.append(
            H2RTraitLossCalibrationRecord(
                **base,
                horizon=horizon,
                total_normalized_barrier_increase=increase,
                h1_full_state_source_prepared=True,
                anchor_barrier=anchor,
                canonical_interval_width=width,
                projection_supported=True,
                baseline_realised_high_trait_present=baseline_present,
                trajectory_seed=trajectory_seed,
                barrier_first_generation=schedule[0],
                barrier_final_generation=schedule[-1],
                trait_loss_time_post_baseline=post_baseline_loss,
            )
        )
    return tuple(records)


def _summarise_schedules(
    records: Sequence[H2RTraitLossCalibrationRecord],
    seeds: Sequence[int],
    horizons: Sequence[int],
    increases: Sequence[float],
) -> tuple[H2RTraitLossScheduleSummary, ...]:
    summaries: list[H2RTraitLossScheduleSummary] = []
    for horizon in horizons:
        for increase in increases:
            subset = tuple(
                record
                for record in records
                if record.horizon == horizon and record.total_normalized_barrier_increase == increase
            )
            seed_eligible = []
            seed_losses = []
            seed_probabilities: list[float | None] = []
            for seed in seeds:
                by_seed = tuple(record for record in subset if record.master_seed == seed)
                eligible = tuple(record for record in by_seed if record.eligible_for_trait_loss_denominator)
                losses = sum(record.trait_loss_observed_post_baseline is True for record in eligible)
                seed_eligible.append(len(eligible))
                seed_losses.append(losses)
                seed_probabilities.append(None if not eligible else losses / len(eligible))
            eligible_all = tuple(record for record in subset if record.eligible_for_trait_loss_denominator)
            summaries.append(
                H2RTraitLossScheduleSummary(
                    horizon=horizon,
                    total_normalized_barrier_increase=increase,
                    total_attempted_records=len(subset),
                    h1_full_state_source_prepared_count=sum(record.h1_full_state_source_prepared for record in subset),
                    projection_supported_count=sum(record.projection_supported is True for record in subset),
                    baseline_eligible_count=len(eligible_all),
                    trait_loss_count=sum(record.trait_loss_observed_post_baseline is True for record in eligible_all),
                    seed_block_baseline_eligible_counts=tuple(seed_eligible),
                    seed_block_trait_loss_counts=tuple(seed_losses),
                    seed_block_trait_loss_probabilities=tuple(seed_probabilities),
                )
            )
    return tuple(summaries)


def _select_schedule(summaries: Sequence[H2RTraitLossScheduleSummary]) -> dict[str, object]:
    selectable = tuple(
        summary.selectable_candidate
        for summary in summaries
        if summary.selectable_candidate is not None
    )
    if not selectable:
        return {
            "selected": None,
            "reason": "no candidate had a baseline-eligible H1-projected trajectory in every calibration seed block",
            "selection_uses_warning_outcomes": False,
            "target_trait_loss_probability_interval": list(TARGET_TRAIT_LOSS_PROBABILITY_INTERVAL),
        }
    selection = select_trait_loss_only_calibration(selectable)
    result = selection.as_dict()
    result["target_trait_loss_probability_interval"] = list(TARGET_TRAIT_LOSS_PROBABILITY_INTERVAL)
    return result


def write_h2r_trait_loss_calibration_artifacts(
    cells: Iterable[H2RTraitLossCalibrationCell],
    *,
    csv_path: str | Path,
    json_path: str | Path,
    manifest_path: str | Path,
) -> None:
    """Write the complete source, schedule, and trait-loss-only calibration ledger."""
    values = tuple(cells)
    if not values:
        raise ValueError("cells must be nonempty")
    csv_target, json_target, manifest_target = Path(csv_path), Path(json_path), Path(manifest_path)
    for target in (csv_target, json_target, manifest_target):
        target.parent.mkdir(parents=True, exist_ok=True)
    rows = [row for cell in values for row in cell.csv_rows()]
    with csv_target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted({key for row in rows for key in row}))
        writer.writeheader()
        writer.writerows(rows)
    with json_target.open("w", encoding="utf-8") as handle:
        json.dump([cell.as_dict() for cell in values], handle, indent=2, sort_keys=True)
    with manifest_target.open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "campaign": "h2r_trait_loss_only_calibration_v1",
                "h2r_protocol": h2r_protocol_manifest(),
                "primary_domain": domain_manifest(),
                "calibration_measure": "post-baseline realised high-trait loss only",
                "forbidden_calibration_outputs": [
                    "H_alpha",
                    "H_gamma",
                    "relative warning time",
                    "warning lead count",
                    "lead-time magnitude",
                ],
                "deterioration_schedule": {
                    "form": "barrier[g] = anchor + canonical_interval_width * normalized_increase * g / horizon",
                    "first_scheduled_generation": 1,
                    "baseline_is_anchor": True,
                },
                "cell_count": len(values),
                "schedule_rows": len(rows),
            },
            handle,
            indent=2,
            sort_keys=True,
        )


def _target_spec(spec: ExperimentSpec, domain: MutationH1DomainCell, master_seed: int) -> ExperimentSpec:
    return ExperimentSpec(
        experiment_id=spec.experiment_id,
        profile=spec.profile,
        total_area=spec.total_area,
        patch_count=spec.patch_count,
        generations=spec.generations,
        replicates=spec.replicates,
        master_seed=master_seed,
        area_reference_values=(domain.area_reference,),
        interaction_feedback_values=(domain.interaction_feedback,),
        interaction_barrier_values=(0.5,),
        migration_rate=spec.migration_rate,
        h_alpha_warning_threshold=spec.h_alpha_warning_threshold,
        h_gamma_warning_threshold=spec.h_gamma_warning_threshold,
        fst_warning_threshold=spec.fst_warning_threshold,
        allele_loss_threshold=spec.allele_loss_threshold,
        base_parameters=spec.base_parameters,
    )


def _mutation_h1_runner(mutation_rate: float):
    from causal_model.symmetric_allele_mutation_closure import patched_h1_mutation_runner

    return patched_h1_mutation_runner(mutation_rate)


def _anchor_cell(parameters: ParameterCell, anchor: float) -> ParameterCell:
    return ParameterCell(
        cell_index=parameters.cell_index,
        area_reference=parameters.area_reference,
        interaction_feedback=parameters.interaction_feedback,
        interaction_barrier=anchor,
    )


def _with_horizon(parameters, horizon: int, seed: int):
    from dataclasses import replace

    return replace(parameters, generations=horizon, random_seed=seed)


def _trajectory_seed(calibration_seed: int) -> int:
    # Common random numbers are deliberately shared among schedule candidates
    # from a given full-state source; no candidate gets a seed advantage.
    return (calibration_seed * 1_000_003 + 130_001) % (2**31 - 1)


def _validate_master_seeds(values: Sequence[int]) -> tuple[int, ...]:
    seeds = tuple(int(value) for value in values)
    if len(seeds) < 2 or len(seeds) != len(set(seeds)) or any(seed < 0 for seed in seeds):
        raise ValueError("master_seeds must contain at least two distinct non-negative values")
    return seeds


def _validate_horizons(values: Sequence[int]) -> tuple[int, ...]:
    result = tuple(int(value) for value in values)
    if not result or any(value < 1 for value in result) or len(result) != len(set(result)):
        raise ValueError("horizons must be distinct positive integers")
    return result


def _validate_increases(values: Sequence[float]) -> tuple[float, ...]:
    result = tuple(float(value) for value in values)
    if not result or any(value <= 0.0 for value in result) or len(result) != len(set(result)):
        raise ValueError("total_normalized_barrier_increases must be distinct positive values")
    return result
