"""Independent finite validation of the frozen H2-R relative-warning domain.

The preceding calibration selected one (mutation rate, A_ref, kappa, schedule)
configuration using realised trait loss only and independent calibration seeds.
This runner uses fresh master seeds and never reselects a cell or schedule.

For every attempted source it retains one of four states:
1. H1 full-state source unavailable;
2. projection unsupported;
3. trajectory available but a warning/loss comparison censored;
4. a valid same-replicate first-passage comparison.

For each available trajectory, all predeclared endpoints are reported:
H-alpha and H-gamma at relative declines r = 0.05, 0.10, 0.20.  No endpoint,
relative threshold, cell, or schedule is selected after looking at results.

Finite outputs are Type S evidence for the declared symmetric-mutation,
full-state-transfer, equal-isolated, ramp-and-hold closure.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Iterable, Mapping, Sequence

# Applies the existing keyword-only canonical certificate adapter to the
# imported H1 source-reconstruction helper without altering its logic.
import causal_model.mutation_primary_h1_h2_h3_runtime  # noqa: F401
from causal_model.finite_h1_boundary_resolution_audit import (
    FiniteH1BoundaryResolutionCell,
    FiniteH1BoundaryResolutionReplicate,
    run_finite_h1_boundary_resolution_audit,
)
from causal_model.finite_h1_fragment_projection_audit import _scenario_map, project_full_state
from causal_model.h2_relative_warning_contract import (
    DEFAULT_RELATIVE_DECLINE_FRACTIONS,
    RelativeWarningComparison,
    RelativeWarningDefinition,
    compare_relative_warning,
    h2r_protocol_manifest,
)
from causal_model.h2r_ramp_hold_trait_loss_calibration import (
    _anchor_cell,
    _target_spec,
    _trajectory_seed,
    _validate_master_seeds,
    ramp_and_hold_barrier_schedule,
)
from causal_model.h2r_validation_domain import (
    SELECTED_VALIDATION_DOMAIN,
    h2r_validation_domain_manifest,
)
from causal_model.multipatch_criticality_dynamics import tau_trait_realised
from causal_model.multipatch_criticality_experiments import (
    SCENARIO_EQUAL_ISOLATED,
    ExperimentSpec,
    LandscapeScenario,
    parameters_for_cell,
)
from causal_model.mutation_primary_h1_h2_h3_chain import _prepare_mutation_high_state
from causal_model.symmetric_allele_mutation_closure import (
    patched_h1_mutation_runner,
    simulate_with_symmetric_allele_mutation,
)

DEFAULT_H2R_VALIDATION_MASTER_SEEDS = (20261110, 20261111, 20261112, 20261113, 20261114)


def _definitions() -> tuple[RelativeWarningDefinition, ...]:
    return tuple(
        RelativeWarningDefinition(diversity_id, decline)
        for diversity_id in ("H_alpha", "H_gamma")
        for decline in DEFAULT_RELATIVE_DECLINE_FRACTIONS
    )


@dataclass(frozen=True)
class H2RValidationTrajectory:
    """A realised selected-schedule trajectory with raw diversity series."""

    trait_loss_time_post_baseline: int | None
    h_alpha_series: tuple[float, ...]
    h_gamma_series: tuple[float, ...]
    comparisons: tuple[RelativeWarningComparison, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "trait_loss_time_post_baseline": self.trait_loss_time_post_baseline,
            "h_alpha_series": list(self.h_alpha_series),
            "h_gamma_series": list(self.h_gamma_series),
            "comparisons": [comparison.as_dict() for comparison in self.comparisons],
        }


@dataclass(frozen=True)
class H2RValidationRecord:
    """One attempted fresh-seed source and, when possible, selected-schedule outcome."""

    master_seed: int
    replicate_index: int
    calibration_seed: int
    h1_resolution_supported: bool | None
    h1_full_state_source_prepared: bool
    anchor_barrier: float | None
    canonical_interval_width: float | None
    projection_supported: bool | None
    baseline_realised_high_trait_present: bool | None
    trajectory_seed: int | None
    barrier_first_generation: float | None
    barrier_at_hold: float | None
    outcome: H2RValidationTrajectory | None

    @property
    def trajectory_available(self) -> bool:
        return self.outcome is not None

    def as_dict(self) -> dict[str, object]:
        return {
            **asdict(self),
            "trajectory_available": self.trajectory_available,
            "outcome": None if self.outcome is None else self.outcome.as_dict(),
        }


@dataclass(frozen=True)
class H2RValidationResult:
    master_seeds: tuple[int, ...]
    records: tuple[H2RValidationRecord, ...]
    summary: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "validation_domain": SELECTED_VALIDATION_DOMAIN.as_dict(),
            "master_seeds": list(self.master_seeds),
            "records": [record.as_dict() for record in self.records],
            "summary": dict(self.summary),
        }


def run_h2r_independent_relative_validation(
    spec: ExperimentSpec,
    *,
    master_seeds: Sequence[int] = DEFAULT_H2R_VALIDATION_MASTER_SEEDS,
    endpoint_padding_fraction: float = 0.5,
    stage_generations: int = 30,
    hold_generations: int = 30,
    nested_barrier_points: Sequence[int] = (25, 49, 97),
    interaction_separation_threshold: float = 0.05,
    maximum_normalized_bracket_width: float = 0.03,
) -> H2RValidationResult:
    """Evaluate all predeclared relative-warning endpoints in fresh seeds."""
    seeds = _validate_master_seeds(master_seeds)
    if spec.replicates < 1:
        raise ValueError("spec.replicates must be positive")
    domain = SELECTED_VALIDATION_DOMAIN
    records: list[H2RValidationRecord] = []
    for master_seed in seeds:
        target_spec = _target_spec(spec, domain, master_seed)
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
            raise RuntimeError("targeted H2-R validation must yield exactly one H1 parameter cell")
        source_cell = calibration[0]
        isolated = _scenario_map(target_spec)[SCENARIO_EQUAL_ISOLATED]
        for source_record in source_cell.replicates:
            records.append(
                _run_record(
                    spec=target_spec,
                    source_cell=source_cell,
                    source_record=source_record,
                    isolated=isolated,
                    master_seed=master_seed,
                    endpoint_padding_fraction=endpoint_padding_fraction,
                    stage_generations=stage_generations,
                    hold_generations=hold_generations,
                    interaction_separation_threshold=interaction_separation_threshold,
                )
            )
    result_records = tuple(records)
    return H2RValidationResult(seeds, result_records, _summarise(result_records, seeds))


def _run_record(
    *,
    spec: ExperimentSpec,
    source_cell: FiniteH1BoundaryResolutionCell,
    source_record: FiniteH1BoundaryResolutionReplicate,
    isolated: LandscapeScenario,
    master_seed: int,
    endpoint_padding_fraction: float,
    stage_generations: int,
    hold_generations: int,
    interaction_separation_threshold: float,
) -> H2RValidationRecord:
    domain = SELECTED_VALIDATION_DOMAIN
    base = dict(
        master_seed=master_seed,
        replicate_index=source_record.replicate_index,
        calibration_seed=source_record.seed,
        h1_resolution_supported=source_record.resolution_stable_h1_loop_mechanism_supported,
    )
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
    if prepared is None:
        return H2RValidationRecord(
            **base,
            h1_full_state_source_prepared=False,
            anchor_barrier=None,
            canonical_interval_width=None,
            projection_supported=None,
            baseline_realised_high_trait_present=None,
            trajectory_seed=None,
            barrier_first_generation=None,
            barrier_at_hold=None,
            outcome=None,
        )

    source, anchor = prepared
    interval = source_cell.canonical_bistable_barrier_interval
    if interval is None or interval[1] <= interval[0]:
        raise RuntimeError("prepared H1 source requires a positive canonical interval width")
    width = interval[1] - interval[0]
    seed = _trajectory_seed(source_record.seed)
    template = parameters_for_cell(spec, isolated, _anchor_cell(source_cell.parameters, anchor), seed=seed)
    projected, invariants = project_full_state(source, template)
    if not invariants.projection_supported:
        return H2RValidationRecord(
            **base,
            h1_full_state_source_prepared=True,
            anchor_barrier=anchor,
            canonical_interval_width=width,
            projection_supported=False,
            baseline_realised_high_trait_present=None,
            trajectory_seed=seed,
            barrier_first_generation=None,
            barrier_at_hold=None,
            outcome=None,
        )

    barriers = ramp_and_hold_barrier_schedule(
        anchor_barrier=anchor,
        canonical_interval_width=width,
        schedule=domain.schedule,
    )
    result = simulate_with_symmetric_allele_mutation(
        replace(projected, generations=domain.schedule.total_generations, random_seed=seed),
        mutation_rate=domain.mutation_rate,
        interaction_barrier_schedule=barriers,
    )
    baseline_present = any(item.realised_high_trait_occupied for item in result.snapshots[0].trait_occupancy)
    trait_time = tau_trait_realised(result)
    trait_time = None if trait_time in {None, 0} else trait_time
    alpha = tuple(snapshot.h_alpha for snapshot in result.snapshots)
    gamma = tuple(snapshot.h_gamma for snapshot in result.snapshots)
    comparisons = tuple(
        compare_relative_warning(
            alpha if definition.diversity_id == "H_alpha" else gamma,
            trait_loss_time=trait_time,
            definition=definition,
        )
        for definition in _definitions()
    )
    return H2RValidationRecord(
        **base,
        h1_full_state_source_prepared=True,
        anchor_barrier=anchor,
        canonical_interval_width=width,
        projection_supported=True,
        baseline_realised_high_trait_present=baseline_present,
        trajectory_seed=seed,
        barrier_first_generation=barriers[0],
        barrier_at_hold=barriers[-1],
        outcome=H2RValidationTrajectory(trait_time, alpha, gamma, comparisons),
    )


def _summarise(records: Sequence[H2RValidationRecord], seeds: Sequence[int]) -> dict[str, object]:
    definitions = _definitions()
    source_prepared = tuple(record for record in records if record.h1_full_state_source_prepared)
    projection = tuple(record for record in source_prepared if record.projection_supported is True)
    trajectories = tuple(record for record in projection if record.outcome is not None)
    endpoint_summaries = [_summarise_endpoint(records, seeds, definition) for definition in definitions]
    return {
        "denominators": {
            "attempted_seed_replicates": len(records),
            "h1_full_state_source_prepared_count": len(source_prepared),
            "projection_supported_count": len(projection),
            "trajectory_available_count": len(trajectories),
            "trait_loss_observed_count": sum(
                record.outcome is not None and record.outcome.trait_loss_time_post_baseline is not None
                for record in trajectories
            ),
        },
        "endpoint_summaries": endpoint_summaries,
        "interpretation": {
            "selection_repeated": False,
            "cell_or_schedule_search_after_warning": False,
            "endpoint_family": [
                {"diversity_id": definition.diversity_id, "relative_decline_fraction": definition.relative_decline_fraction}
                for definition in definitions
            ],
            "censoring_rule": "records without an observed warning or post-baseline trait-loss event are retained as censored",
            "evidence_label": "Type S finite evidence only",
        },
    }


def _summarise_endpoint(
    records: Sequence[H2RValidationRecord],
    seeds: Sequence[int],
    definition: RelativeWarningDefinition,
) -> dict[str, object]:
    comparisons = [
        _find_comparison(record, definition)
        for record in records
        if record.outcome is not None
    ]
    comparisons = [comparison for comparison in comparisons if comparison is not None]
    valid = [comparison for comparison in comparisons if comparison.valid_pair]
    leads = [comparison for comparison in valid if comparison.warning_leads is True]
    ties = [comparison for comparison in valid if comparison.warning_leads is False and comparison.lead_time_trait_minus_warning == 0]
    lags = [comparison for comparison in valid if comparison.warning_leads is False and comparison.lead_time_trait_minus_warning is not None and comparison.lead_time_trait_minus_warning < 0]
    seed_rows = []
    for seed in seeds:
        by_seed = [
            _find_comparison(record, definition)
            for record in records
            if record.master_seed == seed and record.outcome is not None
        ]
        by_seed = [comparison for comparison in by_seed if comparison is not None]
        valid_seed = [comparison for comparison in by_seed if comparison.valid_pair]
        leads_seed = sum(comparison.warning_leads is True for comparison in valid_seed)
        seed_rows.append({
            "master_seed": seed,
            "trajectory_available_count": len(by_seed),
            "valid_pair_count": len(valid_seed),
            "warning_lead_count": leads_seed,
            "warning_lead_probability": None if not valid_seed else leads_seed / len(valid_seed),
        })
    return {
        "definition": asdict(definition),
        "trajectory_available_count": len(comparisons),
        "baseline_eligible_count": sum(comparison.baseline_eligible for comparison in comparisons),
        "warning_observed_count": sum(comparison.warning_time is not None for comparison in comparisons),
        "trait_loss_observed_count": sum(comparison.trait_loss_time is not None for comparison in comparisons),
        "valid_pair_count": len(valid),
        "censored_count": sum(comparison.censored for comparison in comparisons),
        "warning_lead_count": len(leads),
        "warning_tie_count": len(ties),
        "warning_lag_count": len(lags),
        "warning_lead_probability_among_valid_pairs": None if not valid else len(leads) / len(valid),
        "seed_blocks": seed_rows,
    }


def _find_comparison(
    record: H2RValidationRecord,
    definition: RelativeWarningDefinition,
) -> RelativeWarningComparison | None:
    if record.outcome is None:
        return None
    for comparison in record.outcome.comparisons:
        if comparison.definition == definition:
            return comparison
    raise RuntimeError("trajectory lacks a predeclared relative-warning definition")


def write_h2r_independent_validation_artifacts(
    result: H2RValidationResult,
    *,
    csv_path: str | Path,
    json_path: str | Path,
    manifest_path: str | Path,
) -> None:
    """Write raw trajectories, comparison rows, and the locked validation manifest."""
    csv_target, json_target, manifest_target = Path(csv_path), Path(json_path), Path(manifest_path)
    for target in (csv_target, json_target, manifest_target):
        target.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for record in result.records:
        if record.outcome is None:
            rows.append(_csv_row(record, None))
        else:
            for comparison in record.outcome.comparisons:
                rows.append(_csv_row(record, comparison))
    with csv_target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted({key for row in rows for key in row}))
        writer.writeheader()
        writer.writerows(rows)
    with json_target.open("w", encoding="utf-8") as handle:
        json.dump(result.as_dict(), handle, indent=2, sort_keys=True)
    with manifest_target.open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "campaign": "h2r_independent_relative_warning_validation_v1",
                "h2r_protocol": h2r_protocol_manifest(),
                "validation_domain": h2r_validation_domain_manifest(),
                "outcome_definitions": {
                    "trait_loss": "post-baseline realised high-trait loss",
                    "relative_warning": "H_x(t) <= (1-r) H_x(0), t > 0",
                    "relative_decline_fractions": list(DEFAULT_RELATIVE_DECLINE_FRACTIONS),
                    "diversity_metrics": ["H_alpha", "H_gamma"],
                    "censoring": "missing warning or trait-loss event yields no ordering claim",
                },
                "selection_repeated": False,
                "finite_evidence_label": "Type S",
                "master_seeds": list(result.master_seeds),
                "summary": result.summary,
            },
            handle,
            indent=2,
            sort_keys=True,
        )


def _csv_row(record: H2RValidationRecord, comparison: RelativeWarningComparison | None) -> dict[str, object]:
    domain = SELECTED_VALIDATION_DOMAIN
    base = {
        "mutation_rate": domain.mutation_rate,
        "area_reference": domain.area_reference,
        "interaction_feedback": domain.interaction_feedback,
        "schedule_id": domain.schedule.schedule_id,
        "master_seed": record.master_seed,
        "replicate_index": record.replicate_index,
        "calibration_seed": record.calibration_seed,
        "h1_resolution_supported": record.h1_resolution_supported,
        "h1_full_state_source_prepared": record.h1_full_state_source_prepared,
        "projection_supported": record.projection_supported,
        "baseline_realised_high_trait_present": record.baseline_realised_high_trait_present,
        "trajectory_available": record.trajectory_available,
        "trait_loss_time_post_baseline": None if record.outcome is None else record.outcome.trait_loss_time_post_baseline,
    }
    if comparison is None:
        return base
    base.update({
        "diversity_id": comparison.definition.diversity_id,
        "relative_decline_fraction": comparison.definition.relative_decline_fraction,
        "baseline": comparison.baseline,
        "warning_threshold": comparison.warning_threshold,
        "warning_time": comparison.warning_time,
        "trait_loss_time": comparison.trait_loss_time,
        "baseline_eligible": comparison.baseline_eligible,
        "valid_pair": comparison.valid_pair,
        "censored": comparison.censored,
        "warning_leads": comparison.warning_leads,
        "lead_time_trait_minus_warning": comparison.lead_time_trait_minus_warning,
    })
    return base
