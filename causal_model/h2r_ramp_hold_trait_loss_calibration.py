"""Ramp-and-hold trait-loss-only calibration for the conditional H2-R hypothesis.

Calibration v1 used a linear barrier ramp over its complete 60- or 120-generation
observation window.  It selected no schedule: all 12 frozen primary cells failed
the all-seed-block trait-loss availability rule.  This v2 runner keeps the same
H1 source reconstruction, mutation closure, equal-isolated projection, trait-loss
endpoint, selection rule, and forbidden genetic-warning outputs.  It changes only
the *time shape* of exogenous interaction-support deterioration:

    ramp from the H1 interior anchor for 30 generations
    then hold the reached barrier for 90 or 210 generations.

The total normalized barrier increase remains in the original {0.15, 0.30, 0.45}
family.  Therefore v2 asks whether the same declared deterioration magnitude
needs time at its final support level for realised trait loss to be observable.

This is a calibration runner, not an H2-R validation runner.  It intentionally
never computes H-alpha, H-gamma, relative-warning times, warning leads, or
lead-time magnitudes.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Iterable, Mapping, Sequence

# Import first so the chain module receives its keyword-only certificate adapter.
import causal_model.mutation_primary_h1_h2_h3_runtime  # noqa: F401
from causal_model.finite_h1_boundary_resolution_audit import (
    FiniteH1BoundaryResolutionCell,
    FiniteH1BoundaryResolutionReplicate,
    run_finite_h1_boundary_resolution_audit,
)
from causal_model.finite_h1_fragment_projection_audit import _scenario_map, project_full_state
from causal_model.h2_relative_warning_contract import (
    DeteriorationCalibrationCandidate,
    TARGET_TRAIT_LOSS_PROBABILITY_INTERVAL,
    h2r_protocol_manifest,
    select_trait_loss_only_calibration,
)
from causal_model.h2r_trait_loss_calibration import (
    _anchor_cell,
    _target_spec,
    _trajectory_seed,
    _validate_master_seeds,
)
from causal_model.multipatch_criticality_dynamics import tau_trait_realised
from causal_model.multipatch_criticality_experiments import (
    SCENARIO_EQUAL_ISOLATED,
    ExperimentSpec,
    LandscapeScenario,
    parameters_for_cell,
)
from causal_model.mutation_h1_primary_domain import (
    MutationH1DomainCell,
    domain_manifest,
    primary_analysis_cells,
)
from causal_model.mutation_primary_h1_h2_h3_chain import _prepare_mutation_high_state
from causal_model.symmetric_allele_mutation_closure import (
    patched_h1_mutation_runner,
    simulate_with_symmetric_allele_mutation,
)

PREVIOUS_LINEAR_CALIBRATION_RUN_ID = 28493522149
DEFAULT_H2R_RAMP_HOLD_MASTER_SEEDS = (20261010, 20261011, 20261012, 20261013, 20261014)
DEFAULT_RAMP_GENERATIONS = 30
DEFAULT_HOLD_GENERATIONS = (90, 210)
DEFAULT_TOTAL_NORMALIZED_BARRIER_INCREASES = (0.15, 0.30, 0.45)


@dataclass(frozen=True, order=True)
class RampHoldSchedule:
    """One predeclared monotone external interaction-support deterioration path."""

    ramp_generations: int
    hold_generations: int
    total_normalized_barrier_increase: float

    def __post_init__(self) -> None:
        if self.ramp_generations < 1:
            raise ValueError("ramp_generations must be positive")
        if self.hold_generations < 1:
            raise ValueError("hold_generations must be positive")
        if self.total_normalized_barrier_increase <= 0.0:
            raise ValueError("total_normalized_barrier_increase must be positive")

    @property
    def total_generations(self) -> int:
        return self.ramp_generations + self.hold_generations

    @property
    def schedule_id(self) -> str:
        increase = str(self.total_normalized_barrier_increase).replace(".", "p")
        return f"ramp{self.ramp_generations}_hold{self.hold_generations}_d{increase}"

    def as_dict(self) -> dict[str, object]:
        return {
            "schedule_id": self.schedule_id,
            "ramp_generations": self.ramp_generations,
            "hold_generations": self.hold_generations,
            "total_generations": self.total_generations,
            "total_normalized_barrier_increase": self.total_normalized_barrier_increase,
        }


DEFAULT_RAMP_HOLD_SCHEDULES = tuple(
    RampHoldSchedule(DEFAULT_RAMP_GENERATIONS, hold, increase)
    for hold in DEFAULT_HOLD_GENERATIONS
    for increase in DEFAULT_TOTAL_NORMALIZED_BARRIER_INCREASES
)


@dataclass(frozen=True)
class RampHoldCalibrationRecord:
    """One attempted trajectory; genetic-warning fields are deliberately absent."""

    mutation_rate: float
    area_reference: float
    interaction_feedback: float
    master_seed: int
    replicate_index: int
    calibration_seed: int
    schedule: RampHoldSchedule
    h1_resolution_supported: bool | None
    h1_full_state_source_prepared: bool
    anchor_barrier: float | None
    canonical_interval_width: float | None
    projection_supported: bool | None
    baseline_realised_high_trait_present: bool | None
    trajectory_seed: int | None
    barrier_first_generation: float | None
    barrier_at_hold: float | None
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
        result["schedule"] = self.schedule.as_dict()
        result["eligible_for_trait_loss_denominator"] = self.eligible_for_trait_loss_denominator
        result["trait_loss_observed_post_baseline"] = self.trait_loss_observed_post_baseline
        return result


@dataclass(frozen=True)
class RampHoldScheduleSummary:
    """Trait-loss-only calibration evidence for one ramp-and-hold candidate."""

    schedule: RampHoldSchedule
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
    def selection_candidate(self) -> DeteriorationCalibrationCandidate | None:
        if any(value is None for value in self.seed_block_trait_loss_probabilities):
            return None
        return DeteriorationCalibrationCandidate(
            horizon=self.schedule.total_generations,
            total_normalized_barrier_increase=self.schedule.total_normalized_barrier_increase,
            trait_loss_probability_by_seed_block=tuple(
                float(value) for value in self.seed_block_trait_loss_probabilities
            ),
        )

    def as_dict(self) -> dict[str, object]:
        result = asdict(self)
        result["schedule"] = self.schedule.as_dict()
        result["pooled_trait_loss_probability"] = self.pooled_trait_loss_probability
        result["selection_candidate_available"] = self.selection_candidate is not None
        result["selection_uses_warning_outcomes"] = False
        return result


@dataclass(frozen=True)
class RampHoldCalibrationCell:
    """Full v2 calibration ledger for one frozen primary mutation-H1 cell."""

    domain_cell: MutationH1DomainCell
    master_seeds: tuple[int, ...]
    records: tuple[RampHoldCalibrationRecord, ...]
    schedules: tuple[RampHoldScheduleSummary, ...]
    selection: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "domain_cell": self.domain_cell.as_dict(),
            "master_seeds": list(self.master_seeds),
            "records": [record.as_dict() for record in self.records],
            "schedules": [summary.as_dict() for summary in self.schedules],
            "selection": dict(self.selection),
        }

    def csv_rows(self) -> tuple[dict[str, object], ...]:
        base = self.domain_cell.as_dict()
        base["master_seeds"] = ",".join(str(seed) for seed in self.master_seeds)
        rows: list[dict[str, object]] = []
        for summary in self.schedules:
            row = dict(base)
            row.update(summary.as_dict())
            row["selection"] = json.dumps(dict(self.selection), sort_keys=True)
            rows.append(row)
        return tuple(rows)


def ramp_and_hold_barrier_schedule(
    *,
    anchor_barrier: float,
    canonical_interval_width: float,
    schedule: RampHoldSchedule,
) -> tuple[float, ...]:
    """Return a nondecreasing post-baseline barrier sequence.

    Baseline snapshot zero remains exactly at ``anchor_barrier``.  The first
    simulated generation advances one ramp increment.  The last ramp value is
    then repeated unchanged for the declared positive hold duration.
    """
    if canonical_interval_width <= 0.0:
        raise ValueError("canonical_interval_width must be positive")
    total = canonical_interval_width * schedule.total_normalized_barrier_increase
    ramp = tuple(
        anchor_barrier + total * generation / schedule.ramp_generations
        for generation in range(1, schedule.ramp_generations + 1)
    )
    return ramp + (ramp[-1],) * schedule.hold_generations


def run_h2r_ramp_hold_trait_loss_calibration(
    spec: ExperimentSpec,
    *,
    master_seeds: Sequence[int] = DEFAULT_H2R_RAMP_HOLD_MASTER_SEEDS,
    schedules: Sequence[RampHoldSchedule] = DEFAULT_RAMP_HOLD_SCHEDULES,
    endpoint_padding_fraction: float = 0.5,
    stage_generations: int = 30,
    hold_generations: int = 30,
    nested_barrier_points: Sequence[int] = (25, 49, 97),
    interaction_separation_threshold: float = 0.05,
    maximum_normalized_bracket_width: float = 0.03,
) -> tuple[RampHoldCalibrationCell, ...]:
    """Run v2 schedule calibration using realised trait loss only."""
    seeds = _validate_master_seeds(master_seeds)
    candidate_schedules = _validate_schedules(schedules)
    if spec.replicates < 1:
        raise ValueError("spec.replicates must be positive")

    output: list[RampHoldCalibrationCell] = []
    for domain in primary_analysis_cells():
        records: list[RampHoldCalibrationRecord] = []
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
                raise RuntimeError("targeted H2-R ramp-hold calibration must yield one H1 parameter cell")
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
                        schedules=candidate_schedules,
                        endpoint_padding_fraction=endpoint_padding_fraction,
                        stage_generations=stage_generations,
                        hold_generations=hold_generations,
                        interaction_separation_threshold=interaction_separation_threshold,
                    )
                )
        summaries = _summarise_schedules(records, seeds, candidate_schedules)
        output.append(
            RampHoldCalibrationCell(
                domain_cell=domain,
                master_seeds=seeds,
                records=tuple(records),
                schedules=summaries,
                selection=_select_schedule(summaries),
            )
        )
    return tuple(output)


def _records_from_source(
    *,
    domain: MutationH1DomainCell,
    spec: ExperimentSpec,
    source_cell: FiniteH1BoundaryResolutionCell,
    source_record: FiniteH1BoundaryResolutionReplicate,
    isolated: LandscapeScenario,
    master_seed: int,
    schedules: Sequence[RampHoldSchedule],
    endpoint_padding_fraction: float,
    stage_generations: int,
    hold_generations: int,
    interaction_separation_threshold: float,
) -> tuple[RampHoldCalibrationRecord, ...]:
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
            RampHoldCalibrationRecord(
                **base,
                schedule=schedule,
                h1_full_state_source_prepared=False,
                anchor_barrier=None,
                canonical_interval_width=None,
                projection_supported=None,
                baseline_realised_high_trait_present=None,
                trajectory_seed=None,
                barrier_first_generation=None,
                barrier_at_hold=None,
                trait_loss_time_post_baseline=None,
            )
            for schedule in schedules
        )

    source, anchor = prepared
    interval = source_cell.canonical_bistable_barrier_interval
    if interval is None or interval[1] <= interval[0]:
        raise RuntimeError("prepared H1 source requires a positive canonical bistable interval width")
    interval_width = interval[1] - interval[0]
    trajectory_seed = _trajectory_seed(source_record.seed)
    template = parameters_for_cell(
        spec,
        isolated,
        _anchor_cell(source_cell.parameters, anchor),
        seed=trajectory_seed,
    )
    projected, invariants = project_full_state(source, template)
    if not invariants.projection_supported:
        return tuple(
            RampHoldCalibrationRecord(
                **base,
                schedule=schedule,
                h1_full_state_source_prepared=True,
                anchor_barrier=anchor,
                canonical_interval_width=interval_width,
                projection_supported=False,
                baseline_realised_high_trait_present=None,
                trajectory_seed=trajectory_seed,
                barrier_first_generation=None,
                barrier_at_hold=None,
                trait_loss_time_post_baseline=None,
            )
            for schedule in schedules
        )

    records: list[RampHoldCalibrationRecord] = []
    for schedule in schedules:
        barriers = ramp_and_hold_barrier_schedule(
            anchor_barrier=anchor,
            canonical_interval_width=interval_width,
            schedule=schedule,
        )
        result = simulate_with_symmetric_allele_mutation(
            replace(projected, generations=schedule.total_generations, random_seed=trajectory_seed),
            mutation_rate=domain.mutation_rate,
            interaction_barrier_schedule=barriers,
        )
        baseline_present = any(item.realised_high_trait_occupied for item in result.snapshots[0].trait_occupancy)
        raw_loss_time = tau_trait_realised(result)
        records.append(
            RampHoldCalibrationRecord(
                **base,
                schedule=schedule,
                h1_full_state_source_prepared=True,
                anchor_barrier=anchor,
                canonical_interval_width=interval_width,
                projection_supported=True,
                baseline_realised_high_trait_present=baseline_present,
                trajectory_seed=trajectory_seed,
                barrier_first_generation=barriers[0],
                barrier_at_hold=barriers[-1],
                trait_loss_time_post_baseline=None if raw_loss_time in {None, 0} else raw_loss_time,
            )
        )
    return tuple(records)


def _summarise_schedules(
    records: Sequence[RampHoldCalibrationRecord],
    seeds: Sequence[int],
    schedules: Sequence[RampHoldSchedule],
) -> tuple[RampHoldScheduleSummary, ...]:
    summaries: list[RampHoldScheduleSummary] = []
    for schedule in schedules:
        subset = tuple(record for record in records if record.schedule == schedule)
        seed_eligible: list[int] = []
        seed_losses: list[int] = []
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
            RampHoldScheduleSummary(
                schedule=schedule,
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


def _select_schedule(summaries: Sequence[RampHoldScheduleSummary]) -> dict[str, object]:
    lookup: dict[tuple[int, float], RampHoldScheduleSummary] = {}
    candidates: list[DeteriorationCalibrationCandidate] = []
    for summary in summaries:
        candidate = summary.selection_candidate
        if candidate is None:
            continue
        key = (candidate.horizon, candidate.total_normalized_barrier_increase)
        if key in lookup:
            raise RuntimeError("ramp-hold schedule identity must be unique by total horizon and normalized increase")
        lookup[key] = summary
        candidates.append(candidate)
    if not candidates:
        return {
            "selected": None,
            "reason": "no candidate had a baseline-eligible H1-projected trajectory in every calibration seed block",
            "selection_uses_warning_outcomes": False,
            "target_trait_loss_probability_interval": list(TARGET_TRAIT_LOSS_PROBABILITY_INTERVAL),
        }
    selection = select_trait_loss_only_calibration(candidates)
    result = selection.as_dict()
    if selection.selected is not None:
        key = (selection.selected.horizon, selection.selected.total_normalized_barrier_increase)
        result["selected_schedule"] = lookup[key].schedule.as_dict()
    else:
        result["selected_schedule"] = None
    result["target_trait_loss_probability_interval"] = list(TARGET_TRAIT_LOSS_PROBABILITY_INTERVAL)
    return result


def write_h2r_ramp_hold_calibration_artifacts(
    cells: Iterable[RampHoldCalibrationCell],
    *,
    csv_path: str | Path,
    json_path: str | Path,
    manifest_path: str | Path,
) -> None:
    """Write complete v2 calibration ledgers and a protocol manifest."""
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
                "campaign": "h2r_ramp_hold_trait_loss_only_calibration_v2",
                "predecessor_linear_calibration_run_id": PREVIOUS_LINEAR_CALIBRATION_RUN_ID,
                "h2r_protocol": h2r_protocol_manifest(),
                "primary_domain": domain_manifest(),
                "calibration_measure": "post-baseline realised high-trait loss only",
                "forbidden_calibration_outputs": [
                    "H_alpha", "H_gamma", "relative warning time", "warning lead count", "lead-time magnitude"
                ],
                "schedule_family": {
                    "form": "ramp barrier from anchor over 30 generations, then hold final barrier",
                    "baseline_is_anchor": True,
                    "ramp_generations": DEFAULT_RAMP_GENERATIONS,
                    "hold_generations": list(DEFAULT_HOLD_GENERATIONS),
                    "total_normalized_barrier_increases": list(DEFAULT_TOTAL_NORMALIZED_BARRIER_INCREASES),
                    "candidates": [schedule.as_dict() for schedule in DEFAULT_RAMP_HOLD_SCHEDULES],
                },
                "cell_count": len(values),
                "schedule_rows": len(rows),
            },
            handle,
            indent=2,
            sort_keys=True,
        )


def _validate_schedules(values: Sequence[RampHoldSchedule]) -> tuple[RampHoldSchedule, ...]:
    schedules = tuple(values)
    if not schedules:
        raise ValueError("schedules must be nonempty")
    if len(schedules) != len(set(schedules)):
        raise ValueError("schedules must be unique")
    keys = {(schedule.total_generations, schedule.total_normalized_barrier_increase) for schedule in schedules}
    if len(keys) != len(schedules):
        raise ValueError("schedule total generation and normalized increase pairs must be unique")
    return schedules
