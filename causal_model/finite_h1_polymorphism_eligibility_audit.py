"""Genetic-eligibility audit for H1 full-state high branches.

A branch-conditioned H2 early-warning claim requires genetic variation at the
moment fragmentation begins.  A warning threshold crossed at generation zero is
not an early warning of the subsequent fragmentation trajectory: it is a
pre-existing state.  Likewise, an H3 genetic-diversity contrast cannot be
estimated from a source already fixed for one allele.

This Type S audit therefore sits between conservation-preserving H1 state
projection and any H2/H3 dynamic campaign.  It reuses the declared H1
full-state high-source and fragment projection construction, then classifies
whether every resulting baseline has:

* polymorphism (epsilon < p < 1 - epsilon),
* positive H-alpha and H-gamma, and
* neither H-alpha nor H-gamma already below its predeclared warning threshold.

Records that fail eligibility are retained as *baseline-ineligible*, not
recoded as failed genetic warnings or failed fragmentation effects.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from causal_model.finite_h1_fragment_projection_audit import (
    DEFAULT_MASTER_SEEDS,
    FragmentProjectionCell,
    FragmentProjectionReplicate,
    run_finite_h1_fragment_projection_audit,
)
from causal_model.multipatch_criticality_experiments import (
    SCENARIO_EQUAL_ISOLATED,
    SCENARIO_EQUAL_MIGRATING,
    SCENARIO_ONE_LARGE,
    ExperimentSpec,
    ParameterCell,
)

SCENARIO_IDS = (SCENARIO_ONE_LARGE, SCENARIO_EQUAL_ISOLATED, SCENARIO_EQUAL_MIGRATING)


@dataclass(frozen=True)
class GeneticBaseline:
    """Genetic state and eligibility of a projected landscape at generation zero."""

    scenario_id: str
    high_allele_frequency_mean: float
    h_alpha: float
    h_gamma: float
    fst: float | None
    polymorphic: bool
    h_alpha_warning_preexisting: bool
    h_gamma_warning_preexisting: bool
    h2_dynamic_warning_eligible: bool
    h3_genetic_contrast_eligible: bool

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class PolymorphismEligibilityReplicate:
    """One H1 high-source replicate with baseline eligibility by landscape."""

    master_seed: int
    replicate_index: int
    calibration_seed: int
    h1_full_state_prepared: bool | None
    anchor_barrier: float | None
    baselines: Mapping[str, GeneticBaseline] | None

    def as_dict(self) -> dict[str, object]:
        return {
            "master_seed": self.master_seed,
            "replicate_index": self.replicate_index,
            "calibration_seed": self.calibration_seed,
            "h1_full_state_prepared": self.h1_full_state_prepared,
            "anchor_barrier": self.anchor_barrier,
            "baselines": None if self.baselines is None else {key: value.as_dict() for key, value in self.baselines.items()},
        }


@dataclass(frozen=True)
class PolymorphismEligibilityCell:
    """H1-conditioned genetic eligibility across independent master seeds."""

    experiment_id: str
    profile: str
    pair_index: int
    parameters: ParameterCell
    master_seeds: tuple[int, ...]
    polymorphism_epsilon: float
    h_alpha_warning_threshold: float
    h_gamma_warning_threshold: float
    projection_design: Mapping[str, object]
    replicates: tuple[PolymorphismEligibilityReplicate, ...]
    summary: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "pair_index": self.pair_index,
            "parameters": asdict(self.parameters),
            "master_seeds": list(self.master_seeds),
            "eligibility_rule": {
                "polymorphism_epsilon": self.polymorphism_epsilon,
                "h_alpha_warning_threshold": self.h_alpha_warning_threshold,
                "h_gamma_warning_threshold": self.h_gamma_warning_threshold,
                "h2_dynamic_warning_eligible": "polymorphic and neither H-alpha nor H-gamma warning is already present at generation zero",
                "h3_genetic_contrast_eligible": "polymorphic with positive H-alpha and H-gamma at generation zero",
            },
            "projection_design": dict(self.projection_design),
            "replicates": [record.as_dict() for record in self.replicates],
            "summary": dict(self.summary),
        }

    def to_csv_row(self) -> dict[str, object]:
        row: dict[str, object] = {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "pair_index": self.pair_index,
            "master_seeds": ",".join(str(value) for value in self.master_seeds),
            "replicate_count": len(self.replicates),
            "polymorphism_epsilon": self.polymorphism_epsilon,
            "h_alpha_warning_threshold": self.h_alpha_warning_threshold,
            "h_gamma_warning_threshold": self.h_gamma_warning_threshold,
        }
        row.update(asdict(self.parameters))
        row.update(_flatten(self.summary))
        return row


def run_finite_h1_polymorphism_eligibility_audit(
    spec: ExperimentSpec,
    *,
    master_seeds: Sequence[int] = DEFAULT_MASTER_SEEDS,
    endpoint_padding_fraction: float = 0.5,
    stage_generations: int = 30,
    hold_generations: int = 30,
    nested_barrier_points: Sequence[int] = (25, 49, 97),
    interaction_separation_threshold: float = 0.05,
    maximum_normalized_bracket_width: float = 0.03,
    polymorphism_epsilon: float = 1e-12,
) -> tuple[PolymorphismEligibilityCell, ...]:
    """Classify whether H1 high full states can support genetic H2/H3 endpoints."""
    if not 0.0 < polymorphism_epsilon < 0.5:
        raise ValueError("polymorphism_epsilon must lie strictly between 0 and 0.5")
    projection_cells = run_finite_h1_fragment_projection_audit(
        spec,
        master_seeds=master_seeds,
        endpoint_padding_fraction=endpoint_padding_fraction,
        stage_generations=stage_generations,
        hold_generations=hold_generations,
        nested_barrier_points=nested_barrier_points,
        interaction_separation_threshold=interaction_separation_threshold,
        maximum_normalized_bracket_width=maximum_normalized_bracket_width,
    )
    output: list[PolymorphismEligibilityCell] = []
    for source_cell in projection_cells:
        records = tuple(
            _classify_record(record, spec, polymorphism_epsilon)
            for record in source_cell.replicates
        )
        output.append(
            PolymorphismEligibilityCell(
                experiment_id=spec.experiment_id,
                profile=spec.profile,
                pair_index=source_cell.pair_index,
                parameters=source_cell.parameters,
                master_seeds=source_cell.master_seeds,
                polymorphism_epsilon=polymorphism_epsilon,
                h_alpha_warning_threshold=spec.h_alpha_warning_threshold,
                h_gamma_warning_threshold=spec.h_gamma_warning_threshold,
                projection_design={
                    "endpoint_padding_fraction": source_cell.endpoint_padding_fraction,
                    "stage_generations": source_cell.stage_generations,
                    "hold_generations": source_cell.hold_generations,
                    "nested_barrier_points": list(source_cell.nested_barrier_points),
                    "interaction_separation_threshold": source_cell.interaction_separation_threshold,
                    "projection_source": "finite_h1_fragment_projection_v1",
                },
                replicates=records,
                summary=_summarise(records, source_cell.master_seeds),
            )
        )
    return tuple(output)


def write_finite_h1_polymorphism_eligibility_artifacts(
    cells: Iterable[PolymorphismEligibilityCell],
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


def _classify_record(
    record: FragmentProjectionReplicate,
    spec: ExperimentSpec,
    epsilon: float,
) -> PolymorphismEligibilityReplicate:
    if record.h1_full_state_hold_supported is not True or record.scenario_projections is None:
        return PolymorphismEligibilityReplicate(
            master_seed=record.master_seed,
            replicate_index=record.replicate_index,
            calibration_seed=record.calibration_seed,
            h1_full_state_prepared=record.h1_full_state_hold_supported,
            anchor_barrier=record.anchor_barrier,
            baselines=None,
        )
    baselines = {
        scenario_id: _genetic_baseline(
            scenario_id,
            record.scenario_projections[scenario_id].target_initial_population,
            record.scenario_projections[scenario_id].target_initial_high_allele_frequency,
            spec,
            epsilon,
        )
        for scenario_id in SCENARIO_IDS
    }
    return PolymorphismEligibilityReplicate(
        master_seed=record.master_seed,
        replicate_index=record.replicate_index,
        calibration_seed=record.calibration_seed,
        h1_full_state_prepared=True,
        anchor_barrier=record.anchor_barrier,
        baselines=baselines,
    )


def _genetic_baseline(
    scenario_id: str,
    population: Sequence[int],
    frequencies: Sequence[float],
    spec: ExperimentSpec,
    epsilon: float,
) -> GeneticBaseline:
    if not population or len(population) != len(frequencies):
        raise ValueError("population and frequencies must be nonempty and aligned")
    total = sum(population)
    if total <= 0:
        raise ValueError("population total must be positive")
    weights = tuple(value / total for value in population)
    p_mean = sum(weight * value for weight, value in zip(weights, frequencies))
    h_alpha = sum(weight * 2.0 * value * (1.0 - value) for weight, value in zip(weights, frequencies))
    h_gamma = 2.0 * p_mean * (1.0 - p_mean)
    fst = None if h_gamma <= 0.0 else 1.0 - h_alpha / h_gamma
    polymorphic = epsilon < p_mean < 1.0 - epsilon and h_alpha > epsilon and h_gamma > epsilon
    alpha_preexisting = h_alpha <= spec.h_alpha_warning_threshold
    gamma_preexisting = h_gamma <= spec.h_gamma_warning_threshold
    return GeneticBaseline(
        scenario_id=scenario_id,
        high_allele_frequency_mean=p_mean,
        h_alpha=h_alpha,
        h_gamma=h_gamma,
        fst=fst,
        polymorphic=polymorphic,
        h_alpha_warning_preexisting=alpha_preexisting,
        h_gamma_warning_preexisting=gamma_preexisting,
        h2_dynamic_warning_eligible=polymorphic and not alpha_preexisting and not gamma_preexisting,
        h3_genetic_contrast_eligible=polymorphic,
    )


def _summarise(
    records: Sequence[PolymorphismEligibilityReplicate],
    master_seeds: Sequence[int],
) -> dict[str, object]:
    total = len(records)
    prepared = tuple(record for record in records if record.h1_full_state_prepared is True)
    by_scenario: dict[str, object] = {}
    for scenario_id in SCENARIO_IDS:
        values = tuple(record.baselines[scenario_id] for record in prepared if record.baselines is not None)
        by_scenario[scenario_id] = {
            "h1_full_state_prepared_count": len(values),
            "baseline_polymorphic_count": sum(value.polymorphic for value in values),
            "baseline_polymorphic_probability_across_all_seed_replicates": sum(value.polymorphic for value in values) / total,
            "baseline_polymorphic_probability_conditional_on_h1_full_state": None if not values else sum(value.polymorphic for value in values) / len(values),
            "h_alpha_warning_preexisting_count": sum(value.h_alpha_warning_preexisting for value in values),
            "h_gamma_warning_preexisting_count": sum(value.h_gamma_warning_preexisting for value in values),
            "h2_dynamic_warning_eligible_count": sum(value.h2_dynamic_warning_eligible for value in values),
            "h2_dynamic_warning_eligible_probability_across_all_seed_replicates": sum(value.h2_dynamic_warning_eligible for value in values) / total,
            "h3_genetic_contrast_eligible_count": sum(value.h3_genetic_contrast_eligible for value in values),
            "h3_genetic_contrast_eligible_probability_across_all_seed_replicates": sum(value.h3_genetic_contrast_eligible for value in values) / total,
            "high_allele_frequency_mean": _summary(value.high_allele_frequency_mean for value in values),
            "h_alpha": _summary(value.h_alpha for value in values),
            "h_gamma": _summary(value.h_gamma for value in values),
        }
    return {
        "denominators": {
            "total_seed_replicates": total,
            "h1_full_state_prepared_count": len(prepared),
            "h1_full_state_prepared_probability": len(prepared) / total,
        },
        "by_scenario": by_scenario,
        "by_master_seed": {
            str(seed): _seed_summary(tuple(record for record in records if record.master_seed == seed))
            for seed in master_seeds
        },
    }


def _seed_summary(records: Sequence[PolymorphismEligibilityReplicate]) -> dict[str, object]:
    if not records:
        return {"replicate_count": 0}
    prepared = tuple(record for record in records if record.h1_full_state_prepared is True)
    return {
        "replicate_count": len(records),
        "h1_full_state_prepared_probability": len(prepared) / len(records),
        "equal_isolated_h2_dynamic_warning_eligible_probability": None
        if not prepared
        else sum(record.baselines[SCENARIO_EQUAL_ISOLATED].h2_dynamic_warning_eligible for record in prepared if record.baselines is not None) / len(prepared),
        "equal_isolated_h3_genetic_contrast_eligible_probability": None
        if not prepared
        else sum(record.baselines[SCENARIO_EQUAL_ISOLATED].h3_genetic_contrast_eligible for record in prepared if record.baselines is not None) / len(prepared),
    }


def _summary(values: Iterable[float]) -> dict[str, float | None]:
    observed = tuple(float(value) for value in values)
    if not observed:
        return {"mean": None, "minimum": None, "maximum": None}
    return {"mean": sum(observed) / len(observed), "minimum": min(observed), "maximum": max(observed)}


def _flatten(mapping: Mapping[str, object], prefix: str = "") -> dict[str, object]:
    output: dict[str, object] = {}
    for key, value in mapping.items():
        name = f"{prefix}.{key}" if prefix else key
        if isinstance(value, Mapping):
            output.update(_flatten(value, name))
        else:
            output[name] = value
    return output
