"""Find mutation rates that preserve both finite H1 memory and genetic variation.

The zero-mutation finite closure supports robust H1 route memory but its high
branch is fixed for the high allele.  This audit introduces an explicit
symmetric allele-state mutation closure and asks a restricted preliminary
question for each predeclared mutation rate:

1. does the finite H1 high/low full-state hold still pass? and
2. is the held high state genetically eligible for later H2/H3 endpoints?

A rate is *screen-supported* only when both predicates hold in the same
seed-replicate.  No H2 timing or H3 fragmentation outcome is inferred here.
This is Type S evidence for the mutation closure, not a mutation-rate estimate
or a theorem about polymorphism.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import median
from typing import Iterable, Mapping, Sequence

from causal_model.finite_h1_continuation_state_audit import (
    DEFAULT_MASTER_SEEDS,
    ContinuationStateCell,
    ContinuationStateReplicate,
    run_finite_h1_continuation_state_audit,
)
from causal_model.multipatch_criticality_experiments import ExperimentSpec, ParameterCell
from causal_model.symmetric_allele_mutation_closure import patched_h1_mutation_runner, validate_symmetric_allele_mutation_rate

DEFAULT_MUTATION_RATES = (0.0, 0.05, 0.10, 0.15, 0.20)


@dataclass(frozen=True)
class MutationWindowReplicate:
    """One full-state H1 high-hold record under a declared mutation rate."""

    mutation_rate: float
    master_seed: int
    replicate_index: int
    calibration_seed: int
    h1_full_state_hold_supported: bool | None
    high_terminal_interaction_mean: float | None
    low_terminal_interaction_mean: float | None
    high_allele_frequency_mean: float | None
    h_alpha: float | None
    h_gamma: float | None
    polymorphic: bool | None
    h2_genetic_baseline_eligible: bool | None
    screen_supported: bool | None

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class MutationWindowCell:
    """Mutation-rate-specific H1 and polymorphism evidence for one parameter pair."""

    experiment_id: str
    profile: str
    mutation_rate: float
    pair_index: int
    parameters: ParameterCell
    master_seeds: tuple[int, ...]
    polymorphism_epsilon: float
    h_alpha_warning_threshold: float
    h_gamma_warning_threshold: float
    design: Mapping[str, object]
    replicates: tuple[MutationWindowReplicate, ...]
    summary: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "mutation_rate": self.mutation_rate,
            "pair_index": self.pair_index,
            "parameters": asdict(self.parameters),
            "master_seeds": list(self.master_seeds),
            "eligibility_rule": {
                "polymorphism_epsilon": self.polymorphism_epsilon,
                "h_alpha_warning_threshold": self.h_alpha_warning_threshold,
                "h_gamma_warning_threshold": self.h_gamma_warning_threshold,
                "h2_genetic_baseline_eligible": "epsilon < p < 1-epsilon and H-alpha/H-gamma are above their warning thresholds",
                "screen_supported": "finite H1 full-state hold and high-branch genetic baseline eligibility in the same record",
            },
            "design": dict(self.design),
            "replicates": [record.as_dict() for record in self.replicates],
            "summary": dict(self.summary),
        }

    def to_csv_row(self) -> dict[str, object]:
        row: dict[str, object] = {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "mutation_rate": self.mutation_rate,
            "pair_index": self.pair_index,
            "master_seeds": ",".join(str(seed) for seed in self.master_seeds),
            "replicate_count": len(self.replicates),
            "polymorphism_epsilon": self.polymorphism_epsilon,
            "h_alpha_warning_threshold": self.h_alpha_warning_threshold,
            "h_gamma_warning_threshold": self.h_gamma_warning_threshold,
        }
        row.update(asdict(self.parameters))
        row.update(_flatten(self.summary))
        return row


def run_finite_h1_mutation_window_audit(
    spec: ExperimentSpec,
    *,
    mutation_rates: Sequence[float] = DEFAULT_MUTATION_RATES,
    master_seeds: Sequence[int] = DEFAULT_MASTER_SEEDS,
    endpoint_padding_fraction: float = 0.5,
    stage_generations: int = 30,
    hold_generations: int = 30,
    nested_barrier_points: Sequence[int] = (25, 49, 97),
    interaction_separation_threshold: float = 0.05,
    maximum_normalized_bracket_width: float = 0.03,
    polymorphism_epsilon: float = 1e-12,
) -> tuple[MutationWindowCell, ...]:
    """Screen declared mutation rates for joint H1-memory/genetic eligibility."""
    rates = _validate_rates(mutation_rates)
    if not 0.0 < polymorphism_epsilon < 0.5:
        raise ValueError("polymorphism_epsilon must lie strictly between 0 and 0.5")
    output: list[MutationWindowCell] = []
    for rate in rates:
        with patched_h1_mutation_runner(rate):
            source_cells = run_finite_h1_continuation_state_audit(
                spec,
                master_seeds=master_seeds,
                endpoint_padding_fraction=endpoint_padding_fraction,
                stage_generations=stage_generations,
                hold_generations=hold_generations,
                nested_barrier_points=nested_barrier_points,
                interaction_separation_threshold=interaction_separation_threshold,
                maximum_normalized_bracket_width=maximum_normalized_bracket_width,
            )
        for source in source_cells:
            records = tuple(
                _classify_record(
                    rate,
                    record,
                    epsilon=polymorphism_epsilon,
                    h_alpha_threshold=spec.h_alpha_warning_threshold,
                    h_gamma_threshold=spec.h_gamma_warning_threshold,
                )
                for record in source.replicates
            )
            output.append(
                MutationWindowCell(
                    experiment_id=spec.experiment_id,
                    profile=spec.profile,
                    mutation_rate=rate,
                    pair_index=source.pair_index,
                    parameters=source.parameters,
                    master_seeds=source.master_seeds,
                    polymorphism_epsilon=polymorphism_epsilon,
                    h_alpha_warning_threshold=spec.h_alpha_warning_threshold,
                    h_gamma_warning_threshold=spec.h_gamma_warning_threshold,
                    design={
                        "mutation_timing": "after selection and migration, before finite drift",
                        "mutation_map": "p_mut = mu + (1 - 2 mu) p",
                        "endpoint_padding_fraction": source.endpoint_padding_fraction,
                        "stage_generations": source.stage_generations,
                        "hold_generations": source.hold_generations,
                        "nested_barrier_points": list(source.nested_barrier_points),
                        "interaction_separation_threshold": source.interaction_separation_threshold,
                    },
                    replicates=records,
                    summary=_summarise(records, source.master_seeds),
                )
            )
    return tuple(output)


def write_finite_h1_mutation_window_artifacts(
    cells: Iterable[MutationWindowCell],
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
    mutation_rate: float,
    record: ContinuationStateReplicate,
    *,
    epsilon: float,
    h_alpha_threshold: float,
    h_gamma_threshold: float,
) -> MutationWindowReplicate:
    high = record.rising_high_hold
    low = record.falling_low_hold
    if high is None or low is None:
        return MutationWindowReplicate(
            mutation_rate=mutation_rate,
            master_seed=record.master_seed,
            replicate_index=record.replicate_index,
            calibration_seed=record.calibration_seed,
            h1_full_state_hold_supported=record.full_state_hold_supported,
            high_terminal_interaction_mean=None,
            low_terminal_interaction_mean=None,
            high_allele_frequency_mean=None,
            h_alpha=None,
            h_gamma=None,
            polymorphic=None,
            h2_genetic_baseline_eligible=None,
            screen_supported=None,
        )
    p = high.terminal_high_allele_frequency_mean
    h = 2.0 * p * (1.0 - p)
    polymorphic = epsilon < p < 1.0 - epsilon and h > epsilon
    genetic_eligible = polymorphic and h > h_alpha_threshold and h > h_gamma_threshold
    h1_supported = record.full_state_hold_supported
    return MutationWindowReplicate(
        mutation_rate=mutation_rate,
        master_seed=record.master_seed,
        replicate_index=record.replicate_index,
        calibration_seed=record.calibration_seed,
        h1_full_state_hold_supported=h1_supported,
        high_terminal_interaction_mean=high.terminal_interaction_mean,
        low_terminal_interaction_mean=low.terminal_interaction_mean,
        high_allele_frequency_mean=p,
        h_alpha=h,
        h_gamma=h,
        polymorphic=polymorphic,
        h2_genetic_baseline_eligible=genetic_eligible,
        screen_supported=(h1_supported is True and genetic_eligible),
    )


def _validate_rates(values: Sequence[float]) -> tuple[float, ...]:
    rates = tuple(validate_symmetric_allele_mutation_rate(value) for value in values)
    if not rates or len(set(rates)) != len(rates):
        raise ValueError("mutation_rates must be nonempty and distinct")
    return rates


def _summarise(records: Sequence[MutationWindowReplicate], master_seeds: Sequence[int]) -> dict[str, object]:
    total = len(records)
    h1 = tuple(record for record in records if record.h1_full_state_hold_supported is True)
    polymorphic = tuple(record for record in records if record.polymorphic is True)
    eligible = tuple(record for record in records if record.h2_genetic_baseline_eligible is True)
    screen = tuple(record for record in records if record.screen_supported is True)
    return {
        "denominators": {
            "total_seed_replicates": total,
            "h1_full_state_hold_supported_count": len(h1),
            "h1_full_state_hold_supported_probability": len(h1) / total,
            "high_branch_polymorphic_count": len(polymorphic),
            "high_branch_polymorphic_probability": len(polymorphic) / total,
            "h2_genetic_baseline_eligible_count": len(eligible),
            "h2_genetic_baseline_eligible_probability": len(eligible) / total,
            "screen_supported_count": len(screen),
            "screen_supported_probability": len(screen) / total,
        },
        "high_branch_terminal_interaction": _summary(
            record.high_terminal_interaction_mean for record in records if record.high_terminal_interaction_mean is not None
        ),
        "low_branch_terminal_interaction": _summary(
            record.low_terminal_interaction_mean for record in records if record.low_terminal_interaction_mean is not None
        ),
        "high_branch_allele_frequency": _summary(
            record.high_allele_frequency_mean for record in records if record.high_allele_frequency_mean is not None
        ),
        "high_branch_H_alpha": _summary(record.h_alpha for record in records if record.h_alpha is not None),
        "by_master_seed": {
            str(seed): _seed_summary(tuple(record for record in records if record.master_seed == seed))
            for seed in master_seeds
        },
    }


def _seed_summary(records: Sequence[MutationWindowReplicate]) -> dict[str, object]:
    if not records:
        return {"replicate_count": 0}
    total = len(records)
    return {
        "replicate_count": total,
        "h1_full_state_hold_supported_probability": sum(record.h1_full_state_hold_supported is True for record in records) / total,
        "h2_genetic_baseline_eligible_probability": sum(record.h2_genetic_baseline_eligible is True for record in records) / total,
        "screen_supported_probability": sum(record.screen_supported is True for record in records) / total,
    }


def _summary(values: Iterable[float]) -> dict[str, float | None]:
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
