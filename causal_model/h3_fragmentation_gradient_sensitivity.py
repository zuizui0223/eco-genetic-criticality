"""Post-review H3 fragmentation-gradient sensitivity analysis.

This module extends the *frozen* mutation-primary H1/H3 closure with a new,
separately declared patch-count gradient.  It does not modify the canonical
H1/H3 evidence ledger.  Each successfully prepared H1 high full state is
projected to equal isolated landscapes with several patch counts at fixed total
area, preserving the same source across the repeated patch-count outcomes.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Iterable, Sequence

from causal_model.finite_h1_boundary_resolution_audit import run_finite_h1_boundary_resolution_audit
from causal_model.finite_h1_fragment_projection_audit import _mean, project_full_state
from causal_model.multipatch_criticality_experiments import (
    ExperimentSpec,
    LandscapeScenario,
    parameters_for_cell,
    summarise_replicate,
)
from causal_model.mutation_h1_primary_domain import primary_analysis_cells
from causal_model.mutation_primary_h1_h2_h3_chain import _prepare_mutation_high_state
from causal_model.symmetric_allele_mutation_closure import (
    patched_h1_mutation_runner,
    simulate_with_symmetric_allele_mutation,
)

FRESH_GRADIENT_MASTER_SEEDS = (20260820, 20260821, 20260822, 20260823, 20260824)
FRAGMENT_PATCH_COUNTS = (1, 2, 3, 4, 6, 8, 12, 16)
CANONICAL_SCIENTIFIC_COMMIT = "dd8ee379d0d3518194c767d16402042525bc00dc"


@dataclass(frozen=True)
class FragmentationGradientRecord:
    primary_cell_index: int
    mutation_rate: float
    area_reference: float
    interaction_feedback: float
    master_seed: int
    replicate_index: int
    calibration_seed: int
    h1_resolution_supported: bool | None
    source_prepared: bool
    anchor_barrier: float | None
    patch_count: int
    patch_area: float
    theoretical_k_fragment: float
    theoretical_k_above_four: bool
    outcome_seed: int | None
    projection_supported: bool | None
    final_interaction_mean: float | None
    final_effective_size_mean: float | None
    realised_high_trait_mass_mean: float | None
    potential_high_trait_viable: bool | None
    realised_high_trait_persists: bool | None
    h_alpha: float | None
    h_gamma: float | None

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class FragmentationGradientArtifact:
    primary_cell_index: int
    protocol: dict[str, object]
    records: tuple[FragmentationGradientRecord, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "primary_cell_index": self.primary_cell_index,
            "protocol": dict(self.protocol),
            "record_count": len(self.records),
            "records": [record.as_dict() for record in self.records],
        }


def equal_isolated_scenario(total_area: float, patch_count: int) -> LandscapeScenario:
    """Build an equal isolated landscape at fixed total area."""
    if total_area <= 0.0:
        raise ValueError("total_area must be positive")
    if patch_count < 1:
        raise ValueError("patch_count must be positive")
    patch_area = total_area / patch_count
    return LandscapeScenario(
        scenario_id=f"equal_isolated_n{patch_count}",
        patch_areas=tuple(patch_area for _ in range(patch_count)),
        migration_rate=0.0,
    )


def theoretical_fragment_k(
    *, interaction_feedback: float, patch_area: float, area_reference: float
) -> float:
    if area_reference <= 0.0:
        raise ValueError("area_reference must be positive")
    return interaction_feedback * patch_area / area_reference


def run_h3_fragmentation_gradient_sensitivity(
    spec: ExperimentSpec,
    *,
    primary_cell_indices: Sequence[int] | None = None,
    master_seeds: Sequence[int] = FRESH_GRADIENT_MASTER_SEEDS,
    patch_counts: Sequence[int] = FRAGMENT_PATCH_COUNTS,
    endpoint_padding_fraction: float = 0.5,
    stage_generations: int = 30,
    hold_generations: int = 30,
    nested_barrier_points: Sequence[int] = (25, 49, 97),
    interaction_separation_threshold: float = 0.05,
    maximum_normalized_bracket_width: float = 0.03,
) -> tuple[FragmentationGradientArtifact, ...]:
    """Run the preregistered fixed-total-area patch-count gradient.

    Source preparation is replayed independently for the fresh seed family.
    Patch-count outcomes from one prepared source are repeated measures.
    """
    seeds = _validate_unique_nonnegative(master_seeds, "master_seeds")
    counts = _validate_patch_counts(patch_counts)
    domains = tuple(primary_analysis_cells())
    if primary_cell_indices is None:
        indices = tuple(range(len(domains)))
    else:
        indices = tuple(int(index) for index in primary_cell_indices)
        if not indices or len(indices) != len(set(indices)):
            raise ValueError("primary_cell_indices must be nonempty and unique")
        if any(index < 0 or index >= len(domains) for index in indices):
            raise ValueError("primary_cell_indices contains an out-of-range index")

    artifacts: list[FragmentationGradientArtifact] = []
    for cell_index in indices:
        domain = domains[cell_index]
        records: list[FragmentationGradientRecord] = []
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
                    nested_barrier_points=tuple(nested_barrier_points),
                    interaction_separation_threshold=interaction_separation_threshold,
                    maximum_normalized_bracket_width=maximum_normalized_bracket_width,
                )
            if len(calibration) != 1:
                raise RuntimeError("targeted gradient calibration must yield exactly one parameter cell")
            source_cell = calibration[0]
            for calibration_record in source_cell.replicates:
                prepared = _prepare_mutation_high_state(
                    domain.mutation_rate,
                    target_spec,
                    source_cell,
                    calibration_record,
                    endpoint_padding_fraction=endpoint_padding_fraction,
                    stage_generations=stage_generations,
                    hold_generations=hold_generations,
                    interaction_separation_threshold=interaction_separation_threshold,
                )
                h1_supported = calibration_record.resolution_stable_h1_loop_mechanism_supported
                if prepared is None:
                    for patch_count in counts:
                        patch_area = spec.total_area / patch_count
                        k_fragment = theoretical_fragment_k(
                            interaction_feedback=domain.interaction_feedback,
                            patch_area=patch_area,
                            area_reference=domain.area_reference,
                        )
                        records.append(
                            FragmentationGradientRecord(
                                primary_cell_index=cell_index,
                                mutation_rate=domain.mutation_rate,
                                area_reference=domain.area_reference,
                                interaction_feedback=domain.interaction_feedback,
                                master_seed=master_seed,
                                replicate_index=calibration_record.replicate_index,
                                calibration_seed=calibration_record.seed,
                                h1_resolution_supported=h1_supported,
                                source_prepared=False,
                                anchor_barrier=None,
                                patch_count=patch_count,
                                patch_area=patch_area,
                                theoretical_k_fragment=k_fragment,
                                theoretical_k_above_four=k_fragment > 4.0,
                                outcome_seed=None,
                                projection_supported=None,
                                final_interaction_mean=None,
                                final_effective_size_mean=None,
                                realised_high_trait_mass_mean=None,
                                potential_high_trait_viable=None,
                                realised_high_trait_persists=None,
                                h_alpha=None,
                                h_gamma=None,
                            )
                        )
                    continue

                source, anchor = prepared
                anchor_cell = replace(source_cell.parameters, interaction_barrier=anchor)
                for patch_count in counts:
                    scenario = equal_isolated_scenario(spec.total_area, patch_count)
                    outcome_seed = _outcome_seed(calibration_record.seed, patch_count)
                    template = parameters_for_cell(
                        target_spec,
                        scenario,
                        anchor_cell,
                        seed=outcome_seed,
                    )
                    projected, invariant = project_full_state(source, template)
                    if not invariant.projection_supported:
                        summary = None
                    else:
                        result = simulate_with_symmetric_allele_mutation(
                            replace(
                                projected,
                                generations=spec.generations,
                                random_seed=outcome_seed,
                            ),
                            mutation_rate=domain.mutation_rate,
                        )
                        summary = summarise_replicate(
                            result,
                            replicate_index=calibration_record.replicate_index,
                            seed=outcome_seed,
                            h_alpha_warning_threshold=spec.h_alpha_warning_threshold,
                            h_gamma_warning_threshold=spec.h_gamma_warning_threshold,
                            fst_warning_threshold=spec.fst_warning_threshold,
                            allele_loss_threshold=spec.allele_loss_threshold,
                        )
                    patch_area = spec.total_area / patch_count
                    k_fragment = theoretical_fragment_k(
                        interaction_feedback=domain.interaction_feedback,
                        patch_area=patch_area,
                        area_reference=domain.area_reference,
                    )
                    records.append(
                        FragmentationGradientRecord(
                            primary_cell_index=cell_index,
                            mutation_rate=domain.mutation_rate,
                            area_reference=domain.area_reference,
                            interaction_feedback=domain.interaction_feedback,
                            master_seed=master_seed,
                            replicate_index=calibration_record.replicate_index,
                            calibration_seed=calibration_record.seed,
                            h1_resolution_supported=h1_supported,
                            source_prepared=True,
                            anchor_barrier=anchor,
                            patch_count=patch_count,
                            patch_area=patch_area,
                            theoretical_k_fragment=k_fragment,
                            theoretical_k_above_four=k_fragment > 4.0,
                            outcome_seed=outcome_seed,
                            projection_supported=invariant.projection_supported,
                            final_interaction_mean=None if summary is None else _mean(summary.final_q_by_patch),
                            final_effective_size_mean=None if summary is None else _mean(summary.final_effective_size_by_patch),
                            realised_high_trait_mass_mean=None if summary is None else summary.realised_high_trait_mass_mean,
                            potential_high_trait_viable=None if summary is None else summary.potential_high_trait_viable,
                            realised_high_trait_persists=None if summary is None else summary.realised_high_trait_persists,
                            h_alpha=None if summary is None else summary.h_alpha,
                            h_gamma=None if summary is None else summary.h_gamma,
                        )
                    )
        protocol = {
            "analysis": "post-review H3 fragmentation-gradient sensitivity",
            "canonical_scientific_commit": CANONICAL_SCIENTIFIC_COMMIT,
            "fresh_master_seeds": list(seeds),
            "patch_counts": list(counts),
            "total_area": spec.total_area,
            "outcome_generations": spec.generations,
            "source_replicates_per_master_seed": spec.replicates,
            "endpoint_padding_fraction": endpoint_padding_fraction,
            "stage_generations": stage_generations,
            "hold_generations": hold_generations,
            "nested_barrier_points": list(nested_barrier_points),
            "interaction_separation_threshold": interaction_separation_threshold,
            "maximum_normalized_bracket_width": maximum_normalized_bracket_width,
            "warning_endpoints_evaluated": False,
            "evidence_label": "new supplementary finite Type S sensitivity evidence",
        }
        artifacts.append(FragmentationGradientArtifact(cell_index, protocol, tuple(records)))
    return tuple(artifacts)


def write_fragmentation_gradient_artifacts(
    artifacts: Iterable[FragmentationGradientArtifact],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> None:
    values = tuple(artifacts)
    if not values:
        raise ValueError("artifacts must be nonempty")
    rows = [record.as_dict() for artifact in values for record in artifact.records]
    csv_target = Path(csv_path)
    json_target = Path(json_path)
    csv_target.parent.mkdir(parents=True, exist_ok=True)
    json_target.parent.mkdir(parents=True, exist_ok=True)
    with csv_target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    json_target.write_text(
        json.dumps([artifact.as_dict() for artifact in values], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _outcome_seed(calibration_seed: int, patch_count: int) -> int:
    return (int(calibration_seed) * 1_000_003 + 170_003 + int(patch_count) * 101) % (2**31 - 1)


def _validate_unique_nonnegative(values: Sequence[int], label: str) -> tuple[int, ...]:
    result = tuple(int(value) for value in values)
    if not result or len(result) != len(set(result)) or any(value < 0 for value in result):
        raise ValueError(f"{label} must be nonempty, unique, nonnegative integers")
    return result


def _validate_patch_counts(values: Sequence[int]) -> tuple[int, ...]:
    result = tuple(int(value) for value in values)
    if not result or len(result) != len(set(result)) or any(value < 1 for value in result):
        raise ValueError("patch_counts must be nonempty, unique, positive integers")
    if tuple(sorted(result)) != result:
        raise ValueError("patch_counts must be in ascending order")
    if result[0] != 1:
        raise ValueError("patch_counts must begin with the paired one-patch reference")
    return result
