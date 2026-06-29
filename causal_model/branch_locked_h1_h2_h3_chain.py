"""Branch-locked, same-replicate finite H1--H2--H3 chain audit.

The H1 validation series established a finite collapse--recovery loop that is
robust to stage duration, endpoint range, nested barrier-grid refinement, and a
predeclared master-seed ensemble.  That is evidence about H1 itself; it does
*not* establish that H2 warning order and H3 fragmentation consequences occur in
the same finite replicate.

This module runs those statements in a fixed order:

1. run the declared 97-point H1 boundary-resolution calibration for each master
   seed and replicate;
2. take the conservative interior of that replicate's finite H1 bracket as a
   branch-anchor barrier;
3. initialize canonical low/high H1 starts at that anchor, using the exact same
   derived seed as the H1 calibration;
4. evaluate H2 only in branch-locked replicates and only conditionally on
   observed trait loss, with censored first-passage pairs retained explicitly;
5. compare the high-start one-large, equal-isolated, and equal-migrating
   landscapes in the same replicate; and
6. report strictly same-replicate H1 -> H2 -> H3 predicates separately for
   H-alpha, H-gamma, and their predeclared union.

Migration in this finite multipatch closure is allele-frequency mixing.  The
module reports FST and allele-frequency modulation relative to isolation, but
never labels migration an ecological or demographic rescue.

All output is Type S evidence for this declared numerical closure.  Neither a
closed finite loop nor complete same-replicate support proves a general
bifurcation theorem or a universal ecological sequence.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from statistics import median
from typing import Iterable, Mapping, Sequence

from causal_model.canonical_h1_bifurcation import CanonicalH1Certificate, canonical_h1_certificate
from causal_model.finite_h1_boundary_resolution_audit import (
    DEFAULT_NESTED_BARRIER_POINTS,
    FiniteH1BoundaryResolutionCell,
    FiniteH1BoundaryResolutionReplicate,
    run_finite_h1_boundary_resolution_audit,
)
from causal_model.multipatch_criticality_dynamics import FirstPassageEvent, simulate
from causal_model.multipatch_criticality_experiments import (
    SCENARIO_EQUAL_ISOLATED,
    SCENARIO_EQUAL_MIGRATING,
    SCENARIO_ONE_LARGE,
    ExperimentSpec,
    LandscapeScenario,
    ParameterCell,
    ReplicateSummary,
    default_scenarios,
    parameters_for_cell,
    summarise_replicate,
)

DEFAULT_CHAIN_MASTER_SEEDS = (20260630, 20260631, 20260632, 20260633, 20260634)
_BRANCHES = ("low_start", "high_start")
_SCENARIOS = (SCENARIO_ONE_LARGE, SCENARIO_EQUAL_ISOLATED, SCENARIO_EQUAL_MIGRATING)
_H2_WARNING_ATTRIBUTES = (("H_alpha", "tau_H_alpha"), ("H_gamma", "tau_H_gamma"))


@dataclass(frozen=True)
class H2FirstPassageComparison:
    """One warning-versus-realised-trait comparison with explicit censoring."""

    warning_id: str
    warning_time: int | None
    trait_loss_time: int | None
    warning_threshold: float | int | None
    warning_aggregation_rule: str | None
    trait_threshold: float | int | None
    trait_aggregation_rule: str | None
    trait_loss_observed: bool
    warning_observed: bool
    valid_pair: bool
    censored: bool
    warning_leads: bool | None
    lead_time_trait_minus_warning: int | None

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class BranchLockedOutcome:
    """One landscape outcome under an H1 anchor and declared low/high start."""

    scenario_id: str
    branch_id: str
    summary: ReplicateSummary
    terminal_interaction_mean: float
    terminal_local_effective_size_mean: float
    terminal_allele_frequency_mean: float
    terminal_potential_high_trait_viable: bool
    h2_comparisons: tuple[H2FirstPassageComparison, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "scenario_id": self.scenario_id,
            "branch_id": self.branch_id,
            "summary": self.summary.as_dict(),
            "terminal_interaction_mean": self.terminal_interaction_mean,
            "terminal_local_effective_size_mean": self.terminal_local_effective_size_mean,
            "terminal_allele_frequency_mean": self.terminal_allele_frequency_mean,
            "terminal_potential_high_trait_viable": self.terminal_potential_high_trait_viable,
            "h2_comparisons": [comparison.as_dict() for comparison in self.h2_comparisons],
        }


@dataclass(frozen=True)
class HighBranchH3Contrast:
    """High-start fragmentation contrast relative to the matched one-large run."""

    interaction_difference_isolated_minus_large: float
    local_effective_size_difference_isolated_minus_large: float
    realised_high_trait_mass_difference_isolated_minus_large: float
    h_alpha_difference_isolated_minus_large: float
    h_gamma_difference_isolated_minus_large: float
    fragmentation_pattern_supported: bool

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class MigrationMixingContrast:
    """High-start migration modulation, explicitly not a rescue predicate."""

    interaction_difference_migrating_minus_isolated: float
    local_effective_size_difference_migrating_minus_isolated: float
    realised_high_trait_mass_difference_migrating_minus_isolated: float
    allele_frequency_mean_difference_migrating_minus_isolated: float
    fst_difference_migrating_minus_isolated: float | None
    fst_lower_with_migration: bool | None
    h_alpha_difference_migrating_minus_isolated: float
    h_gamma_difference_migrating_minus_isolated: float

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class BranchLockedChainReplicate:
    """One seed- and replicate-paired H1 calibration, H2, H3, and chain record."""

    master_seed: int
    replicate_index: int
    seed: int
    h1_resolution_stable_mechanism_supported: bool | None
    h1_anchor_lower_barrier: float | None
    h1_anchor_upper_barrier: float | None
    h1_anchor_barrier: float | None
    h1_anchor_bracket_width: float | None
    canonical_anchor_h1: CanonicalH1Certificate | None
    low_initial_interaction: float | None
    high_initial_interaction: float | None
    outcomes: Mapping[str, Mapping[str, BranchLockedOutcome]] | None
    h1_branch_lock_supported: bool | None
    high_branch_h3: HighBranchH3Contrast | None
    migration_mixing: MigrationMixingContrast | None
    isolated_high_h_alpha_leads_trait: bool | None
    isolated_high_h_gamma_leads_trait: bool | None
    same_replicate_chain_h_alpha_supported: bool | None
    same_replicate_chain_h_gamma_supported: bool | None
    same_replicate_chain_any_warning_supported: bool | None

    def as_dict(self) -> dict[str, object]:
        return {
            "master_seed": self.master_seed,
            "replicate_index": self.replicate_index,
            "seed": self.seed,
            "h1_resolution_stable_mechanism_supported": self.h1_resolution_stable_mechanism_supported,
            "h1_anchor": {
                "lower_barrier": self.h1_anchor_lower_barrier,
                "upper_barrier": self.h1_anchor_upper_barrier,
                "anchor_barrier": self.h1_anchor_barrier,
                "conservative_interior_width": self.h1_anchor_bracket_width,
            },
            "canonical_anchor_h1": None if self.canonical_anchor_h1 is None else asdict(self.canonical_anchor_h1),
            "low_initial_interaction": self.low_initial_interaction,
            "high_initial_interaction": self.high_initial_interaction,
            "outcomes": None
            if self.outcomes is None
            else {
                scenario_id: {branch_id: outcome.as_dict() for branch_id, outcome in branch_map.items()}
                for scenario_id, branch_map in self.outcomes.items()
            },
            "h1_branch_lock_supported": self.h1_branch_lock_supported,
            "high_branch_h3": None if self.high_branch_h3 is None else self.high_branch_h3.as_dict(),
            "migration_mixing": None if self.migration_mixing is None else self.migration_mixing.as_dict(),
            "isolated_high_h_alpha_leads_trait": self.isolated_high_h_alpha_leads_trait,
            "isolated_high_h_gamma_leads_trait": self.isolated_high_h_gamma_leads_trait,
            "same_replicate_chain_h_alpha_supported": self.same_replicate_chain_h_alpha_supported,
            "same_replicate_chain_h_gamma_supported": self.same_replicate_chain_h_gamma_supported,
            "same_replicate_chain_any_warning_supported": self.same_replicate_chain_any_warning_supported,
        }


@dataclass(frozen=True)
class BranchLockedChainCell:
    """Multi-master-seed same-replicate H1--H2--H3 evidence for one parameter pair."""

    experiment_id: str
    profile: str
    pair_index: int
    parameters: ParameterCell
    canonical_bistable_barrier_interval: tuple[float, float] | None
    master_seeds: tuple[int, ...]
    h1_endpoint_padding_fraction: float
    h1_stage_generations: int
    h1_nested_barrier_points: tuple[int, ...]
    interaction_separation_threshold: float
    replicates: tuple[BranchLockedChainReplicate, ...]
    summary: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "pair_index": self.pair_index,
            "parameters": asdict(self.parameters),
            "canonical_bistable_barrier_interval": None
            if self.canonical_bistable_barrier_interval is None
            else list(self.canonical_bistable_barrier_interval),
            "master_seeds": list(self.master_seeds),
            "h1_calibration": {
                "endpoint_padding_fraction": self.h1_endpoint_padding_fraction,
                "stage_generations": self.h1_stage_generations,
                "nested_barrier_points": list(self.h1_nested_barrier_points),
                "interaction_separation_threshold": self.interaction_separation_threshold,
            },
            "replicate_count": len(self.replicates),
            "replicates": [replicate.as_dict() for replicate in self.replicates],
            "summary": dict(self.summary),
        }

    def to_csv_row(self) -> dict[str, object]:
        row: dict[str, object] = {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "pair_index": self.pair_index,
            "master_seed_count": len(self.master_seeds),
            "master_seeds": ",".join(str(seed) for seed in self.master_seeds),
            "replicate_count": len(self.replicates),
            "h1_endpoint_padding_fraction": self.h1_endpoint_padding_fraction,
            "h1_stage_generations": self.h1_stage_generations,
            "h1_nested_barrier_points": ",".join(str(value) for value in self.h1_nested_barrier_points),
            "interaction_separation_threshold": self.interaction_separation_threshold,
            "canonical_bistable_barrier_lower": None
            if self.canonical_bistable_barrier_interval is None
            else self.canonical_bistable_barrier_interval[0],
            "canonical_bistable_barrier_upper": None
            if self.canonical_bistable_barrier_interval is None
            else self.canonical_bistable_barrier_interval[1],
        }
        row.update(asdict(self.parameters))
        row.update(_flatten_mapping(self.summary))
        return row


def run_branch_locked_h1_h2_h3_chain_audit(
    spec: ExperimentSpec,
    *,
    master_seeds: Sequence[int] = DEFAULT_CHAIN_MASTER_SEEDS,
    h1_endpoint_padding_fraction: float = 0.5,
    h1_stage_generations: int = 30,
    h1_nested_barrier_points: Sequence[int] = DEFAULT_NESTED_BARRIER_POINTS,
    interaction_separation_threshold: float = 0.05,
    h1_maximum_normalized_bracket_width: float = 0.03,
) -> tuple[BranchLockedChainCell, ...]:
    """Run H1 calibration, then same-replicate H2 and high-branch H3 tests.

    The H1 calibration is an explicit precondition rather than an outcome filter.
    A calibration replicate without a resolution-stable finite mechanism remains
    in the returned record with H2/H3 values unavailable.  H2 first-passage
    comparisons are valid only when both the warning and realised trait-loss
    event occur; censoring is retained, never converted to an endpoint time.
    """
    seeds = _validate_master_seeds(master_seeds)
    if h1_endpoint_padding_fraction <= 0.0:
        raise ValueError("h1_endpoint_padding_fraction must be positive")
    if h1_stage_generations < 1:
        raise ValueError("h1_stage_generations must be positive")
    if interaction_separation_threshold < 0.0:
        raise ValueError("interaction_separation_threshold must be non-negative")
    scenarios = _scenario_map(spec)
    h1_campaigns = tuple(
        run_finite_h1_boundary_resolution_audit(
            replace(spec, master_seed=master_seed),
            endpoint_padding_fraction=h1_endpoint_padding_fraction,
            stage_generations=h1_stage_generations,
            nested_barrier_points=h1_nested_barrier_points,
            interaction_separation_threshold=interaction_separation_threshold,
            maximum_normalized_bracket_width=h1_maximum_normalized_bracket_width,
        )
        for master_seed in seeds
    )
    count = len(h1_campaigns[0])
    if any(len(campaign) != count for campaign in h1_campaigns):
        raise RuntimeError("H1 calibration campaigns produced unequal parameter-cell counts")

    output: list[BranchLockedChainCell] = []
    for index in range(count):
        calibration_cells = tuple(campaign[index] for campaign in h1_campaigns)
        reference = calibration_cells[0]
        if any(_calibration_identity(cell) != _calibration_identity(reference) for cell in calibration_cells[1:]):
            raise RuntimeError("H1 calibration cells are not aligned across master seeds")
        replicates = tuple(
            _run_calibrated_replicate(master_seed, seed_spec=replace(spec, master_seed=master_seed), calibration=cell, scenarios=scenarios,
                interaction_separation_threshold=interaction_separation_threshold)
            for master_seed, cell in zip(seeds, calibration_cells, strict=True)
            for _unused in (0,)
            for _record in ()
        )
        # The expression above is intentionally replaced immediately below: one
        # cell needs all of its replicate records, not a nested seed-cell shell.
        records: list[BranchLockedChainReplicate] = []
        for master_seed, calibration in zip(seeds, calibration_cells, strict=True):
            seed_spec = replace(spec, master_seed=master_seed)
            records.extend(
                _run_calibrated_replicate(
                    master_seed,
                    seed_spec=seed_spec,
                    calibration=calibration,
                    scenarios=scenarios,
                    interaction_separation_threshold=interaction_separation_threshold,
                    calibration_replicate=calibration_replicate,
                )
                for calibration_replicate in calibration.replicates
            )
        del replicates
        output.append(
            BranchLockedChainCell(
                experiment_id=spec.experiment_id,
                profile=spec.profile,
                pair_index=reference.pair_index,
                parameters=reference.parameters,
                canonical_bistable_barrier_interval=reference.canonical_bistable_barrier_interval,
                master_seeds=seeds,
                h1_endpoint_padding_fraction=h1_endpoint_padding_fraction,
                h1_stage_generations=h1_stage_generations,
                h1_nested_barrier_points=tuple(int(value) for value in h1_nested_barrier_points),
                interaction_separation_threshold=interaction_separation_threshold,
                replicates=tuple(records),
                summary=_summarise_cell(records, seeds),
            )
        )
    return tuple(output)


def write_branch_locked_h1_h2_h3_chain_artifacts(
    cells: Iterable[BranchLockedChainCell],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> None:
    """Write an aligned ledger and complete same-replicate raw evidence."""
    values = tuple(cells)
    if not values:
        raise ValueError("cells must be nonempty")
    csv_target = Path(csv_path)
    json_target = Path(json_path)
    csv_target.parent.mkdir(parents=True, exist_ok=True)
    json_target.parent.mkdir(parents=True, exist_ok=True)
    rows = [cell.to_csv_row() for cell in values]
    fields = sorted({field for row in rows for field in row})
    with csv_target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with json_target.open("w", encoding="utf-8") as handle:
        json.dump([cell.as_dict() for cell in values], handle, indent=2, sort_keys=True)


def _run_calibrated_replicate(
    master_seed: int,
    *,
    seed_spec: ExperimentSpec,
    calibration: FiniteH1BoundaryResolutionCell,
    scenarios: Mapping[str, LandscapeScenario],
    interaction_separation_threshold: float,
    calibration_replicate: FiniteH1BoundaryResolutionReplicate,
) -> BranchLockedChainReplicate:
    observation = calibration_replicate.observations[-1]
    if calibration_replicate.resolution_stable_h1_loop_mechanism_supported is not True:
        return _unavailable_record(master_seed, calibration_replicate)
    if not observation.finite_h1_loop_bracket_mechanism_supported:
        return _unavailable_record(master_seed, calibration_replicate)
    collapse = observation.rising_collapse_bracket
    recovery = observation.falling_recovery_bracket
    if collapse is None or recovery is None:
        raise RuntimeError("resolution-stable H1 calibration lacks finite boundary brackets")
    anchor_lower = recovery.upper_barrier
    anchor_upper = collapse.lower_barrier
    if not anchor_lower < anchor_upper:
        return _unavailable_record(
            master_seed,
            calibration_replicate,
            anchor_lower=anchor_lower,
            anchor_upper=anchor_upper,
        )
    anchor = (anchor_lower + anchor_upper) / 2.0
    anchor_cell = replace(calibration.parameters, interaction_barrier=anchor)
    one_large = scenarios[SCENARIO_ONE_LARGE]
    base = parameters_for_cell(seed_spec, one_large, anchor_cell, seed=calibration_replicate.seed)
    canonical = canonical_h1_certificate(
        feedback_strength=anchor_cell.interaction_feedback,
        area=one_large.patch_areas[0],
        area_reference=anchor_cell.area_reference,
        barrier=anchor,
        trait_parameters=base,
    )
    if not canonical.branch_dependent_high_trait_mode:
        return _unavailable_record(
            master_seed,
            calibration_replicate,
            anchor_lower=anchor_lower,
            anchor_upper=anchor_upper,
            anchor=anchor,
            canonical=canonical,
        )
    assert canonical.low_stable_branch is not None
    assert canonical.high_stable_branch is not None
    low_initial = canonical.low_stable_branch.interaction
    high_initial = canonical.high_stable_branch.interaction
    outcomes = {
        scenario_id: {
            "low_start": _run_outcome(seed_spec, scenario, anchor_cell, calibration_replicate.seed, calibration_replicate.replicate_index, "low_start", low_initial),
            "high_start": _run_outcome(seed_spec, scenario, anchor_cell, calibration_replicate.seed, calibration_replicate.replicate_index, "high_start", high_initial),
        }
        for scenario_id, scenario in scenarios.items()
    }
    large_low = outcomes[SCENARIO_ONE_LARGE]["low_start"]
    large_high = outcomes[SCENARIO_ONE_LARGE]["high_start"]
    branch_lock = (
        large_high.terminal_interaction_mean - large_low.terminal_interaction_mean > interaction_separation_threshold
        and large_high.terminal_potential_high_trait_viable
        and not large_low.terminal_potential_high_trait_viable
    )
    if not branch_lock:
        return BranchLockedChainReplicate(
            master_seed=master_seed,
            replicate_index=calibration_replicate.replicate_index,
            seed=calibration_replicate.seed,
            h1_resolution_stable_mechanism_supported=True,
            h1_anchor_lower_barrier=anchor_lower,
            h1_anchor_upper_barrier=anchor_upper,
            h1_anchor_barrier=anchor,
            h1_anchor_bracket_width=anchor_upper - anchor_lower,
            canonical_anchor_h1=canonical,
            low_initial_interaction=low_initial,
            high_initial_interaction=high_initial,
            outcomes=outcomes,
            h1_branch_lock_supported=False,
            high_branch_h3=None,
            migration_mixing=None,
            isolated_high_h_alpha_leads_trait=None,
            isolated_high_h_gamma_leads_trait=None,
            same_replicate_chain_h_alpha_supported=False,
            same_replicate_chain_h_gamma_supported=False,
            same_replicate_chain_any_warning_supported=False,
        )
    large_high = outcomes[SCENARIO_ONE_LARGE]["high_start"]
    isolated_high = outcomes[SCENARIO_EQUAL_ISOLATED]["high_start"]
    migrating_high = outcomes[SCENARIO_EQUAL_MIGRATING]["high_start"]
    h3 = HighBranchH3Contrast(
        interaction_difference_isolated_minus_large=isolated_high.terminal_interaction_mean - large_high.terminal_interaction_mean,
        local_effective_size_difference_isolated_minus_large=(
            isolated_high.terminal_local_effective_size_mean - large_high.terminal_local_effective_size_mean
        ),
        realised_high_trait_mass_difference_isolated_minus_large=(
            isolated_high.summary.realised_high_trait_mass_mean - large_high.summary.realised_high_trait_mass_mean
        ),
        h_alpha_difference_isolated_minus_large=isolated_high.summary.h_alpha - large_high.summary.h_alpha,
        h_gamma_difference_isolated_minus_large=isolated_high.summary.h_gamma - large_high.summary.h_gamma,
        fragmentation_pattern_supported=(
            isolated_high.terminal_interaction_mean < large_high.terminal_interaction_mean
            and isolated_high.terminal_local_effective_size_mean < large_high.terminal_local_effective_size_mean
            and isolated_high.summary.realised_high_trait_mass_mean < large_high.summary.realised_high_trait_mass_mean
        ),
    )
    migration = _migration_contrast(migrating_high, isolated_high)
    by_warning = {comparison.warning_id: comparison for comparison in isolated_high.h2_comparisons}
    alpha_lead = by_warning["H_alpha"].warning_leads
    gamma_lead = by_warning["H_gamma"].warning_leads
    return BranchLockedChainReplicate(
        master_seed=master_seed,
        replicate_index=calibration_replicate.replicate_index,
        seed=calibration_replicate.seed,
        h1_resolution_stable_mechanism_supported=True,
        h1_anchor_lower_barrier=anchor_lower,
        h1_anchor_upper_barrier=anchor_upper,
        h1_anchor_barrier=anchor,
        h1_anchor_bracket_width=anchor_upper - anchor_lower,
        canonical_anchor_h1=canonical,
        low_initial_interaction=low_initial,
        high_initial_interaction=high_initial,
        outcomes=outcomes,
        h1_branch_lock_supported=True,
        high_branch_h3=h3,
        migration_mixing=migration,
        isolated_high_h_alpha_leads_trait=alpha_lead,
        isolated_high_h_gamma_leads_trait=gamma_lead,
        same_replicate_chain_h_alpha_supported=h3.fragmentation_pattern_supported and alpha_lead is True,
        same_replicate_chain_h_gamma_supported=h3.fragmentation_pattern_supported and gamma_lead is True,
        same_replicate_chain_any_warning_supported=(
            h3.fragmentation_pattern_supported and (alpha_lead is True or gamma_lead is True)
        ),
    )


def _unavailable_record(
    master_seed: int,
    calibration_replicate: FiniteH1BoundaryResolutionReplicate,
    *,
    anchor_lower: float | None = None,
    anchor_upper: float | None = None,
    anchor: float | None = None,
    canonical: CanonicalH1Certificate | None = None,
) -> BranchLockedChainReplicate:
    return BranchLockedChainReplicate(
        master_seed=master_seed,
        replicate_index=calibration_replicate.replicate_index,
        seed=calibration_replicate.seed,
        h1_resolution_stable_mechanism_supported=calibration_replicate.resolution_stable_h1_loop_mechanism_supported,
        h1_anchor_lower_barrier=anchor_lower,
        h1_anchor_upper_barrier=anchor_upper,
        h1_anchor_barrier=anchor,
        h1_anchor_bracket_width=None if anchor_lower is None or anchor_upper is None else anchor_upper - anchor_lower,
        canonical_anchor_h1=canonical,
        low_initial_interaction=None,
        high_initial_interaction=None,
        outcomes=None,
        h1_branch_lock_supported=None,
        high_branch_h3=None,
        migration_mixing=None,
        isolated_high_h_alpha_leads_trait=None,
        isolated_high_h_gamma_leads_trait=None,
        same_replicate_chain_h_alpha_supported=None,
        same_replicate_chain_h_gamma_supported=None,
        same_replicate_chain_any_warning_supported=None,
    )


def _run_outcome(
    spec: ExperimentSpec,
    scenario: LandscapeScenario,
    cell: ParameterCell,
    seed: int,
    replicate_index: int,
    branch_id: str,
    initial_interaction: float,
) -> BranchLockedOutcome:
    parameters = parameters_for_cell(spec, scenario, cell, seed=seed)
    parameters = replace(parameters, initial_interaction=tuple(initial_interaction for _ in parameters.patch_areas))
    result = simulate(parameters)
    summary = summarise_replicate(
        result,
        replicate_index=replicate_index,
        seed=seed,
        h_alpha_warning_threshold=spec.h_alpha_warning_threshold,
        h_gamma_warning_threshold=spec.h_gamma_warning_threshold,
        fst_warning_threshold=spec.fst_warning_threshold,
        allele_loss_threshold=spec.allele_loss_threshold,
    )
    final = result.snapshots[-1]
    return BranchLockedOutcome(
        scenario_id=scenario.scenario_id,
        branch_id=branch_id,
        summary=summary,
        terminal_interaction_mean=_mean(final.interaction),
        terminal_local_effective_size_mean=_mean(final.effective_size),
        terminal_allele_frequency_mean=_mean(final.allele_frequency),
        terminal_potential_high_trait_viable=any(item.high_trait_component_present for item in final.trait_space),
        h2_comparisons=_h2_comparisons(summary),
    )


def _h2_comparisons(summary: ReplicateSummary) -> tuple[H2FirstPassageComparison, ...]:
    trait_event = _event(summary.events, "tau_trait_realised")
    return tuple(
        _h2_comparison(
            warning_id,
            getattr(summary, attribute),
            _event(summary.events, attribute),
            summary.tau_trait_realised,
            trait_event,
        )
        for warning_id, attribute in _H2_WARNING_ATTRIBUTES
    )


def _h2_comparison(
    warning_id: str,
    warning_time: int | None,
    warning_event: FirstPassageEvent | None,
    trait_loss_time: int | None,
    trait_event: FirstPassageEvent | None,
) -> H2FirstPassageComparison:
    valid = warning_time is not None and trait_loss_time is not None
    return H2FirstPassageComparison(
        warning_id=warning_id,
        warning_time=warning_time,
        trait_loss_time=trait_loss_time,
        warning_threshold=None if warning_event is None else warning_event.threshold,
        warning_aggregation_rule=None if warning_event is None else warning_event.aggregation_rule,
        trait_threshold=None if trait_event is None else trait_event.threshold,
        trait_aggregation_rule=None if trait_event is None else trait_event.aggregation_rule,
        trait_loss_observed=trait_loss_time is not None,
        warning_observed=warning_time is not None,
        valid_pair=valid,
        censored=not valid,
        warning_leads=None if not valid else warning_time < trait_loss_time,
        lead_time_trait_minus_warning=None if not valid else trait_loss_time - warning_time,
    )


def _migration_contrast(migrating: BranchLockedOutcome, isolated: BranchLockedOutcome) -> MigrationMixingContrast:
    fst_difference = None if migrating.summary.fst is None or isolated.summary.fst is None else migrating.summary.fst - isolated.summary.fst
    return MigrationMixingContrast(
        interaction_difference_migrating_minus_isolated=(
            migrating.terminal_interaction_mean - isolated.terminal_interaction_mean
        ),
        local_effective_size_difference_migrating_minus_isolated=(
            migrating.terminal_local_effective_size_mean - isolated.terminal_local_effective_size_mean
        ),
        realised_high_trait_mass_difference_migrating_minus_isolated=(
            migrating.summary.realised_high_trait_mass_mean - isolated.summary.realised_high_trait_mass_mean
        ),
        allele_frequency_mean_difference_migrating_minus_isolated=(
            migrating.terminal_allele_frequency_mean - isolated.terminal_allele_frequency_mean
        ),
        fst_difference_migrating_minus_isolated=fst_difference,
        fst_lower_with_migration=None if fst_difference is None else fst_difference < 0.0,
        h_alpha_difference_migrating_minus_isolated=migrating.summary.h_alpha - isolated.summary.h_alpha,
        h_gamma_difference_migrating_minus_isolated=migrating.summary.h_gamma - isolated.summary.h_gamma,
    )


def _summarise_cell(
    records: Sequence[BranchLockedChainReplicate],
    master_seeds: Sequence[int],
) -> dict[str, object]:
    if not records:
        raise ValueError("records must be nonempty")
    locked = tuple(record for record in records if record.h1_branch_lock_supported is True)
    h1_calibrated = tuple(record for record in records if record.h1_resolution_stable_mechanism_supported is True)
    h3_records = tuple(record for record in locked if record.high_branch_h3 is not None)
    summary = {
        "denominators": {
            "total_seed_replicates": len(records),
            "h1_resolution_stable_mechanism_count": len(h1_calibrated),
            "h1_resolution_stable_mechanism_probability": len(h1_calibrated) / len(records),
            "h1_branch_lock_count": len(locked),
            "h1_branch_lock_probability": len(locked) / len(records),
            "h1_calibration_unavailable_count": len(records) - len(h1_calibrated),
            "branch_lock_failed_after_calibration_count": len(h1_calibrated) - len(locked),
        },
        "h2_isolated_high": {
            warning_id: _summarise_h2_warning(locked, warning_id)
            for warning_id, _attribute in _H2_WARNING_ATTRIBUTES
        },
        "h3_high_branch_fragmentation": _summarise_h3(h3_records),
        "migration_as_allele_frequency_mixing": _summarise_migration(h3_records),
        "same_replicate_chain": _summarise_chain(records, locked),
        "by_master_seed": {
            str(master_seed): _summarise_seed(tuple(record for record in records if record.master_seed == master_seed))
            for master_seed in master_seeds
        },
    }
    return summary


def _summarise_h2_warning(
    locked: Sequence[BranchLockedChainReplicate],
    warning_id: str,
) -> dict[str, object]:
    comparisons = tuple(
        comparison
        for record in locked
        for comparison in _isolated_high(record).h2_comparisons
        if comparison.warning_id == warning_id
    )
    valid = tuple(value for value in comparisons if value.valid_pair)
    trait_loss = tuple(value for value in comparisons if value.trait_loss_observed)
    lead_count = sum(value.warning_leads is True for value in valid)
    return {
        "h1_branch_locked_replicates": len(locked),
        "trait_loss_observed_count": len(trait_loss),
        "trait_loss_censored_count": len(comparisons) - len(trait_loss),
        "valid_warning_trait_pair_count": len(valid),
        "censored_warning_or_trait_pair_count": sum(value.censored for value in comparisons),
        "lead_count": lead_count,
        "lead_probability_conditional_on_valid_pair": None if not valid else lead_count / len(valid),
        "lead_probability_across_h1_branch_locked_replicates": None if not comparisons else lead_count / len(comparisons),
        "lead_time_trait_minus_warning": _summary(
            value.lead_time_trait_minus_warning
            for value in valid
            if value.lead_time_trait_minus_warning is not None
        ),
    }


def _summarise_h3(records: Sequence[BranchLockedChainReplicate]) -> dict[str, object]:
    if not records:
        return {"branch_locked_replicates": 0, "fragmentation_pattern_supported_probability": None}
    values = tuple(_required(record.high_branch_h3) for record in records)
    return {
        "branch_locked_replicates": len(values),
        "fragmentation_pattern_supported_probability": _probability(value.fragmentation_pattern_supported for value in values),
        "interaction_difference_isolated_minus_large": _summary(
            value.interaction_difference_isolated_minus_large for value in values
        ),
        "local_effective_size_difference_isolated_minus_large": _summary(
            value.local_effective_size_difference_isolated_minus_large for value in values
        ),
        "realised_high_trait_mass_difference_isolated_minus_large": _summary(
            value.realised_high_trait_mass_difference_isolated_minus_large for value in values
        ),
        "h_alpha_difference_isolated_minus_large": _summary(
            value.h_alpha_difference_isolated_minus_large for value in values
        ),
        "h_gamma_difference_isolated_minus_large": _summary(
            value.h_gamma_difference_isolated_minus_large for value in values
        ),
    }


def _summarise_migration(records: Sequence[BranchLockedChainReplicate]) -> dict[str, object]:
    if not records:
        return {"branch_locked_replicates": 0, "fst_lower_with_migration_probability": None}
    values = tuple(_required(record.migration_mixing) for record in records)
    fst = tuple(value.fst_lower_with_migration for value in values if value.fst_lower_with_migration is not None)
    return {
        "branch_locked_replicates": len(values),
        "interaction_difference_migrating_minus_isolated": _summary(
            value.interaction_difference_migrating_minus_isolated for value in values
        ),
        "local_effective_size_difference_migrating_minus_isolated": _summary(
            value.local_effective_size_difference_migrating_minus_isolated for value in values
        ),
        "realised_high_trait_mass_difference_migrating_minus_isolated": _summary(
            value.realised_high_trait_mass_difference_migrating_minus_isolated for value in values
        ),
        "allele_frequency_mean_difference_migrating_minus_isolated": _summary(
            value.allele_frequency_mean_difference_migrating_minus_isolated for value in values
        ),
        "fst_difference_migrating_minus_isolated": _summary(
            value.fst_difference_migrating_minus_isolated for value in values if value.fst_difference_migrating_minus_isolated is not None
        ),
        "fst_lower_with_migration_probability": None if not fst else _probability(fst),
        "h_alpha_difference_migrating_minus_isolated": _summary(
            value.h_alpha_difference_migrating_minus_isolated for value in values
        ),
        "h_gamma_difference_migrating_minus_isolated": _summary(
            value.h_gamma_difference_migrating_minus_isolated for value in values
        ),
        "interpretation": "migration is allele-frequency mixing in this closure; no demographic rescue or recolonisation is modelled",
    }


def _summarise_chain(
    records: Sequence[BranchLockedChainReplicate],
    locked: Sequence[BranchLockedChainReplicate],
) -> dict[str, object]:
    def _support(attribute: str) -> dict[str, object]:
        observed = tuple(getattr(record, attribute) for record in records)
        applicable = tuple(value for value in observed if value is not None)
        return {
            "support_count": sum(value is True for value in applicable),
            "support_probability_across_all_seed_replicates": sum(value is True for value in applicable) / len(records),
            "support_probability_conditional_on_h1_branch_lock": None
            if not locked
            else sum(getattr(record, attribute) is True for record in locked) / len(locked),
            "unavailable_h1_calibration_count": sum(value is None for value in observed),
        }
    return {
        "h_alpha": _support("same_replicate_chain_h_alpha_supported"),
        "h_gamma": _support("same_replicate_chain_h_gamma_supported"),
        "any_predeclared_genetic_warning": _support("same_replicate_chain_any_warning_supported"),
    }


def _summarise_seed(records: Sequence[BranchLockedChainReplicate]) -> dict[str, object]:
    if not records:
        return {"replicate_count": 0}
    locked = tuple(record for record in records if record.h1_branch_lock_supported is True)
    return {
        "replicate_count": len(records),
        "h1_branch_lock_probability": len(locked) / len(records),
        "same_replicate_chain_h_alpha_probability": _probability_or_none(
            record.same_replicate_chain_h_alpha_supported for record in records
        ),
        "same_replicate_chain_h_gamma_probability": _probability_or_none(
            record.same_replicate_chain_h_gamma_supported for record in records
        ),
        "same_replicate_chain_any_warning_probability": _probability_or_none(
            record.same_replicate_chain_any_warning_supported for record in records
        ),
    }


def _scenario_map(spec: ExperimentSpec) -> Mapping[str, LandscapeScenario]:
    values = tuple(default_scenarios(spec))
    by_id = {scenario.scenario_id: scenario for scenario in values}
    if tuple(by_id) != _SCENARIOS:
        missing = sorted(set(_SCENARIOS).difference(by_id))
        extra = sorted(set(by_id).difference(_SCENARIOS))
        raise ValueError(f"unexpected default H3 scenarios; missing={missing}, extra={extra}")
    return by_id


def _calibration_identity(cell: FiniteH1BoundaryResolutionCell) -> str:
    return json.dumps(
        {
            "pair_index": cell.pair_index,
            "parameters": asdict(cell.parameters),
            "interval": cell.canonical_bistable_barrier_interval,
            "padding": cell.endpoint_padding_fraction,
            "stage_generations": cell.stage_generations,
            "grids": cell.nested_barrier_points,
        },
        sort_keys=True,
    )


def _validate_master_seeds(values: Sequence[int]) -> tuple[int, ...]:
    seeds = tuple(int(value) for value in values)
    if len(seeds) < 2:
        raise ValueError("master_seeds must contain at least two independent seeds")
    if len(set(seeds)) != len(seeds):
        raise ValueError("master_seeds must be distinct")
    if any(seed < 0 for seed in seeds):
        raise ValueError("master_seeds must be non-negative")
    return seeds


def _isolated_high(record: BranchLockedChainReplicate) -> BranchLockedOutcome:
    return _required(_required(record.outcomes)[SCENARIO_EQUAL_ISOLATED]["high_start"])


def _event(events: Sequence[FirstPassageEvent], name: str) -> FirstPassageEvent | None:
    return next((event for event in events if event.name == name), None)


def _mean(values: Sequence[float]) -> float:
    if not values:
        raise ValueError("values must be nonempty")
    return sum(float(value) for value in values) / len(values)


def _probability(values: Iterable[bool]) -> float:
    observed = tuple(values)
    if not observed:
        raise ValueError("values must be nonempty")
    return sum(observed) / len(observed)


def _probability_or_none(values: Iterable[bool | None]) -> float | None:
    observed = tuple(value for value in values if value is not None)
    return None if not observed else sum(value is True for value in observed) / len(observed)


def _summary(values: Iterable[float | int | None]) -> dict[str, float | None]:
    observed = tuple(float(value) for value in values if value is not None)
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
        raise RuntimeError("unexpected missing branch-locked chain value")
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
