"""Exact certificates for deliberately narrow H2 and H3 canonical maps.

H2 uses an expected-heterozygosity recursion with constant effective size and a
monotone realised-trait decline.  It proves an event ordering only for that
specified deterministic expectation map, not for all stochastic realisations.

H3 uses binary local interaction support and thresholded post-arrival
establishment.  It proves when equal subdivision removes local support at fixed
total area, and when declared external arrivals cross a rescue threshold.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite


@dataclass(frozen=True)
class CanonicalH2LeadCertificate:
    """First-passage ordering in the declared deterministic H2 expectation map."""

    initial_heterozygosity: float
    initial_trait_abundance: float
    effective_population_size: int
    trait_retention: float
    heterozygosity_threshold: float
    trait_threshold: float
    horizon: int
    expected_heterozygosity: tuple[float, ...]
    realised_trait_abundance: tuple[float, ...]
    genetic_warning_time: int | None
    trait_loss_time: int | None
    strict_expected_genetic_lead_certified: bool


@dataclass(frozen=True)
class CanonicalH3FragmentationCertificate:
    """Fixed-total-area local-support certificate for equal isolated fragments."""

    total_area: float
    patch_count: int
    local_support_threshold: float
    one_large_supports_high_trait: bool
    equal_fragment_area: float
    every_equal_fragment_lacks_support: bool
    fragmentation_removes_local_support_certified: bool


@dataclass(frozen=True)
class CanonicalH3RescueCertificate:
    """Threshold certificate for a recipient after external arrival."""

    residents_after_survival: int
    external_immigrants: int
    resident_high_trait_individuals: int
    immigrant_high_trait_individuals: int
    establishment_threshold: int
    high_trait_establishment_threshold: int
    local_interaction_support: bool
    post_arrival_population: int
    post_arrival_high_trait_individuals: int
    demographic_establishment_certified: bool
    recolonisation_certified: bool
    rescue_certified: bool
    high_trait_rescue_certified: bool


def _first_at_or_below(values: tuple[float, ...], threshold: float) -> int | None:
    return next((time for time, value in enumerate(values) if value <= threshold), None)


def _validate_probability(value: float, name: str, *, positive: bool = False) -> float:
    value = float(value)
    lower = 0.0 if not positive else 0.0
    if not isfinite(value) or not lower <= value <= 1.0 or (positive and value == 0.0):
        interval = "(0, 1]" if positive else "[0, 1]"
        raise ValueError(f"{name} must lie in {interval}")
    return value


def canonical_h2_expected_lead_certificate(
    *,
    initial_heterozygosity: float,
    initial_trait_abundance: float,
    effective_population_size: int,
    trait_retention: float,
    heterozygosity_threshold: float,
    trait_threshold: float,
    horizon: int,
) -> CanonicalH2LeadCertificate:
    """Certify H2 event ordering for a specified expectation map.

    The model is

    ``E[H[t+1]] = (1 - 1/(2*N_e)) E[H[t]]``
    ``T[t+1] = r_T T[t]``.

    With positive finite ``N_e`` and ``0 < r_T <= 1``, both paths are monotone.
    Their recorded first-passage times are therefore exact for this map.  A
    strict genetic lead is certified only when both events occur within the
    declared horizon and ``tau_H < tau_T``.  Censored events remain ``None``.
    """
    h0 = _validate_probability(initial_heterozygosity, "initial_heterozygosity", positive=True)
    h_star = _validate_probability(heterozygosity_threshold, "heterozygosity_threshold")
    if h_star >= h0:
        raise ValueError("heterozygosity_threshold must be below initial_heterozygosity")
    if not isfinite(initial_trait_abundance) or initial_trait_abundance <= 0.0:
        raise ValueError("initial_trait_abundance must be finite and positive")
    if not isfinite(trait_threshold) or not 0.0 < trait_threshold < initial_trait_abundance:
        raise ValueError("trait_threshold must be positive and below initial_trait_abundance")
    if not isinstance(effective_population_size, int) or effective_population_size < 1:
        raise ValueError("effective_population_size must be a positive integer")
    retention = _validate_probability(trait_retention, "trait_retention", positive=True)
    if not isinstance(horizon, int) or horizon < 1:
        raise ValueError("horizon must be a positive integer")

    heterozygosity_retention = 1.0 - 1.0 / (2.0 * effective_population_size)
    heterozygosity = [h0]
    trait = [float(initial_trait_abundance)]
    for _ in range(horizon):
        heterozygosity.append(heterozygosity[-1] * heterozygosity_retention)
        trait.append(trait[-1] * retention)
    h_path = tuple(heterozygosity)
    trait_path = tuple(trait)
    genetic_time = _first_at_or_below(h_path, h_star)
    trait_time = _first_at_or_below(trait_path, trait_threshold)
    return CanonicalH2LeadCertificate(
        initial_heterozygosity=h0,
        initial_trait_abundance=float(initial_trait_abundance),
        effective_population_size=effective_population_size,
        trait_retention=retention,
        heterozygosity_threshold=h_star,
        trait_threshold=float(trait_threshold),
        horizon=horizon,
        expected_heterozygosity=h_path,
        realised_trait_abundance=trait_path,
        genetic_warning_time=genetic_time,
        trait_loss_time=trait_time,
        strict_expected_genetic_lead_certified=(
            genetic_time is not None and trait_time is not None and genetic_time < trait_time
        ),
    )


def canonical_h3_fragmentation_certificate(
    *,
    total_area: float,
    patch_count: int,
    local_support_threshold: float,
) -> CanonicalH3FragmentationCertificate:
    """Certify loss of binary local support after equal isolated subdivision.

    In this canonical H3 closure a patch can maintain the high-trait mode exactly
    when ``patch_area >= local_support_threshold``.  Thus, when one large patch
    clears the threshold but every equal isolated fragment does not, fixed-total-
    area fragmentation removes the specified maintenance mechanism.
    """
    for name, value in (("total_area", total_area), ("local_support_threshold", local_support_threshold)):
        if not isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    if not isinstance(patch_count, int) or patch_count < 1:
        raise ValueError("patch_count must be a positive integer")
    equal_area = float(total_area) / patch_count
    one_large = float(total_area) >= float(local_support_threshold)
    fragments_lack = equal_area < float(local_support_threshold)
    return CanonicalH3FragmentationCertificate(
        total_area=float(total_area),
        patch_count=patch_count,
        local_support_threshold=float(local_support_threshold),
        one_large_supports_high_trait=one_large,
        equal_fragment_area=equal_area,
        every_equal_fragment_lacks_support=fragments_lack,
        fragmentation_removes_local_support_certified=one_large and fragments_lack,
    )


def canonical_h3_rescue_certificate(
    *,
    residents_after_survival: int,
    external_immigrants: int,
    resident_high_trait_individuals: int,
    immigrant_high_trait_individuals: int,
    establishment_threshold: int,
    high_trait_establishment_threshold: int,
    local_interaction_support: bool,
) -> CanonicalH3RescueCertificate:
    """Certify thresholded recolonisation or rescue after declared arrivals.

    The canonical post-arrival rule is: a recipient establishes when its total
    post-arrival population reaches ``establishment_threshold``.  Its high-trait
    component establishes when local interaction support is present and its
    post-arrival high-trait count reaches ``high_trait_establishment_threshold``.
    External immigrants are deliberately separate from residents, so self-loop
    transport cannot be labelled a rescue.
    """
    integers = (
        ("residents_after_survival", residents_after_survival),
        ("external_immigrants", external_immigrants),
        ("resident_high_trait_individuals", resident_high_trait_individuals),
        ("immigrant_high_trait_individuals", immigrant_high_trait_individuals),
        ("establishment_threshold", establishment_threshold),
        ("high_trait_establishment_threshold", high_trait_establishment_threshold),
    )
    for name, value in integers:
        if not isinstance(value, int) or value < 0:
            raise ValueError(f"{name} must be a non-negative integer")
    if establishment_threshold < 1 or high_trait_establishment_threshold < 1:
        raise ValueError("establishment thresholds must be positive")
    if resident_high_trait_individuals > residents_after_survival:
        raise ValueError("resident high-trait abundance cannot exceed residents")
    if immigrant_high_trait_individuals > external_immigrants:
        raise ValueError("immigrant high-trait abundance cannot exceed immigrants")

    population = residents_after_survival + external_immigrants
    high_trait = resident_high_trait_individuals + immigrant_high_trait_individuals
    establishes = population >= establishment_threshold
    recolonisation = residents_after_survival == 0 and establishes
    rescue = 0 < residents_after_survival < establishment_threshold and establishes
    high_trait_rescue = establishes and local_interaction_support and high_trait >= high_trait_establishment_threshold
    return CanonicalH3RescueCertificate(
        residents_after_survival=residents_after_survival,
        external_immigrants=external_immigrants,
        resident_high_trait_individuals=resident_high_trait_individuals,
        immigrant_high_trait_individuals=immigrant_high_trait_individuals,
        establishment_threshold=establishment_threshold,
        high_trait_establishment_threshold=high_trait_establishment_threshold,
        local_interaction_support=bool(local_interaction_support),
        post_arrival_population=population,
        post_arrival_high_trait_individuals=high_trait,
        demographic_establishment_certified=establishes,
        recolonisation_certified=recolonisation,
        rescue_certified=rescue,
        high_trait_rescue_certified=high_trait_rescue,
    )
