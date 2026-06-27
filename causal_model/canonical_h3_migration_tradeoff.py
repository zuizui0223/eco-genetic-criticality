"""A two-patch theorem showing that rescue and genetic homogenisation can coexist.

The genetic layer is the symmetric frequency-mixing map
p1'=(1-m)p1 + m p2 and p2'=m p1 + (1-m)p2.
For 0<m<1/2, non-identical patch frequencies become strictly closer after one
step.  A separate thresholded arrival certificate records demographic rescue;
this avoids claiming that allele-frequency mixing alone rescues census size.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from causal_model.canonical_h2_h3_certificates import CanonicalH3RescueCertificate


@dataclass(frozen=True)
class CanonicalH3MigrationTradeoffCertificate:
    """Exact two-patch frequency homogenisation plus declared rescue status."""

    migration_rate: float
    initial_frequency_patch_1: float
    initial_frequency_patch_2: float
    post_migration_frequency_patch_1: float
    post_migration_frequency_patch_2: float
    initial_frequency_difference: float
    post_migration_frequency_difference: float
    mean_frequency_conserved: bool
    strict_allelic_homogenisation_certified: bool
    demographic_rescue_certified: bool
    high_trait_rescue_certified: bool
    rescue_homogenisation_tradeoff_certified: bool


def _frequency(value: float, name: str) -> float:
    value = float(value)
    if not isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return value


def canonical_h3_migration_tradeoff_certificate(
    *,
    frequency_patch_1: float,
    frequency_patch_2: float,
    migration_rate: float,
    rescue: CanonicalH3RescueCertificate | None = None,
    tolerance: float = 1e-12,
) -> CanonicalH3MigrationTradeoffCertificate:
    """Certify when symmetric migration strictly reduces allele-frequency contrast.

    If ``p1 != p2`` and ``0 < m < 1/2``, then

    ``|p1' - p2'| = (1 - 2m) |p1 - p2| < |p1 - p2|``.

    Mean allele frequency is conserved exactly.  When a separate external-arrival
    rescue certificate is supplied, the result can therefore show the canonical
    coexistence of demographic/high-trait rescue and allelic homogenisation.
    """
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative")
    p1 = _frequency(frequency_patch_1, "frequency_patch_1")
    p2 = _frequency(frequency_patch_2, "frequency_patch_2")
    m = _frequency(migration_rate, "migration_rate")
    next_1 = (1.0 - m) * p1 + m * p2
    next_2 = m * p1 + (1.0 - m) * p2
    initial_difference = abs(p1 - p2)
    final_difference = abs(next_1 - next_2)
    conserved = abs((p1 + p2) - (next_1 + next_2)) <= tolerance
    homogenises = initial_difference > tolerance and 0.0 < m < 0.5 and final_difference < initial_difference - tolerance
    demographic_rescue = bool(rescue and rescue.rescue_certified)
    high_trait_rescue = bool(rescue and rescue.high_trait_rescue_certified)
    return CanonicalH3MigrationTradeoffCertificate(
        migration_rate=m,
        initial_frequency_patch_1=p1,
        initial_frequency_patch_2=p2,
        post_migration_frequency_patch_1=next_1,
        post_migration_frequency_patch_2=next_2,
        initial_frequency_difference=initial_difference,
        post_migration_frequency_difference=final_difference,
        mean_frequency_conserved=conserved,
        strict_allelic_homogenisation_certified=homogenises,
        demographic_rescue_certified=demographic_rescue,
        high_trait_rescue_certified=high_trait_rescue,
        rescue_homogenisation_tradeoff_certified=homogenises and (demographic_rescue or high_trait_rescue),
    )
