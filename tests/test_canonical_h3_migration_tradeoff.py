import pytest

from causal_model.canonical_h2_h3_certificates import canonical_h3_rescue_certificate
from causal_model.canonical_h3_migration_tradeoff import canonical_h3_migration_tradeoff_certificate


def test_symmetric_migration_strictly_homogenises_nonidentical_patch_frequencies_below_half():
    certificate = canonical_h3_migration_tradeoff_certificate(
        frequency_patch_1=0.9,
        frequency_patch_2=0.1,
        migration_rate=0.25,
    )

    assert certificate.post_migration_frequency_patch_1 == pytest.approx(0.7)
    assert certificate.post_migration_frequency_patch_2 == pytest.approx(0.3)
    assert certificate.initial_frequency_difference == pytest.approx(0.8)
    assert certificate.post_migration_frequency_difference == pytest.approx(0.4)
    assert certificate.mean_frequency_conserved
    assert certificate.strict_allelic_homogenisation_certified


def test_tradeoff_certifies_when_external_arrivals_rescue_a_recipient():
    rescue = canonical_h3_rescue_certificate(
        residents_after_survival=1,
        external_immigrants=3,
        resident_high_trait_individuals=0,
        immigrant_high_trait_individuals=2,
        establishment_threshold=4,
        high_trait_establishment_threshold=2,
        local_interaction_support=True,
    )
    certificate = canonical_h3_migration_tradeoff_certificate(
        frequency_patch_1=0.9,
        frequency_patch_2=0.1,
        migration_rate=0.25,
        rescue=rescue,
    )

    assert certificate.demographic_rescue_certified
    assert certificate.high_trait_rescue_certified
    assert certificate.rescue_homogenisation_tradeoff_certified


def test_no_strict_homogenisation_for_identical_frequencies_or_boundary_rate():
    identical = canonical_h3_migration_tradeoff_certificate(
        frequency_patch_1=0.5,
        frequency_patch_2=0.5,
        migration_rate=0.25,
    )
    boundary = canonical_h3_migration_tradeoff_certificate(
        frequency_patch_1=0.9,
        frequency_patch_2=0.1,
        migration_rate=0.5,
    )

    assert not identical.strict_allelic_homogenisation_certified
    assert not boundary.strict_allelic_homogenisation_certified
    assert boundary.post_migration_frequency_difference == pytest.approx(0.0)


def test_invalid_frequency_is_rejected():
    with pytest.raises(ValueError):
        canonical_h3_migration_tradeoff_certificate(
            frequency_patch_1=1.1,
            frequency_patch_2=0.1,
            migration_rate=0.25,
        )
