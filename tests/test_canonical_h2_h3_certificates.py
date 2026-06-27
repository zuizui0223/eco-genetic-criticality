import pytest

from causal_model.canonical_h2_h3_certificates import (
    canonical_h2_expected_lead_certificate,
    canonical_h3_fragmentation_certificate,
    canonical_h3_rescue_certificate,
)


def test_h2_certifies_strict_expected_genetic_lead_when_drift_crosses_first():
    certificate = canonical_h2_expected_lead_certificate(
        initial_heterozygosity=1.0,
        initial_trait_abundance=10.0,
        effective_population_size=2,
        trait_retention=0.9,
        heterozygosity_threshold=0.5,
        trait_threshold=5.0,
        horizon=10,
    )

    assert certificate.genetic_warning_time == 3
    assert certificate.trait_loss_time == 7
    assert certificate.strict_expected_genetic_lead_certified


def test_h2_keeps_censored_event_none_and_does_not_claim_lead():
    certificate = canonical_h2_expected_lead_certificate(
        initial_heterozygosity=1.0,
        initial_trait_abundance=10.0,
        effective_population_size=1000,
        trait_retention=1.0,
        heterozygosity_threshold=0.5,
        trait_threshold=5.0,
        horizon=10,
    )

    assert certificate.genetic_warning_time is None
    assert certificate.trait_loss_time is None
    assert not certificate.strict_expected_genetic_lead_certified


def test_h3_certifies_fixed_total_area_loss_of_local_support_after_fragmentation():
    certificate = canonical_h3_fragmentation_certificate(
        total_area=100.0,
        patch_count=4,
        local_support_threshold=30.0,
    )

    assert certificate.one_large_supports_high_trait
    assert certificate.equal_fragment_area == 25.0
    assert certificate.every_equal_fragment_lacks_support
    assert certificate.fragmentation_removes_local_support_certified


def test_h3_does_not_claim_fragmentation_loss_when_each_fragment_remains_supported():
    certificate = canonical_h3_fragmentation_certificate(
        total_area=100.0,
        patch_count=4,
        local_support_threshold=25.0,
    )

    assert certificate.one_large_supports_high_trait
    assert not certificate.every_equal_fragment_lacks_support
    assert not certificate.fragmentation_removes_local_support_certified


def test_h3_certifies_external_demographic_and_high_trait_rescue_only_when_thresholds_cross():
    certificate = canonical_h3_rescue_certificate(
        residents_after_survival=1,
        external_immigrants=3,
        resident_high_trait_individuals=0,
        immigrant_high_trait_individuals=2,
        establishment_threshold=4,
        high_trait_establishment_threshold=2,
        local_interaction_support=True,
    )

    assert certificate.demographic_establishment_certified
    assert certificate.rescue_certified
    assert not certificate.recolonisation_certified
    assert certificate.high_trait_rescue_certified


def test_h3_recolonisation_requires_external_arrivals_and_local_support_for_high_trait():
    certificate = canonical_h3_rescue_certificate(
        residents_after_survival=0,
        external_immigrants=4,
        resident_high_trait_individuals=0,
        immigrant_high_trait_individuals=2,
        establishment_threshold=4,
        high_trait_establishment_threshold=2,
        local_interaction_support=False,
    )

    assert certificate.recolonisation_certified
    assert not certificate.rescue_certified
    assert not certificate.high_trait_rescue_certified


def test_h2_and_h3_reject_incoherent_inputs():
    with pytest.raises(ValueError):
        canonical_h2_expected_lead_certificate(
            initial_heterozygosity=0.5,
            initial_trait_abundance=10.0,
            effective_population_size=10,
            trait_retention=0.9,
            heterozygosity_threshold=0.5,
            trait_threshold=5.0,
            horizon=10,
        )
    with pytest.raises(ValueError):
        canonical_h3_rescue_certificate(
            residents_after_survival=1,
            external_immigrants=1,
            resident_high_trait_individuals=2,
            immigrant_high_trait_individuals=0,
            establishment_threshold=2,
            high_trait_establishment_threshold=1,
            local_interaction_support=True,
        )
