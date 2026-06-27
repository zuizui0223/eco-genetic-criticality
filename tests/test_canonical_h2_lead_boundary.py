import pytest

from causal_model.canonical_h2_lead_boundary import canonical_h2_lead_boundary
from causal_model.canonical_h2_h3_certificates import canonical_h2_expected_lead_certificate


def test_closed_form_h2_boundary_matches_iterated_first_passages():
    kwargs = {
        "initial_heterozygosity": 1.0,
        "initial_trait_abundance": 10.0,
        "effective_population_size": 2,
        "trait_retention": 0.9,
        "heterozygosity_threshold": 0.5,
        "trait_threshold": 5.0,
    }
    boundary = canonical_h2_lead_boundary(**kwargs)
    iterated = canonical_h2_expected_lead_certificate(**kwargs, horizon=10)

    assert boundary.discrete_genetic_warning_time == iterated.genetic_warning_time == 3
    assert boundary.discrete_trait_loss_time == iterated.trait_loss_time == 7
    assert boundary.strict_expected_genetic_lead_certified


def test_unit_trait_retention_keeps_trait_event_unobserved_in_closed_form():
    boundary = canonical_h2_lead_boundary(
        initial_heterozygosity=1.0,
        initial_trait_abundance=10.0,
        effective_population_size=2,
        trait_retention=1.0,
        heterozygosity_threshold=0.5,
        trait_threshold=5.0,
    )

    assert boundary.discrete_trait_loss_time is None
    assert not boundary.strict_expected_genetic_lead_certified


def test_boundary_rejects_invalid_effective_size_or_threshold_order():
    with pytest.raises(ValueError):
        canonical_h2_lead_boundary(
            initial_heterozygosity=1.0,
            initial_trait_abundance=10.0,
            effective_population_size=0.5,
            trait_retention=0.9,
            heterozygosity_threshold=0.5,
            trait_threshold=5.0,
        )
    with pytest.raises(ValueError):
        canonical_h2_lead_boundary(
            initial_heterozygosity=0.5,
            initial_trait_abundance=10.0,
            effective_population_size=2,
            trait_retention=0.9,
            heterozygosity_threshold=0.5,
            trait_threshold=5.0,
        )
