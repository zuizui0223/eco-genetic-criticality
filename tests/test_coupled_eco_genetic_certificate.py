from causal_model.canonical_h1_bifurcation import canonical_h1_certificate
from causal_model.canonical_h2_h3_certificates import (
    canonical_h2_expected_lead_certificate,
    canonical_h3_fragmentation_certificate,
    canonical_h3_rescue_certificate,
)
from causal_model.coupled_eco_genetic_certificate import coupled_eco_genetic_certificate
from causal_model.multipatch_criticality_dynamics import DynamicsParameters


def _h1():
    return canonical_h1_certificate(
        feedback_strength=8.0,
        area=1.0,
        area_reference=1.0,
        barrier=0.5,
        trait_parameters=DynamicsParameters(patch_areas=(1.0,)),
    )


def _h2_lead():
    return canonical_h2_expected_lead_certificate(
        initial_heterozygosity=1.0,
        initial_trait_abundance=10.0,
        effective_population_size=2,
        trait_retention=0.9,
        heterozygosity_threshold=0.5,
        trait_threshold=5.0,
        horizon=10,
    )


def test_chain_certifies_only_when_each_declared_link_is_present():
    certificate = coupled_eco_genetic_certificate(
        h1=_h1(),
        h2_low_branch=_h2_lead(),
        h3_fragmentation=canonical_h3_fragmentation_certificate(
            total_area=100.0,
            patch_count=4,
            local_support_threshold=30.0,
        ),
        baseline_effective_population_size=2.0,
        interaction_slope=10.0,
    )

    assert certificate.h1_branch_switch_certified
    assert certificate.h3_fragmentation_removes_support_certified
    assert certificate.lower_effective_size_on_low_branch_certified
    assert certificate.h2_expected_genetic_lead_certified
    assert certificate.fragmentation_to_genetic_warning_chain_certified
    assert certificate.low_branch_effective_population_size < certificate.high_branch_effective_population_size


def test_chain_fails_without_a_positive_interaction_effective_size_link():
    certificate = coupled_eco_genetic_certificate(
        h1=_h1(),
        h2_low_branch=_h2_lead(),
        h3_fragmentation=canonical_h3_fragmentation_certificate(
            total_area=100.0,
            patch_count=4,
            local_support_threshold=30.0,
        ),
        baseline_effective_population_size=2.0,
        interaction_slope=0.0,
    )

    assert not certificate.lower_effective_size_on_low_branch_certified
    assert not certificate.fragmentation_to_genetic_warning_chain_certified


def test_external_rescue_marks_the_mechanism_that_interrupts_the_chain():
    certificate = coupled_eco_genetic_certificate(
        h1=_h1(),
        h2_low_branch=_h2_lead(),
        h3_fragmentation=canonical_h3_fragmentation_certificate(
            total_area=100.0,
            patch_count=4,
            local_support_threshold=30.0,
        ),
        baseline_effective_population_size=2.0,
        interaction_slope=10.0,
        h3_rescue=canonical_h3_rescue_certificate(
            residents_after_survival=1,
            external_immigrants=3,
            resident_high_trait_individuals=0,
            immigrant_high_trait_individuals=2,
            establishment_threshold=4,
            high_trait_establishment_threshold=2,
            local_interaction_support=True,
        ),
    )

    assert certificate.fragmentation_to_genetic_warning_chain_certified
    assert certificate.external_demographic_rescue_certified
    assert certificate.external_high_trait_rescue_certified
    assert certificate.rescue_interrupts_fragmentation_chain
