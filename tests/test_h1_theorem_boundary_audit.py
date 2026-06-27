import pytest

from causal_model.canonical_full_simulator_bridge import canonical_full_simulator_parameters
from causal_model.h1_theorem_boundary_audit import audit_h1_theorem_boundary
from causal_model.multipatch_criticality_dynamics import DynamicsParameters, simulate


def test_audit_certifies_the_declared_single_patch_canonical_limit():
    result = simulate(
        canonical_full_simulator_parameters(
            area=1.0,
            area_reference=1.0,
            feedback_strength=8.0,
            barrier=0.5,
            initial_interaction=0.01,
            generations=12,
            carrying_population=50,
        )
    )
    audit = audit_h1_theorem_boundary(result)

    assert audit.patchwise_canonical_update_certified
    assert audit.single_patch_canonical_theorem_limit_certified
    assert audit.maximum_canonical_update_residual == pytest.approx(0.0, abs=1e-15)
    assert audit.departure_labels == ()


def test_audit_detects_density_departure_from_the_canonical_map():
    parameters = DynamicsParameters(
        patch_areas=(1.0,),
        generations=3,
        initial_population=(5,),
        initial_interaction=(0.5,),
        density_capacity=20.0,
        interaction_feedback=8.0,
        interaction_barrier=0.5,
        q_feedback_alpha=1.0,
        q_feedback_beta_trait=0.0,
        q_feedback_gamma_allele=0.0,
        baseline_growth=0.0,
        interaction_growth=0.0,
        high_allele_growth=0.0,
    )
    audit = audit_h1_theorem_boundary(simulate(parameters))

    assert not audit.patchwise_canonical_update_certified
    assert audit.maximum_density_deviation_from_one > 0.0
    assert audit.maximum_canonical_update_residual > 0.0
    assert "density_not_one" in audit.departure_labels


def test_audit_names_trait_allele_migration_and_multiple_patch_departures():
    parameters = DynamicsParameters(
        patch_areas=(1.0, 1.0),
        generations=2,
        initial_population=(40, 40),
        initial_interaction=(0.3, 0.7),
        density_capacity=40.0,
        q_feedback_alpha=0.6,
        q_feedback_beta_trait=0.3,
        q_feedback_gamma_allele=0.1,
        migration_rate=0.2,
        baseline_growth=1.0,
        interaction_growth=0.0,
        high_allele_growth=0.0,
    )
    audit = audit_h1_theorem_boundary(simulate(parameters))

    assert "trait_feedback_enabled" in audit.departure_labels
    assert "allele_feedback_enabled" in audit.departure_labels
    assert "migration_enabled" in audit.departure_labels
    assert "multiple_patches" in audit.departure_labels
    assert not audit.single_patch_canonical_theorem_limit_certified
