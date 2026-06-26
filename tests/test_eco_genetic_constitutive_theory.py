from causal_model.eco_genetic_constitutive_theory import (
    MODEL_INSTANTIATIONS,
    model_by_identifier,
    proposition_ids,
    validate_constitutive_registry,
)


def test_constitutive_registry_is_coherent() -> None:
    assert validate_constitutive_registry() == ()
    assert len(proposition_ids()) == len(set(proposition_ids()))


def test_registry_separates_claim_layers() -> None:
    canonical = model_by_identifier("canonical_logistic_criticality")
    finite_bin = model_by_identifier("finite_bin_coupled_feedback")
    lead = model_by_identifier("conditional_h_alpha_lead")

    assert canonical.strongest_claim == "T"
    assert finite_bin.strongest_claim == "S"
    assert lead.strongest_claim == "C"
    assert "universal theorem" in finite_bin.excludes
    assert "universal lead ordering" in lead.excludes


def test_every_registered_model_declares_domain_and_scope_limit() -> None:
    assert MODEL_INSTANTIATIONS
    for model in MODEL_INSTANTIATIONS:
        assert model.domain
        assert model.propositions
        assert model.excludes
        assert model.supports


def test_hypothesis_coverage_is_explicit() -> None:
    coverage = {hypothesis for model in MODEL_INSTANTIATIONS for hypothesis in model.supports}
    assert coverage == {"H1", "H2", "H3"}
