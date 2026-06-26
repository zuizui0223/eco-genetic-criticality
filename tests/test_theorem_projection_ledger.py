import pytest

from causal_model.theorem_projection_ledger import (
    projection_for,
    summary_by_status,
    theorem_projections,
    validate_projection_ledger,
)


def test_projection_ledger_is_internally_valid():
    validate_projection_ledger()
    assert tuple(item.key for item in theorem_projections()) == (
        "abstract_positive_two_factor_model",
        "colonization_one_step_recruitment_submodel",
        "spatial_pollination_abm",
        "colonization_connectivity_abm",
        "defense_metapopulation_abm",
        "campanula_published_record",
    )


def test_only_declared_factor_models_are_marked_theorem_exact():
    summary = summary_by_status()
    assert summary["exact"] == (
        "abstract_positive_two_factor_model",
        "colonization_one_step_recruitment_submodel",
    )
    assert set(summary["requires_factorization_extension"]) == {
        "spatial_pollination_abm",
        "colonization_connectivity_abm",
        "defense_metapopulation_abm",
    }
    assert summary["not_applicable"] == ("campanula_published_record",)


def test_one_step_colonization_submodel_is_exact_but_multistep_backend_is_not():
    one_step = projection_for("colonization_one_step_recruitment_submodel")
    multistep = projection_for("colonization_connectivity_abm")

    assert one_step.status == "exact"
    assert one_step.theorem_ids == ("N1", "N2", "N3", "N4")
    assert "juvenile recruitment" in one_step.target
    assert "long-run invasion lambda" in one_step.prohibited_conclusion

    assert multistep.status == "requires_factorization_extension"
    assert multistep.theorem_ids == ()
    assert "one-step" in multistep.current_factorisation.lower()
    assert "one_step_to_lambda_discrepancy" in multistep.next_outputs


def test_spatial_and_defense_backends_cannot_claim_direct_channel_identification():
    spatial = projection_for("spatial_pollination_abm")
    defense = projection_for("defense_metapopulation_abm")

    for projection in (spatial, defense):
        assert projection.status == "requires_factorization_extension"
        assert projection.theorem_ids == ()

    assert "cannot" in spatial.permitted_conclusion
    assert "not an exact instance" in defense.permitted_conclusion


def test_published_campanula_record_is_explicitly_not_a_channel_identification_case():
    projection = projection_for("campanula_published_record")
    assert projection.status == "not_applicable"
    assert projection.theorem_ids == ()
    assert "does not" in projection.permitted_conclusion
    assert "calibrated_pollination_or_establishment_proxy" in projection.next_outputs


def test_unknown_projection_is_not_silently_coerced():
    with pytest.raises(KeyError, match="unknown projection"):
        projection_for("made_up_backend")
