import pytest

from causal_model.mutation_h1_primary_domain import (
    PRIMARY_MINIMUM_SEED_BLOCK_SUPPORT,
    VALIDATED_CELLS,
    domain_manifest,
    find_validated_cell,
    is_primary_analysis_cell,
    primary_analysis_cells,
    validate_primary_domain,
)


def test_primary_domain_retains_complete_validation_ledger_and_exactly_twelve_cells():
    validate_primary_domain()
    assert len(VALIDATED_CELLS) == 27
    assert len(primary_analysis_cells()) == 12
    assert domain_manifest()["counts"] == {
        "validated_cells": 27,
        "primary_analysis_cells": 12,
        "excluded_cells": 15,
    }


def test_all_seed_rule_accepts_boundary_cell_and_rejects_high_pooled_but_heterogeneous_cell():
    boundary = find_validated_cell(0.10, 0.8, 6.0)
    heterogeneous = find_validated_cell(0.20, 1.0, 3.0)
    assert boundary.minimum_seed_block_support == PRIMARY_MINIMUM_SEED_BLOCK_SUPPORT
    assert boundary.primary_analysis_eligible is True
    assert heterogeneous.pooled_joint_support == 0.80
    assert heterogeneous.minimum_seed_block_support == 0.60
    assert heterogeneous.primary_analysis_eligible is False


def test_unvalidated_configuration_is_rejected_instead_of_silently_selected():
    with pytest.raises(ValueError, match="not part of independent_validation_v1"):
        is_primary_analysis_cell(0.05, 0.8, 6.0)
