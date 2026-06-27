from causal_model.empirical_measurement_contract import audit_empirical_columns


def test_complete_measurement_contract_marks_all_hypotheses_ready():
    columns = (
        "patch_id",
        "patch_area",
        "interaction_state",
        "trait_value",
        "performance",
        "time",
        "realised_high_trait_abundance",
        "sample_size",
        "high_allele_copies",
        "census_population",
        "source_patch_id",
        "destination_patch_id",
        "dispersal_count",
    )
    audit = audit_empirical_columns(columns)

    assert audit.h1.ready
    assert audit.h2.ready
    assert audit.h3.ready
    assert audit.all_ready


def test_audit_keeps_hypothesis_specific_missing_measurements_visible():
    audit = audit_empirical_columns(("patch_id", "patch_area", "trait_value"))

    assert not audit.h1.ready
    assert set(audit.h1.missing_columns) == {"interaction_state", "performance"}
    assert "time" in audit.h2.missing_columns
    assert "source_patch_id" in audit.h3.missing_columns
    assert not audit.all_ready


def test_audit_normalises_blank_and_duplicate_column_names():
    audit = audit_empirical_columns((" patch_id ", "", "patch_id", "patch_area"))

    assert audit.h1.available_columns == ("patch_area", "patch_id")
