import pytest

from causal_model.network_migration_matrix_theory import (
    common_floor_certificate,
    complete_graph_mixing_matrix,
    destination_lower_bound,
    focal_rescue_certificate,
    mix_allele_frequencies,
    validate_destination_by_source_matrix,
)


def test_complete_graph_matrix_reproduces_global_mean_migration_update():
    frequencies = (0.2, 0.8)
    weights = (1.0, 3.0)
    migration_rate = 0.25
    matrix = complete_graph_mixing_matrix(weights, migration_rate)
    observed = mix_allele_frequencies(frequencies, matrix)
    mean = 0.25 * frequencies[0] + 0.75 * frequencies[1]
    expected = tuple(
        (1.0 - migration_rate) * frequency + migration_rate * mean
        for frequency in frequencies
    )
    assert observed == pytest.approx(expected)


def test_common_floor_is_preserved_by_an_asymmetric_network():
    matrix = (
        (0.7, 0.3, 0.0),
        (0.2, 0.5, 0.3),
        (0.0, 0.6, 0.4),
    )
    certificate = common_floor_certificate(matrix, source_floor=0.35)
    assert certificate.common_floor_preserved
    assert certificate.destination_lower_bounds == pytest.approx((0.35, 0.35, 0.35))


def test_focal_rescue_requires_sufficient_incoming_source_floor_and_weight():
    matrix = (
        (0.7, 0.3),
        (0.2, 0.8),
    )
    source_bounds = (0.1, 0.8)
    assert destination_lower_bound(matrix, source_bounds, destination_index=0) == pytest.approx(0.31)

    rescued = focal_rescue_certificate(
        matrix,
        source_bounds,
        destination_index=0,
        target_floor=0.30,
    )
    not_rescued = focal_rescue_certificate(
        matrix,
        source_bounds,
        destination_index=0,
        target_floor=0.32,
    )
    assert rescued.rescue_certified
    assert not not_rescued.rescue_certified


def test_invalid_matrix_or_incompatible_vector_is_rejected():
    with pytest.raises(ValueError):
        validate_destination_by_source_matrix(((0.7, 0.2), (0.3, 0.7)))
    with pytest.raises(ValueError):
        validate_destination_by_source_matrix(((1.0, 0.0, 0.0),))
    with pytest.raises(ValueError):
        mix_allele_frequencies((0.2,), ((0.5, 0.5), (0.5, 0.5)))
