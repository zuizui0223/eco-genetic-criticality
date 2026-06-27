import csv
import json

import pytest

from causal_model.h3_landscape_presets import (
    equal_complete_network,
    equal_isolated,
    one_large,
)
from causal_model.h3_phase_diagram import (
    H3PhaseDiagramSpec,
    initial_states_for_preset,
    run_h3_phase_diagram,
    standard_h3_landscapes,
    write_h3_phase_diagram_artifacts,
)


def _spec(**overrides):
    values = {
        "total_capacity": 40,
        "patch_count": 4,
        "generations": 2,
        "replicates": 2,
        "adult_survival_grid": (1.0,),
        "emigration_grid": (0.0, 1.0),
        "recruitment_per_adult": 0.0,
        "persistence_threshold": 1,
        "colonisation_threshold": 1,
        "initial_occupancy_fraction": 0.5,
        "initial_high_trait_fraction": 0.5,
        "initial_high_allele_frequency": 0.5,
        "random_seed": 5,
    }
    values.update(overrides)
    return H3PhaseDiagramSpec(**values)


def test_presets_keep_total_capacity_matched_and_distinguish_connectivity():
    large = one_large(40)
    isolated = equal_isolated(40, 4)
    connected = equal_complete_network(40, 4)

    assert [preset.total_capacity for preset in (large, isolated, connected)] == [40, 40, 40]
    assert isolated.kernel[0] == (1.0, 0.0, 0.0, 0.0)
    assert connected.kernel[0][0] == 0.0
    assert connected.kernel[0][1:] == pytest.approx((1 / 3, 1 / 3, 1 / 3))


def test_initial_states_have_matched_composition_proportional_to_capacity():
    spec = _spec()
    large, isolated, _ = standard_h3_landscapes(spec)
    large_state = initial_states_for_preset(large, spec)
    isolated_states = initial_states_for_preset(isolated, spec)

    assert large_state == (type(large_state[0])(20, 10, 20),)
    assert sum(state.population for state in isolated_states) == 20
    assert sum(state.high_trait_abundance for state in isolated_states) == 8
    assert sum(state.high_allele_copies for state in isolated_states) == 16


def test_phase_diagram_has_all_landscape_and_grid_cells_with_visible_denominators():
    spec = _spec()
    rows = run_h3_phase_diagram(spec)

    assert len(rows) == 3 * 1 * 2
    assert {row.scenario_id for row in rows} == {
        "one_large",
        "equal_isolated",
        "equal_complete_network",
    }
    assert all(row.total_capacity == 40 for row in rows)
    assert all(row.replicates == 2 for row in rows)
    assert all(0.0 <= row.metapopulation_extinction_probability <= 1.0 for row in rows)
    assert all(0.0 <= row.realised_high_trait_loss_probability <= 1.0 for row in rows)


def test_phase_diagram_artifacts_are_flat_csv_and_json(tmp_path):
    rows = run_h3_phase_diagram(_spec(emigration_grid=(0.0,)))
    csv_path = tmp_path / "h3.csv"
    json_path = tmp_path / "h3.json"
    write_h3_phase_diagram_artifacts(rows, csv_path=csv_path, json_path=json_path)

    with csv_path.open(encoding="utf-8", newline="") as handle:
        table = tuple(csv.DictReader(handle))
    with json_path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    assert len(table) == len(rows) == len(payload)
    assert table[0]["scenario_id"]
    assert "recolonisation_probability" in payload[0]


def test_invalid_capacity_partition_or_unmatched_landscape_is_rejected():
    with pytest.raises(ValueError):
        _spec(total_capacity=41)
    with pytest.raises(ValueError):
        run_h3_phase_diagram(
            _spec(),
            landscapes=(one_large(20),),
        )
