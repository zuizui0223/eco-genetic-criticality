import csv
import json
from dataclasses import replace

import pytest

from causal_model.canonical_h1_bifurcation import canonical_bistable_barrier_interval
from causal_model.h1_targeted_validation_campaign import (
    build_h1_targeted_design,
    run_h1_targeted_validation_campaign,
    targeted_output_paths,
    write_h1_targeted_validation_campaign,
)
from causal_model.h1_targeted_validation_campaign_cli import main
from causal_model.multipatch_criticality_experiments import quick_profile


def _spec(**overrides):
    values = dict(
        total_area=4.0,
        patch_count=4,
        generations=2,
        replicates=1,
        master_seed=83,
        area_reference_values=(1.0,),
        interaction_feedback_values=(3.5,),
        interaction_barrier_values=(0.4,),
        migration_rate=0.2,
    )
    values.update(overrides)
    return replace(quick_profile(), **values)


def test_targeted_design_places_every_barrier_strictly_inside_its_one_large_interval():
    spec = _spec(
        area_reference_values=(0.8, 1.0),
        interaction_feedback_values=(3.5, 4.5),
    )
    design = build_h1_targeted_design(spec, normalized_positions=(0.25, 0.5, 0.75))

    assert len(design) == 12
    assert [cell.design_cell_index for cell in design] == list(range(12))
    pair_seeds = {}
    for cell in design:
        interval = canonical_bistable_barrier_interval(
            cell.interaction_feedback,
            spec.total_area,
            cell.area_reference,
        )
        assert interval is not None
        lower, upper = interval
        assert lower == cell.bistable_barrier_lower
        assert upper == cell.bistable_barrier_upper
        assert lower < cell.interaction_barrier < upper
        assert cell.normalized_barrier_position == pytest.approx(
            (cell.interaction_barrier - lower) / (upper - lower)
        )
        pair_seeds.setdefault(cell.pair_index, cell.subcampaign_seed)
        assert pair_seeds[cell.pair_index] == cell.subcampaign_seed
    assert len(set(pair_seeds.values())) == 4


def test_targeted_campaign_retains_inside_design_and_all_subcampaign_artifacts(tmp_path):
    campaign = run_h1_targeted_validation_campaign(
        _spec(),
        normalized_positions=(0.5,),
        interaction_separation_threshold=0.0,
        terminal_window=1,
        hysteresis_barrier_points=3,
        hysteresis_stage_generations=1,
    )

    assert campaign.design_cell_count == 1
    assert len(campaign.subcampaigns) == 1
    assert len(campaign.ledger) == 1
    assert campaign.ledger[0]["design"]["normalized_barrier_position"] == 0.5
    assert "finite_validation" in campaign.ledger[0]

    paths = write_h1_targeted_validation_campaign(
        campaign,
        output_dir=tmp_path,
        prefix="targeted",
        audit_arguments={"terminal_window": 1},
        code_revision="test-revision",
    )
    assert all(path.exists() for path in paths.as_dict().values())
    assert (paths.subcampaign_root / "pair_000" / "pair_000.h1_branch.csv").exists()
    with paths.ledger_csv.open(encoding="utf-8", newline="") as handle:
        rows = tuple(csv.DictReader(handle))
    records = json.loads(paths.ledger_json.read_text(encoding="utf-8"))
    manifest = json.loads(paths.manifest_json.read_text(encoding="utf-8"))
    assert len(rows) == len(records) == 1
    assert "design.normalized_barrier_position" in rows[0]
    assert manifest["campaign"] == "one_large_canonical_h1_targeted_v1"
    assert manifest["code_revision"] == "test-revision"
    assert "No finite branch outcome" in manifest["selection_policy"]


def test_targeted_design_rejects_boundary_positions_and_nonbistable_one_large_pairs():
    with pytest.raises(ValueError, match="strictly inside"):
        build_h1_targeted_design(_spec(), normalized_positions=(0.0,))
    with pytest.raises(ValueError, match="no strict bistable interval"):
        build_h1_targeted_design(_spec(interaction_feedback_values=(0.5,)))


def test_targeted_cli_records_normalized_positions_and_top_level_outputs(tmp_path):
    exit_code = main(
        [
            "--profile",
            "quick",
            "--output-dir",
            str(tmp_path),
            "--prefix",
            "cli_targeted",
            "--replicates",
            "1",
            "--generations",
            "1",
            "--inside-position",
            "0.5",
            "--terminal-window",
            "1",
            "--hysteresis-barrier-points",
            "3",
            "--hysteresis-stage-generations",
            "1",
            "--interaction-separation-threshold",
            "0.0",
        ]
    )

    assert exit_code == 0
    manifest = json.loads((tmp_path / "cli_targeted.manifest.json").read_text(encoding="utf-8"))
    assert manifest["normalized_positions"] == [0.5]
    assert manifest["design_cell_count"] == 1
    assert manifest["subcampaign_count"] == 1
    assert (tmp_path / "cli_targeted.subcampaigns" / "pair_000" / "pair_000.ledger.json").exists()


def test_targeted_output_paths_are_stable(tmp_path):
    paths = targeted_output_paths(tmp_path, "h1")
    assert paths.design_json.name == "h1.design.json"
    assert paths.ledger_csv.name == "h1.ledger.csv"
    assert paths.subcampaign_root.name == "h1.subcampaigns"
