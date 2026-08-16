from __future__ import annotations

import argparse
import csv
import html
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Iterable

EXPECTED_CELL_COUNT = 12
EXPECTED_SEEDS_PER_CELL = 5
EXPECTED_REPLICATES_PER_SEED = 20
EXPECTED_PATCH_COUNTS = (1, 2, 3, 4, 6, 8, 12, 16)
PRIMARY_METRICS = (
    "final_interaction_mean",
    "final_effective_size_mean",
    "realised_high_trait_mass_mean",
)
METRIC_LABELS = {
    "final_interaction_mean": "Final interaction",
    "final_effective_size_mean": "Local effective size",
    "realised_high_trait_mass_mean": "Realised high-trait mass",
}


def _read_rows(root: Path) -> list[dict[str, str]]:
    paths = sorted(root.rglob("*.csv"))
    if not paths:
        raise FileNotFoundError(f"no CSV inputs under {root}")
    rows: list[dict[str, str]] = []
    for path in paths:
        with path.open(encoding="utf-8", newline="") as handle:
            rows.extend(csv.DictReader(handle))
    return rows


def _float(row: dict[str, str], key: str) -> float | None:
    value = row[key]
    return None if value in ("", "None") else float(value)


def _bool(row: dict[str, str], key: str) -> bool | None:
    value = row[key]
    if value in ("", "None"):
        return None
    return value.lower() == "true"


def _source_key(row: dict[str, str]) -> tuple[int, int, int]:
    return (int(row["primary_cell_index"]), int(row["master_seed"]), int(row["replicate_index"]))


def _quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    if len(ordered) == 1:
        return ordered[0]
    position = probability * (len(ordered) - 1)
    lo = math.floor(position)
    hi = math.ceil(position)
    if lo == hi:
        return ordered[lo]
    weight = position - lo
    return ordered[lo] * (1.0 - weight) + ordered[hi] * weight


def _median_iqr(values: Iterable[float]) -> tuple[float | None, float | None, float | None]:
    values = [float(value) for value in values if value is not None and math.isfinite(float(value))]
    if not values:
        return None, None, None
    return statistics.median(values), _quantile(values, 0.25), _quantile(values, 0.75)


def _validated_rows(rows: list[dict[str, str]]) -> tuple[list[dict[str, object]], dict[tuple[int, int, int], dict[int, dict[str, object]]]]:
    expected_rows = (
        EXPECTED_CELL_COUNT
        * EXPECTED_SEEDS_PER_CELL
        * EXPECTED_REPLICATES_PER_SEED
        * len(EXPECTED_PATCH_COUNTS)
    )
    if len(rows) != expected_rows:
        raise ValueError(f"expected {expected_rows} rows, found {len(rows)}")

    parsed: list[dict[str, object]] = []
    by_source: dict[tuple[int, int, int], dict[int, dict[str, object]]] = defaultdict(dict)
    for raw in rows:
        row: dict[str, object] = dict(raw)
        for key in ("primary_cell_index", "master_seed", "replicate_index", "calibration_seed", "patch_count"):
            row[key] = int(raw[key])
        for key in (
            "mutation_rate", "area_reference", "interaction_feedback", "patch_area",
            "theoretical_k_fragment", *PRIMARY_METRICS, "h_alpha", "h_gamma",
        ):
            row[key] = _float(raw, key)
        for key in (
            "source_prepared", "projection_supported", "potential_high_trait_viable",
            "realised_high_trait_persists", "theoretical_k_above_four",
        ):
            row[key] = _bool(raw, key)
        key = _source_key(raw)
        patch_count = int(raw["patch_count"])
        if patch_count in by_source[key]:
            raise ValueError(f"duplicate source/patch row: {key}, n={patch_count}")
        by_source[key][patch_count] = row
        parsed.append(row)

    if len(by_source) != EXPECTED_CELL_COUNT * EXPECTED_SEEDS_PER_CELL * EXPECTED_REPLICATES_PER_SEED:
        raise ValueError(f"expected 1,200 source replicates, found {len(by_source)}")
    for key, patch_map in by_source.items():
        if tuple(sorted(patch_map)) != EXPECTED_PATCH_COUNTS:
            raise ValueError(f"patch-count coverage mismatch for {key}: {sorted(patch_map)}")
        source_flags = {bool(item["source_prepared"]) for item in patch_map.values()}
        if len(source_flags) != 1:
            raise ValueError(f"source-preparation flag varies across repeated measures for {key}")
    return parsed, dict(by_source)


def _paired_values(
    patch_map: dict[int, dict[str, object]], metric: str, patch_count: int
) -> tuple[float | None, float | None]:
    reference = patch_map[1]
    target = patch_map[patch_count]
    if not (
        reference["source_prepared"]
        and reference["projection_supported"]
        and target["source_prepared"]
        and target["projection_supported"]
    ):
        return None, None
    ref = reference[metric]
    value = target[metric]
    if ref is None or value is None or float(ref) <= 0.0:
        return None, None
    ratio = float(value) / float(ref)
    return ratio, 1.0 - ratio


def _cell_summary(
    parsed: list[dict[str, object]],
    by_source: dict[tuple[int, int, int], dict[int, dict[str, object]]],
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for cell_index in range(EXPECTED_CELL_COUNT):
        cell_sources = {key: value for key, value in by_source.items() if key[0] == cell_index}
        if len(cell_sources) != 100:
            raise ValueError(f"cell {cell_index}: expected 100 source replicates, found {len(cell_sources)}")
        for patch_count in EXPECTED_PATCH_COUNTS:
            rows = [patch_map[patch_count] for patch_map in cell_sources.values()]
            supported = [row for row in rows if row["source_prepared"] and row["projection_supported"]]
            first = rows[0]
            summary: dict[str, object] = {
                "primary_cell_index": cell_index,
                "mutation_rate": first["mutation_rate"],
                "area_reference": first["area_reference"],
                "interaction_feedback": first["interaction_feedback"],
                "patch_count": patch_count,
                "patch_area": first["patch_area"],
                "theoretical_k_fragment": first["theoretical_k_fragment"],
                "theoretical_k_above_four": first["theoretical_k_above_four"],
                "attempted_sources": len(rows),
                "source_prepared": sum(bool(row["source_prepared"]) for row in rows),
                "projection_supported": len(supported),
                "potential_high_trait_viable_probability": (
                    sum(bool(row["potential_high_trait_viable"]) for row in supported) / len(supported)
                    if supported else None
                ),
                "realised_high_trait_persists_probability": (
                    sum(bool(row["realised_high_trait_persists"]) for row in supported) / len(supported)
                    if supported else None
                ),
            }
            all_three_lower: list[bool] = []
            for source in cell_sources.values():
                ratios = [_paired_values(source, metric, patch_count)[0] for metric in PRIMARY_METRICS]
                if all(value is not None for value in ratios):
                    all_three_lower.append(all(float(value) < 1.0 for value in ratios))
            summary["all_three_below_one_patch_probability"] = (
                sum(all_three_lower) / len(all_three_lower) if all_three_lower else None
            )
            for metric in PRIMARY_METRICS:
                raw_median, raw_q1, raw_q3 = _median_iqr(
                    float(row[metric]) for row in supported if row[metric] is not None
                )
                ratios = []
                reductions = []
                for source in cell_sources.values():
                    ratio, reduction = _paired_values(source, metric, patch_count)
                    if ratio is not None:
                        ratios.append(ratio)
                        reductions.append(reduction)
                ratio_median, ratio_q1, ratio_q3 = _median_iqr(ratios)
                red_median, red_q1, red_q3 = _median_iqr(reductions)
                summary.update(
                    {
                        f"{metric}_median": raw_median,
                        f"{metric}_q1": raw_q1,
                        f"{metric}_q3": raw_q3,
                        f"{metric}_ratio_to_n1_median": ratio_median,
                        f"{metric}_ratio_to_n1_q1": ratio_q1,
                        f"{metric}_ratio_to_n1_q3": ratio_q3,
                        f"{metric}_reduction_from_n1_median": red_median,
                        f"{metric}_reduction_from_n1_q1": red_q1,
                        f"{metric}_reduction_from_n1_q3": red_q3,
                        f"{metric}_paired_n": len(ratios),
                    }
                )
            output.append(summary)
    return output


def _pooled_summary(
    cell_rows: list[dict[str, object]],
    by_source: dict[tuple[int, int, int], dict[int, dict[str, object]]],
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for patch_count in EXPECTED_PATCH_COUNTS:
        rows = [patch_map[patch_count] for patch_map in by_source.values()]
        supported = [row for row in rows if row["source_prepared"] and row["projection_supported"]]
        summary: dict[str, object] = {
            "patch_count": patch_count,
            "attempted_sources": len(rows),
            "source_prepared": sum(bool(row["source_prepared"]) for row in rows),
            "projection_supported": len(supported),
            "theoretical_k_fragment_min": min(float(row["theoretical_k_fragment"]) for row in rows),
            "theoretical_k_fragment_median": statistics.median(float(row["theoretical_k_fragment"]) for row in rows),
            "theoretical_k_fragment_max": max(float(row["theoretical_k_fragment"]) for row in rows),
        }
        all_three_lower = []
        for source in by_source.values():
            ratios = [_paired_values(source, metric, patch_count)[0] for metric in PRIMARY_METRICS]
            if all(value is not None for value in ratios):
                all_three_lower.append(all(float(value) < 1.0 for value in ratios))
        summary["all_three_below_one_patch_probability"] = (
            sum(all_three_lower) / len(all_three_lower) if all_three_lower else None
        )
        for metric in PRIMARY_METRICS:
            ratios = []
            reductions = []
            for source in by_source.values():
                ratio, reduction = _paired_values(source, metric, patch_count)
                if ratio is not None:
                    ratios.append(ratio)
                    reductions.append(reduction)
            ratio_median, ratio_q1, ratio_q3 = _median_iqr(ratios)
            red_median, red_q1, red_q3 = _median_iqr(reductions)
            summary.update(
                {
                    f"{metric}_ratio_to_n1_median": ratio_median,
                    f"{metric}_ratio_to_n1_q1": ratio_q1,
                    f"{metric}_ratio_to_n1_q3": ratio_q3,
                    f"{metric}_reduction_from_n1_median": red_median,
                    f"{metric}_reduction_from_n1_q1": red_q1,
                    f"{metric}_reduction_from_n1_q3": red_q3,
                    f"{metric}_paired_n": len(ratios),
                }
            )
        output.append(summary)
    return output


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _x_position(patch_count: int, left: float, width: float) -> float:
    lo = math.log2(EXPECTED_PATCH_COUNTS[0])
    hi = math.log2(EXPECTED_PATCH_COUNTS[-1])
    return left + width * (math.log2(patch_count) - lo) / (hi - lo)


def _gradient_svg(cell_rows: list[dict[str, object]], pooled_rows: list[dict[str, object]]) -> str:
    width, height = 1260, 480
    panels = (
        ("final_interaction_mean", 45),
        ("final_effective_size_mean", 455),
        ("realised_high_trait_mass_mean", 865),
    )
    plot_w, plot_h, top = 335, 300, 85
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        '<title id="title">Supplementary H3 fragmentation-gradient sensitivity</title>',
        '<desc id="desc">Three panels show each primary cell median retained fraction relative to its paired one-patch projection as isolated patch count increases at fixed total area. Thick black lines are pooled medians across all paired source replicates.</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="630" y="32" text-anchor="middle" font-family="sans-serif" font-size="20" font-weight="bold">Fixed-area fragmentation gradient from paired H1-prepared sources</text>',
    ]
    for panel_index, (metric, left) in enumerate(panels):
        parts.append(f'<text x="{left}" y="62" font-family="sans-serif" font-size="15" font-weight="bold">{chr(65+panel_index)}. {html.escape(METRIC_LABELS[metric])}</text>')
        parts.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="white" stroke="#444"/>')
        for y_value in (0.0, 0.25, 0.5, 0.75, 1.0, 1.25):
            y = top + plot_h - min(y_value, 1.25) / 1.25 * plot_h
            parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#e5e7eb"/>')
            parts.append(f'<text x="{left-8}" y="{y+4:.2f}" text-anchor="end" font-family="sans-serif" font-size="10">{y_value:.2g}</text>')
        one_y = top + plot_h - 1.0 / 1.25 * plot_h
        parts.append(f'<line x1="{left}" y1="{one_y:.2f}" x2="{left+plot_w}" y2="{one_y:.2f}" stroke="#666" stroke-dasharray="5,4"/>')
        for patch_count in EXPECTED_PATCH_COUNTS:
            x = _x_position(patch_count, left, plot_w)
            parts.append(f'<text x="{x:.2f}" y="{top+plot_h+22}" text-anchor="middle" font-family="sans-serif" font-size="10">{patch_count}</text>')
        parts.append(f'<text x="{left+plot_w/2}" y="{top+plot_h+45}" text-anchor="middle" font-family="sans-serif" font-size="11">number of isolated equal patches</text>')
        if panel_index == 0:
            parts.append(f'<text x="{left-34}" y="{top+plot_h/2}" text-anchor="middle" font-family="sans-serif" font-size="11" transform="rotate(-90 {left-34} {top+plot_h/2})">retained fraction relative to paired n=1</text>')

        for cell_index in range(EXPECTED_CELL_COUNT):
            rows = sorted(
                (row for row in cell_rows if int(row["primary_cell_index"]) == cell_index),
                key=lambda row: int(row["patch_count"]),
            )
            points = []
            for row in rows:
                value = row[f"{metric}_ratio_to_n1_median"]
                if value is None:
                    continue
                x = _x_position(int(row["patch_count"]), left, plot_w)
                clipped = max(0.0, min(float(value), 1.25))
                y = top + plot_h - clipped / 1.25 * plot_h
                points.append((x, y))
            if points:
                point_text = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
                parts.append(f'<polyline points="{point_text}" fill="none" stroke="#b8b8b8" stroke-width="1.2"/>')

        pooled_points = []
        for row in sorted(pooled_rows, key=lambda item: int(item["patch_count"])):
            value = row[f"{metric}_ratio_to_n1_median"]
            if value is None:
                continue
            x = _x_position(int(row["patch_count"]), left, plot_w)
            clipped = max(0.0, min(float(value), 1.25))
            y = top + plot_h - clipped / 1.25 * plot_h
            pooled_points.append((x, y))
        point_text = " ".join(f"{x:.2f},{y:.2f}" for x, y in pooled_points)
        parts.append(f'<polyline points="{point_text}" fill="none" stroke="#111" stroke-width="3"/>')
        for x, y in pooled_points:
            parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.5" fill="#111"/>')
    parts.append('<text x="630" y="466" text-anchor="middle" font-family="sans-serif" font-size="10">Grey lines: cell-specific medians (12 frozen primary cells). Black: pooled paired median. Values above 1 indicate an increase relative to the paired one-patch outcome.</text>')
    parts.append('</svg>')
    return "\n".join(parts) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    raw_rows = _read_rows(args.input_root)
    parsed, by_source = _validated_rows(raw_rows)
    cell_rows = _cell_summary(parsed, by_source)
    pooled_rows = _pooled_summary(cell_rows, by_source)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(args.output_dir / "h3_fragmentation_gradient_records.csv", [dict(row) for row in parsed])
    _write_csv(args.output_dir / "h3_fragmentation_gradient_cell_summary.csv", cell_rows)
    _write_csv(args.output_dir / "h3_fragmentation_gradient_pooled_summary.csv", pooled_rows)
    (args.output_dir / "figure_s_fragmentation_gradient.svg").write_text(
        _gradient_svg(cell_rows, pooled_rows), encoding="utf-8"
    )
    metadata = {
        "analysis": "post-review H3 fragmentation-gradient sensitivity",
        "attempted_source_replicates": len(by_source),
        "patch_counts": list(EXPECTED_PATCH_COUNTS),
        "record_count": len(parsed),
        "cell_count": EXPECTED_CELL_COUNT,
        "fresh_seed_family": [20260820, 20260821, 20260822, 20260823, 20260824],
        "evidence_label": "new supplementary finite Type S sensitivity evidence",
        "warning_endpoints_evaluated": False,
    }
    (args.output_dir / "h3_fragmentation_gradient_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"Aggregated {len(parsed)} rows / {len(by_source)} source replicates; "
        f"wrote {len(cell_rows)} cell summaries and {len(pooled_rows)} pooled rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
