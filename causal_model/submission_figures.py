from __future__ import annotations

import csv
import html
import math
from pathlib import Path


def _text(x: float, y: float, value: str, *, size: int = 16, anchor: str = "middle", weight: str = "normal") -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
        f'font-family="sans-serif" font-size="{size}" font-weight="{weight}">{html.escape(value)}</text>'
    )


def write_architecture_figure(path: Path) -> None:
    width, height = 1200, 650
    labels = [
        (120, "patch size"),
        (345, "interaction state q"),
        (585, "potential high-trait viability"),
        (835, "realised occupancy / local Ne"),
        (1080, "allele persistence / diversity"),
    ]
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        '<title id="title">Theorem-to-closure architecture for eco-genetic state separation</title>',
        '<desc id="desc">A causal research architecture links patch size, interaction state, potential trait viability, realised occupancy and effective size, and genetic summaries. Arrow labels distinguish exact, conditional and finite evidence.</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<defs><marker id="arrow" markerWidth="10" markerHeight="10" refX="8" refY="5" orient="auto"><path d="M0,0 L0,10 L9,5 z" fill="#111"/></marker></defs>',
        _text(600, 48, "Fragmentation can separate biological states", size=28, weight="bold"),
        _text(600, 82, "The arrows do not share one evidential status", size=16),
    ]
    y = 235
    for x, label in labels:
        w = 185 if x not in (585, 835, 1080) else 210
        lines.append(f'<rect x="{x-w/2:.1f}" y="{y-55}" width="{w}" height="110" rx="10" fill="white" stroke="#111" stroke-width="2"/>')
        parts = label.split(" / ")
        if len(parts) == 1 and len(label) <= 24:
            lines.append(_text(x, y+6, label, size=17, weight="bold"))
        else:
            if "potential" in label:
                lines.append(_text(x, y-9, "potential high-trait", size=16, weight="bold"))
                lines.append(_text(x, y+17, "viability", size=16, weight="bold"))
            elif "realised" in label:
                lines.append(_text(x, y-9, "realised occupancy", size=16, weight="bold"))
                lines.append(_text(x, y+17, "/ local effective size", size=15))
            elif "allele" in label:
                lines.append(_text(x, y-9, "allele persistence", size=16, weight="bold"))
                lines.append(_text(x, y+17, "/ genetic summaries", size=15))
    arrow_pairs = [(212, 252), (438, 475), (690, 725), (940, 970)]
    evidence = ["declared patchwise closure", "Type C trait margin", "finite Type S closure", "finite genetic summaries"]
    for (x1, x2), label in zip(arrow_pairs, evidence):
        lines.append(f'<line x1="{x1}" y1="{y}" x2="{x2}" y2="{y}" stroke="#111" stroke-width="2" marker-end="url(#arrow)"/>')
        lines.append(_text((x1+x2)/2, y-73, label, size=12))
    lines.extend([
        '<rect x="90" y="390" width="1020" height="150" rx="10" fill="white" stroke="#777" stroke-width="1.5"/>',
        _text(600, 425, "Evidence discipline", size=20, weight="bold"),
        _text(600, 462, "Type T: exact statements for specified maps", size=15),
        _text(600, 491, "Type C: implications after an ecological / life-cycle closure is declared", size=15),
        _text(600, 520, "Type S: finite stochastic evidence — not a theorem or empirical parameter estimate", size=15),
        _text(600, 603, "No single deterioration label substitutes for all five state objects", size=18, weight="bold"),
        '</svg>',
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_analytical_boundaries_figure(path: Path) -> None:
    # Display-only canonical example already used in the repository: kappa=8, A/A_ref=1, theta=0.5.
    # For K=8, theta=0.5 is inside the exact strict-bistability interval.
    kappa = 8.0
    K = 8.0
    theta = 0.5
    q_lo = 0.021247987961365615
    q_mid = 0.5
    q_hi = 0.9787520120386344
    q_minus = (1.0 - math.sqrt(1.0 - 4.0 / K)) / 2.0
    q_plus = (1.0 + math.sqrt(1.0 - 4.0 / K)) / 2.0
    theta_minus = (K*q_minus - math.log(q_minus/(1-q_minus))) / kappa
    theta_plus = (K*q_plus - math.log(q_plus/(1-q_plus))) / kappa

    width, height = 1200, 700
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        '<title id="title">Exact analytical boundaries for interaction feedback and allele-frequency mixing</title>',
        '<desc id="desc">Panel A shows the canonical sigmoid map and identity for K equals eight and theta equals one half, with two stable and one unstable fixed points. Panel B shows the row-stochastic allele-frequency mixing operator and its common-floor and focal-rescue bounds.</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        _text(600, 45, "Exact analytical statements are operator-specific", size=28, weight="bold"),
        _text(300, 88, "A  Canonical interaction map", size=20, weight="bold"),
        _text(900, 88, "B  Network allele-frequency mixing", size=20, weight="bold"),
    ]
    # Panel A axes and curves.
    left, right, top, bottom = 80, 555, 135, 555
    lines.extend([
        f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#111" stroke-width="2"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#111" stroke-width="2"/>',
        _text((left+right)/2, 603, "q(t)", size=15),
        _text(40, 350, "q(t+1)", size=15),
    ])
    pts_map = []
    pts_identity = []
    for i in range(201):
        q = i / 200.0
        z = kappa * (q - theta)
        fq = 1.0 / (1.0 + math.exp(-z))
        x = left + q * (right-left)
        y = bottom - fq * (bottom-top)
        yi = bottom - q * (bottom-top)
        pts_map.append((x,y))
        pts_identity.append((x,yi))
    lines.append('<polyline points="' + ' '.join(f"{x:.1f},{y:.1f}" for x,y in pts_map) + '" fill="none" stroke="#111" stroke-width="2.5"/>')
    lines.append('<polyline points="' + ' '.join(f"{x:.1f},{y:.1f}" for x,y in pts_identity) + '" fill="none" stroke="#777" stroke-width="1.8" stroke-dasharray="7 5"/>')
    for q, label, stable in ((q_lo,"low",True),(q_mid,"middle",False),(q_hi,"high",True)):
        x = left + q*(right-left)
        y = bottom - q*(bottom-top)
        if stable:
            lines.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="#111"/>')
        else:
            lines.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="white" stroke="#111" stroke-width="2"/>')
        lines.append(_text(x, y-15, label, size=12))
    lines.extend([
        _text(315, 635, f"K=8; strict bistability interval θ=({theta_minus:.3f}, {theta_plus:.3f}); display θ=0.500", size=13),
        _text(315, 660, "Exact for the stated one-state sigmoid map", size=13, weight="bold"),
    ])

    # Panel B schematic.
    lines.extend([
        '<defs><marker id="arrow2" markerWidth="10" markerHeight="10" refX="8" refY="5" orient="auto"><path d="M0,0 L0,10 L9,5 z" fill="#111"/></marker></defs>',
        '<circle cx="735" cy="200" r="42" fill="white" stroke="#111" stroke-width="2"/>',
        '<circle cx="735" cy="350" r="42" fill="white" stroke="#111" stroke-width="2"/>',
        '<circle cx="735" cy="500" r="42" fill="white" stroke="#111" stroke-width="2"/>',
        '<rect x="990" y="305" width="130" height="90" rx="10" fill="white" stroke="#111" stroke-width="2"/>',
        _text(735, 205, "source 1", size=14, weight="bold"),
        _text(735, 355, "source 2", size=14, weight="bold"),
        _text(735, 505, "source j", size=14, weight="bold"),
        _text(1055, 345, "destination i", size=14, weight="bold"),
        _text(1055, 370, "p'i", size=15),
        '<line x1="777" y1="200" x2="990" y2="325" stroke="#111" stroke-width="1.8" marker-end="url(#arrow2)"/>',
        '<line x1="777" y1="350" x2="990" y2="350" stroke="#111" stroke-width="1.8" marker-end="url(#arrow2)"/>',
        '<line x1="777" y1="500" x2="990" y2="375" stroke="#111" stroke-width="1.8" marker-end="url(#arrow2)"/>',
        _text(885, 225, "M_i1", size=12),
        _text(885, 338, "M_i2", size=12),
        _text(885, 480, "M_ij", size=12),
        _text(900, 545, "p'i = Σj Mij pj;  Σj Mij = 1", size=17, weight="bold"),
        _text(900, 582, "common floor: pj ≥ pmin for all j  ⇒  p'i ≥ pmin", size=15),
        _text(900, 615, "focal bound: p'i ≥ Σj Mij bj", size=15),
        _text(900, 655, "Composition mixing ≠ demographic or trait rescue", size=14, weight="bold"),
        '</svg>',
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_viability_occupancy_figure(path: Path) -> None:
    width, height = 1100, 650
    patch_counts = [1,2,3,4,6,8,12,16]
    left, right, top, bottom = 100, 1030, 130, 520
    x_positions = [left + i*(right-left)/(len(patch_counts)-1) for i in range(len(patch_counts))]
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        '<title id="title">Potential high-trait viability and realised occupancy separate after fragmentation</title>',
        '<desc id="desc">Potential high-trait viability is present in all 1037 one-patch supported outcomes and absent in all 1037 supported outcomes at every tested subdivision. Realised high-trait occupancy remains approximately 99.6 to 100 percent at the generation-30 endpoint across patch counts.</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        _text(550, 48, "Potential viability and realised occupancy are distinct states", size=27, weight="bold"),
        _text(550, 78, "Fresh preregistered fragmentation gradient; 1,037 common supported sources", size=15),
        f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#111" stroke-width="2"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#111" stroke-width="2"/>',
    ]
    for pct in (0,25,50,75,100):
        y = bottom - pct/100*(bottom-top)
        lines.append(f'<line x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}" stroke="#ddd" stroke-width="1"/>')
        lines.append(_text(left-12, y+5, f"{pct}%", size=12, anchor="end"))
    # realised occupancy band: only supported range is encoded, not fabricated patch-specific values.
    y100 = top
    y996 = bottom - 0.996*(bottom-top)
    band_h = max(2.0, y996-y100)
    lines.append(f'<rect x="{left}" y="{y100:.1f}" width="{right-left}" height="{band_h:.1f}" fill="#d9d9d9" stroke="#777" stroke-width="1"/>')
    lines.append(_text(720, 118, "realised high-trait occupancy at generation 30: ~99.6–100%", size=14))
    # potential viability points.
    viable = [100,0,0,0,0,0,0,0]
    pts=[]
    for x,n,pct in zip(x_positions,patch_counts,viable):
        y = bottom - pct/100*(bottom-top)
        pts.append((x,y))
        lines.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="white" stroke="#111" stroke-width="2"/>')
        lines.append(_text(x, bottom+30, str(n), size=12))
    lines.append('<polyline points="'+' '.join(f"{x:.1f},{y:.1f}" for x,y in pts)+'" fill="none" stroke="#111" stroke-width="2.5"/>')
    lines.extend([
        _text(560, 575, "number of equal isolated patches at fixed total area", size=14),
        _text(58, 340, "supported outcomes", size=14),
        _text(240, 180, "potential viability: 1037/1037", size=13, weight="bold"),
        _text(625, 495, "potential viability: 0/1037 at every tested subdivision", size=13, weight="bold"),
        _text(550, 620, "Finite Type S state separation — not a universal lag law", size=15, weight="bold"),
        '</svg>',
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def validate_gradient_summary(pooled_csv: Path) -> None:
    rows = list(csv.DictReader(pooled_csv.open(encoding="utf-8", newline="")))
    if [int(r["patch_count"]) for r in rows] != [1,2,3,4,6,8,12,16]:
        raise RuntimeError("fragmentation-gradient patch-count contract drifted")
    if any(int(r["projection_supported"]) != 1037 for r in rows):
        raise RuntimeError("common projection-supported denominator drifted")
    expected = {
        2: (0.0017444307534129745, 0.22131147540983603, 0.28291808988011524),
        4: (0.001446627790951059, 0.11270491803278689, 0.30180017351069993),
        16: (0.0012440085611615528, 0.033057851239669415, 0.39388020833333326),
    }
    by_n = {int(r["patch_count"]): r for r in rows}
    for n, values in expected.items():
        got = (
            float(by_n[n]["final_interaction_mean_ratio_to_n1_median"]),
            float(by_n[n]["final_effective_size_mean_ratio_to_n1_median"]),
            float(by_n[n]["realised_high_trait_mass_mean_ratio_to_n1_median"]),
        )
        if got != values:
            raise RuntimeError(f"fragmentation-gradient headline drifted at n={n}: {got}")
