#!/usr/bin/env python3
"""Create descriptive SVG summaries of the checked-in GSE297994 raw counts.

Uses only Python's standard library. Run from any directory with:
    python scripts/visualize_task3.py
"""
from __future__ import annotations

import csv
import html
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "metadata" / "sample_metadata.tsv"
FIGURES = ROOT / "results" / "figures"
SUMMARY = ROOT / "results" / "task3_sample_metrics.tsv"
TABLES = {
    "short_read": ROOT / "data" / "cleaned" / "short_read_raw_counts.txt",
    "long_read": ROOT / "data" / "cleaned" / "long_read_raw_counts.txt",
}
COLORS = {"ex_vivo": "#3977a8", "LPS": "#d4773c"}


def read_metadata() -> dict[tuple[str, int, str], str]:
    result: dict[tuple[str, int, str], str] = {}
    with META.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            match = re.search(r"healthy volunteer (\d+)", row["title"], re.I)
            if not match:
                raise ValueError(f"Cannot identify volunteer in metadata title: {row['title']}")
            key = (row["technology"], int(match.group(1)), row["condition"])
            if key in result:
                raise ValueError(f"Duplicate sample metadata key: {key}")
            result[key] = row["geo_accession"]
    return result


def summarize(path: Path, technology: str, metadata: dict) -> list[dict]:
    metrics = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        id_fields = {"gene_id", "ids"}
        sample_columns = [name for name in header if name not in id_fields]
        sums = {name: 0 for name in sample_columns}
        detected = {name: 0 for name in sample_columns}
        values = {name: [] for name in sample_columns}
        feature_count = 0
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"Malformed row in {path}")
            feature_count += 1
            for name, value in zip(header, row):
                if name in sums:
                    number = float(value)
                    if not math.isfinite(number) or number < 0 or not number.is_integer():
                        raise ValueError(f"Invalid count in {path}: {value}")
                    sums[name] += int(number)
                    detected[name] += number > 0
                    values[name].append(int(number))
    for column in sample_columns:
        match = re.search(r"HC(\d+)(exvivo|LPS)", column, re.I)
        if not match:
            raise ValueError(f"Cannot map count column to metadata: {column}")
        donor = int(match.group(1))
        condition = "ex_vivo" if match.group(2).lower() == "exvivo" else "LPS"
        accession = metadata.get((technology.replace("_", "-"), donor, condition))
        if accession is None:
            raise ValueError(f"No metadata match for {technology}, donor {donor}, {condition}")
        metrics.append({"technology": technology, "sample_column": column,
                        "geo_accession": accession, "donor": donor,
                        "condition": condition, "features": feature_count,
                        "total_counts": sums[column], "detected_features": detected[column],
                        "top100_share_percent": (
                            100 * sum(sorted(values[column], reverse=True)[:100]) / sums[column]
                            if technology == "short_read" and sums[column] else "NA")})
    return metrics


def nice_axis(maximum: float, intervals: int = 5) -> tuple[float, list[float]]:
    raw_step = maximum / intervals if maximum else 1
    magnitude = 10 ** math.floor(math.log10(raw_step))
    step = next(factor * magnitude for factor in (1, 2, 2.5, 5, 10)
                if factor * magnitude >= raw_step)
    axis_max = math.ceil(maximum / step) * step
    ticks = [i * step for i in range(int(round(axis_max / step)) + 1)]
    return axis_max, ticks


def tick_text(value: float) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:g}M"
    if value >= 1_000:
        return f"{value / 1_000:g}k"
    return f"{value:g}"


def svg_chart(technology: str, rows: list[dict]) -> str:
    donors = sorted({row["donor"] for row in rows})
    rows = sorted(rows, key=lambda row: (row["condition"] != "ex_vivo", row["donor"]))
    width, height = 1240, 770
    left, right = 100, 1170
    plot_width = right - left
    top1, panel_h = 145, 215
    top2 = 465
    panels = [("Total raw counts per sample", "total_counts", top1),
              ("Detected rows (count > 0)", "detected_features", top2)]
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="#fff"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#243342}.title{font-size:25px;font-weight:700}.subtitle{font-size:14px;fill:#526577}.panel{font-size:18px;font-weight:700}.axis{font-size:12px;fill:#526577}.sample{font-size:14px}.value{font-size:11px;fill:#34495e}.grid{stroke:#dce3e8;stroke-width:1}.group{font-size:13px;font-weight:700}</style>',
             f'<text x="{left}" y="42" class="title">{html.escape(technology.replace("_", " ").title())} raw-count sample summaries</text>',
             f'<text x="{left}" y="69" class="subtitle">Descriptive sample-level totals and detected rows; conditions are grouped and colored consistently</text>']
    gap = plot_width / len(rows)
    bar_w = min(44, gap * 0.64)
    group_mid = (left + gap * len(donors) / 2, left + gap * (len(donors) + len(rows)) / 2)
    for label, field, top in panels:
        maximum = max(row[field] for row in rows)
        axis_max, ticks = nice_axis(maximum * 1.14)
        baseline = top + panel_h
        parts.append(f'<text x="{left}" y="{top-18}" class="panel">{label}</text>')
        for value in ticks:
            y = baseline - value / axis_max * panel_h
            parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}" class="grid"/>')
            parts.append(f'<text x="{left-12}" y="{y+4:.1f}" text-anchor="end" class="axis">{tick_text(value)}</text>')
        for index, row in enumerate(rows):
            xcenter = left + gap * (index + 0.5)
            value = row[field]
            bar_h = value / axis_max * panel_h
            y = baseline - bar_h
            color = COLORS[row["condition"]]
            parts.append(f'<rect x="{xcenter-bar_w/2:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}" rx="3" fill="{color}"/>')
            if field == "total_counts":
                parts.append(f'<text x="{xcenter:.1f}" y="{y-6:.1f}" text-anchor="middle" class="value">{value/1_000_000:.1f}M</text>')
            parts.append(f'<text x="{xcenter:.1f}" y="{baseline+20}" text-anchor="middle" class="sample">HC{row["donor"]}</text>')
        caption_y = baseline + 47
        for condition, mid in zip(("ex_vivo", "LPS"), group_mid):
            label = "Ex vivo" if condition == "ex_vivo" else "3 h LPS"
            parts.append(f'<rect x="{mid-43:.1f}" y="{caption_y-13}" width="12" height="12" rx="2" fill="{COLORS[condition]}"/>')
            parts.append(f'<text x="{mid-25:.1f}" y="{caption_y-3}" class="group">{label}</text>')
        if len(donors) < len(rows):
            divider_x = left + gap * len(donors)
            parts.append(f'<line x1="{divider_x:.1f}" y1="{top}" x2="{divider_x:.1f}" y2="{baseline+24}" stroke="#9aa8b4" stroke-dasharray="4 4"/>')
    parts.append('</svg>')
    return "\n".join(parts) + "\n"


def svg_top100_chart(rows: list[dict]) -> str:
    rows = [row for row in rows if row["technology"] == "short_read"]
    by_donor = {row["donor"]: {} for row in rows}
    for row in rows:
        by_donor[row["donor"]][row["condition"]] = row
    width, height = 1120, 680
    left, right, top, bottom = 110, 1050, 130, 535
    axis_max, ticks = nice_axis(max(float(r["top100_share_percent"]) for r in rows) * 1.12, intervals=6)
    xgap = (right - left) / 8
    offset = 15
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="#fff"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#243342}.title{font-size:25px;font-weight:700}.subtitle{font-size:14px;fill:#526577}.axis{font-size:13px;fill:#526577}.donor{font-size:14px}.value{font-size:12px;fill:#34495e}.grid{stroke:#dce3e8;stroke-width:1}</style>',
             f'<text x="{left}" y="45" class="title">Top 100 genes account for a larger share of LPS sample counts</text>',
             f'<text x="{left}" y="73" class="subtitle">Within each sample, sum of its 100 highest-count gene rows divided by all gene counts</text>',
             f'<text x="{left}" y="96" class="subtitle">Each line joins the same volunteer across conditions (n = 8 paired volunteers)</text>',
             '<rect x="730" y="104" width="13" height="13" rx="2" fill="#3977a8"/><text x="750" y="115" class="axis">Ex vivo</text>',
             '<rect x="830" y="104" width="13" height="13" rx="2" fill="#d4773c"/><text x="850" y="115" class="axis">3 h LPS</text>']
    for value in ticks:
        y = bottom - value / axis_max * (bottom - top)
        parts.extend([f'<line x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}" class="grid"/>',
                      f'<text x="{left-12}" y="{y+4:.1f}" text-anchor="end" class="axis">{value:g}%</text>'])
    for donor in sorted(by_donor):
        pair = by_donor[donor]
        exv = float(pair["ex_vivo"]["top100_share_percent"])
        lps = float(pair["LPS"]["top100_share_percent"])
        x = left + xgap * (donor - 0.5)
        y1, y2 = bottom - exv / axis_max * (bottom-top), bottom - lps / axis_max * (bottom-top)
        parts.append(f'<line x1="{x-offset}" y1="{y1:.1f}" x2="{x+offset}" y2="{y2:.1f}" stroke="#8493a0" stroke-width="2"/>')
        for xpos, y, value, condition in [(x-offset, y1, exv, "ex_vivo"), (x+offset, y2, lps, "LPS")]:
            color = COLORS[condition]
            parts.append(f'<circle cx="{xpos}" cy="{y:.1f}" r="7" fill="{color}" stroke="#fff" stroke-width="2"/>')
            if condition == "ex_vivo":
                label_x, anchor = xpos - 12, "end"
            else:
                label_x, anchor = xpos + 12, "start"
            parts.append(f'<text x="{label_x}" y="{y+4:.1f}" text-anchor="{anchor}" class="value">{value:.1f}%</text>')
        parts.append(f'<text x="{x}" y="{bottom+29}" text-anchor="middle" class="donor">HC{donor}</text>')
    parts.append('</svg>')
    return "\n".join(parts) + "\n"


def main() -> None:
    metadata = read_metadata()
    all_metrics = []
    for technology, path in TABLES.items():
        metrics = summarize(path, technology, metadata)
        all_metrics.extend(metrics)
        FIGURES.mkdir(parents=True, exist_ok=True)
        (FIGURES / f"{technology}_sample_summary.svg").write_text(svg_chart(technology, metrics), encoding="utf-8")
    (FIGURES / "short_read_top100_gene_share.svg").write_text(svg_top100_chart(all_metrics), encoding="utf-8")
    with SUMMARY.open("w", encoding="utf-8", newline="") as handle:
        fields = ["technology", "sample_column", "geo_accession", "donor", "condition", "features", "total_counts", "detected_features", "top100_share_percent"]
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_metrics)
    print(f"Wrote three SVG figures under {FIGURES.relative_to(ROOT)} and {SUMMARY.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
