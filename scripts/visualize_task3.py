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
        feature_count = 0
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"Malformed row in {path}")
            feature_count += 1
            for name, value in zip(header, row):
                if name in sums:
                    number = float(value)
                    if not math.isfinite(number) or number < 0:
                        raise ValueError(f"Invalid count in {path}: {value}")
                    sums[name] += number
                    detected[name] += number > 0
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
                        "total_counts": sums[column], "detected_features": detected[column]})
    return metrics


def svg_chart(technology: str, rows: list[dict]) -> str:
    width, height = 1100, 690
    left, right = 105, 1030
    plot_width = right - left
    top1, panel_h = 125, 205
    top2 = 405
    gap = plot_width / len(rows)
    bar_w = min(48, gap * 0.62)
    max_total = max(row["total_counts"] for row in rows)
    max_detected = max(row["detected_features"] for row in rows)
    max_total = max_total or 1
    max_detected = max_detected or 1
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="#ffffff"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#243342}.title{font-size:25px;font-weight:700}.subtitle{font-size:14px;fill:#526577}.panel{font-size:18px;font-weight:700}.axis{font-size:12px;fill:#526577}.sample{font-size:12px}.value{font-size:11px;fill:#34495e}.grid{stroke:#dce3e8;stroke-width:1}</style>',
             f'<text x="{left}" y="43" class="title">{html.escape(technology.replace("_", " ").title())} raw-count sample summaries</text>',
             f'<text x="{left}" y="70" class="subtitle">Descriptive totals across all {rows[0]["features"]:,} rows; samples are shown individually</text>']
    for label, top, field, maximum, color in [
            ("Total raw counts per sample", top1, "total_counts", max_total, "#6688a8"),
            ("Detected rows (count > 0)", top2, "detected_features", max_detected, "#56a28c")]:
        baseline = top + panel_h
        parts.append(f'<text x="{left}" y="{top-20}" class="panel">{label}</text>')
        for tick in range(5):
            fraction = tick / 4
            y = baseline - fraction * panel_h
            value = maximum * fraction
            shown = f"{value / 1_000_000:.1f}M" if maximum >= 1_000_000 else f"{value / 1000:.0f}k"
            parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}" class="grid"/>')
            parts.append(f'<text x="{left-12}" y="{y+4:.1f}" text-anchor="end" class="axis">{shown}</text>')
        for index, row in enumerate(rows):
            xcenter = left + gap * (index + 0.5)
            value = row[field]
            bar_h = value / maximum * panel_h
            y = baseline - bar_h
            color = COLORS[row["condition"]] if field == "total_counts" else color
            parts.append(f'<rect x="{xcenter-bar_w/2:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}" rx="3" fill="{color}"/>')
            value_text = f"{value:,.0f}"
            parts.append(f'<text x="{xcenter:.1f}" y="{max(top+10,y-6):.1f}" text-anchor="middle" class="value">{value_text}</text>')
            condition = "Ex vivo" if row["condition"] == "ex_vivo" else "LPS"
            parts.append(f'<text x="{xcenter:.1f}" y="{baseline+19}" text-anchor="middle" class="sample">HC{row["donor"]} {condition}</text>')
    parts.extend(['<rect x="105" y="664" width="12" height="12" fill="#3977a8"/><text x="123" y="675" class="axis">Ex vivo</text>',
                  '<rect x="195" y="664" width="12" height="12" fill="#d4773c"/><text x="213" y="675" class="axis">3 h LPS</text>',
                  '<text x="1030" y="675" text-anchor="end" class="axis">Sample-to-metadata mapping uses volunteer number and condition in the supplied sample metadata.</text>',
                  '</svg>'])
    return "\n".join(parts) + "\n"


def main() -> None:
    metadata = read_metadata()
    all_metrics = []
    for technology, path in TABLES.items():
        metrics = summarize(path, technology, metadata)
        all_metrics.extend(metrics)
        FIGURES.mkdir(parents=True, exist_ok=True)
        (FIGURES / f"{technology}_sample_summary.svg").write_text(svg_chart(technology, metrics), encoding="utf-8")
    with SUMMARY.open("w", encoding="utf-8", newline="") as handle:
        fields = ["technology", "sample_column", "geo_accession", "donor", "condition", "features", "total_counts", "detected_features"]
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_metrics)
    print(f"Wrote two SVG figures under {FIGURES.relative_to(ROOT)} and {SUMMARY.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
