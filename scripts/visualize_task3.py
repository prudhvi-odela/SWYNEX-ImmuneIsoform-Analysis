#!/usr/bin/env python3
"""Create reproducible descriptive figures from the cleaned GSE297994 counts.

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
SVG_NS = 'xmlns="http://www.w3.org/2000/svg"'


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


def read_counts(path: Path) -> tuple[list[str], list[dict]]:
    samples, features = [], []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_indexes = [(i, name) for i, name in enumerate(header) if name not in {"gene_id", "ids"}]
        samples = [name for _, name in sample_indexes]
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"{path.name}:{line_number}: malformed row")
            values = []
            for index, name in sample_indexes:
                value = float(row[index])
                if not math.isfinite(value) or value < 0 or not value.is_integer():
                    raise ValueError(f"Invalid raw count in {path.name}:{line_number}:{name}: {row[index]}")
                values.append(int(value))
            features.append({"id": row[0], "values": values})
    if not samples or not features:
        raise ValueError(f"No sample columns or rows in {path}")
    return samples, features


def sample_info(samples: list[str], technology: str, metadata: dict) -> list[dict]:
    output = []
    for name in samples:
        match = re.search(r"HC(\d+)(exvivo|LPS)", name, re.I)
        if not match:
            raise ValueError(f"Cannot map count column to metadata: {name}")
        donor = int(match.group(1))
        condition = "ex_vivo" if match.group(2).lower() == "exvivo" else "LPS"
        accession = metadata.get((technology.replace("_", "-"), donor, condition))
        if accession is None:
            raise ValueError(f"No metadata match for {technology}, donor {donor}, {condition}")
        output.append({"technology": technology, "sample_column": name, "geo_accession": accession,
                       "donor": donor, "condition": condition})
    return output


def median(values: list[float]) -> float:
    if not values:
        raise ValueError("Cannot calculate a median of no values")
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    low = int(position)
    high = min(low + 1, len(ordered) - 1)
    return ordered[low] + (position - low) * (ordered[high] - ordered[low])


def nice_axis(maximum: float, intervals: int = 5, minimum: float = 0) -> tuple[float, list[float]]:
    span = maximum - minimum
    raw_step = span / intervals if span else 1
    magnitude = 10 ** math.floor(math.log10(abs(raw_step)))
    step = next(factor * magnitude for factor in (1, 2, 2.5, 5, 10)
                if factor * magnitude >= raw_step)
    axis_min = math.floor(minimum / step) * step
    axis_max = math.ceil(maximum / step) * step
    count = int(round((axis_max - axis_min) / step))
    return axis_max, [axis_min + i * step for i in range(count + 1)]


def format_tick(value: float) -> str:
    if abs(value) >= 1_000_000:
        return f"{value / 1_000_000:g}M"
    if abs(value) >= 1_000:
        return f"{value / 1_000:g}k"
    return f"{value:g}"


def analyze(technology: str, samples: list[str], features: list[dict], metadata: dict) -> dict:
    info = sample_info(samples, technology, metadata)
    group_sizes = {condition: sum(s["condition"] == condition for s in info)
                   for condition in ("ex_vivo", "LPS")}
    minimum_group = min(group_sizes.values())
    library_sizes = [sum(feature["values"][j] for feature in features) for j in range(len(samples))]
    if any(size <= 0 for size in library_sizes):
        raise ValueError(f"Empty library in {technology}")
    cpm = [[value / library_sizes[j] * 1_000_000 for j, value in enumerate(feature["values"])]
           for feature in features]
    passing = [i for i, row in enumerate(cpm)
               if sum(value >= 1 for value in row) >= minimum_group]
    if len(passing) < 500:
        raise ValueError(f"Only {len(passing)} features pass the CPM filter for {technology}")

    # DESeq2's default median-of-ratios ('ratio') estimator uses features
    # with positive counts in every sample when calculating geometric means.
    common_positive = [feature["values"] for feature in features
                       if all(value > 0 for value in feature["values"])]
    if not common_positive:
        raise ValueError(f"No all-positive features for median-of-ratios in {technology}")
    ratio_lists = [[] for _ in samples]
    for values in common_positive:
        geometric_mean = math.exp(sum(math.log(value) for value in values) / len(values))
        for j, value in enumerate(values):
            ratio_lists[j].append(value / geometric_mean)
    size_factors = [median(ratios) for ratios in ratio_lists]
    factor_geomean = math.exp(sum(math.log(value) for value in size_factors) / len(size_factors))
    size_factors = [value / factor_geomean for value in size_factors]

    log_normalized = [[math.log2(feature["values"][j] / size_factors[j] + 1)
                       for j in range(len(samples))] for feature in features]
    cpm_medians, normalized_medians = [], []
    for j in range(len(samples)):
        expressed = [i for i in passing if features[i]["values"][j] > 0]
        cpm_medians.append(median([cpm[i][j] for i in expressed]))
        normalized_medians.append(median([log_normalized[i][j] for i in expressed]))

    metrics = []
    for j, item in enumerate(info):
        top_share = "NA"
        if technology == "short_read":
            top = sorted((feature["values"][j] for feature in features), reverse=True)[:100]
            top_share = 100 * sum(top) / library_sizes[j]
        metrics.append({**item, "features": len(features), "library_total": library_sizes[j],
                        "cpm_filter_features": len(passing), "median_expressed_cpm": cpm_medians[j],
                        "size_factor": size_factors[j], "median_log2_normalized_expressed": normalized_medians[j],
                        "top100_share_percent": top_share})
    return {"technology": technology, "samples": info, "features": features, "cpm": cpm,
            "passing": passing, "log_normalized": log_normalized,
            "size_factors": size_factors, "metrics": metrics}


def svg_header(width: int, height: int, title: str, subtitle: str = "") -> list[str]:
    parts = [f'<svg {SVG_NS} width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="#fff"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#243342}.title{font-size:25px;font-weight:700}.subtitle{font-size:14px;fill:#526577}.panel{font-size:18px;font-weight:700}.axis{font-size:12px;fill:#526577}.sample{font-size:12px}.value{font-size:11px;fill:#34495e}.grid{stroke:#dce3e8;stroke-width:1}.group{font-size:13px;font-weight:700}</style>',
             f'<text x="70" y="42" class="title">{html.escape(title)}</text>']
    if subtitle:
        parts.append(f'<text x="70" y="68" class="subtitle">{html.escape(subtitle)}</text>')
    return parts


def svg_top100_chart(result: dict) -> str:
    rows = [m for m in result["metrics"] if m["technology"] == "short_read"]
    by_donor = {m["donor"]: {} for m in rows}
    for row in rows:
        by_donor[row["donor"]][row["condition"]] = row
    width, height = 1120, 680
    left, right, top, bottom = 110, 1050, 130, 535
    axis_max, ticks = nice_axis(max(m["top100_share_percent"] for m in rows) * 1.12, intervals=6)
    xgap, offset = (right - left) / 8, 15
    parts = svg_header(width, height, "Top 100 genes account for a larger share of LPS sample counts",
                       "Within each sample, the 100 highest-count gene rows are divided by all gene counts")
    parts.extend([f'<text x="{left}" y="96" class="subtitle">Each line pairs one volunteer across conditions (n = 8 pairs)</text>',
                  '<rect x="730" y="104" width="13" height="13" rx="2" fill="#3977a8"/><text x="750" y="115" class="axis">Ex vivo</text>',
                  '<rect x="830" y="104" width="13" height="13" rx="2" fill="#d4773c"/><text x="850" y="115" class="axis">3 h LPS</text>'])
    for value in ticks:
        y = bottom - value / axis_max * (bottom - top)
        parts.extend([f'<line x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}" class="grid"/>',
                      f'<text x="{left-12}" y="{y+4:.1f}" text-anchor="end" class="axis">{value:g}%</text>'])
    for donor, pair in sorted(by_donor.items()):
        exv, lps = pair["ex_vivo"]["top100_share_percent"], pair["LPS"]["top100_share_percent"]
        x = left + xgap * (donor - 0.5)
        y1, y2 = bottom - exv / axis_max * (bottom - top), bottom - lps / axis_max * (bottom - top)
        parts.append(f'<line x1="{x-offset}" y1="{y1:.1f}" x2="{x+offset}" y2="{y2:.1f}" stroke="#8493a0" stroke-width="2"/>')
        for xpos, y, value, condition in [(x-offset, y1, exv, "ex_vivo"), (x+offset, y2, lps, "LPS")]:
            parts.append(f'<circle cx="{xpos}" cy="{y:.1f}" r="7" fill="{COLORS[condition]}" stroke="#fff" stroke-width="2"/>')
            lx, anchor = (xpos - 12, "end") if condition == "ex_vivo" else (xpos + 12, "start")
            parts.append(f'<text x="{lx}" y="{y+4:.1f}" text-anchor="{anchor}" class="value">{value:.1f}%</text>')
        parts.append(f'<text x="{x}" y="{bottom+29}" text-anchor="middle" class="sample">HC{donor}</text>')
    parts.append('</svg>')
    return "\n".join(parts) + "\n"


def svg_boxplots(result: dict) -> str:
    rows = sorted(result["metrics"], key=lambda m: (m["condition"] != "ex_vivo", m["donor"]))
    width, height = 1780, 900
    left, right, plot_h = 120, 1730, 250
    panels = [("Expressed-feature CPM (log2 scale)", "cpm", 135),
              ("Median-of-ratios normalized counts (log2(count + 1))", "normalized", 515)]
    parts = svg_header(width, height, "Median-of-ratios normalization aligns short-read sample distributions",
                       "Boxes use genes passing CPM ≥ 1 in at least 8 samples and expressed (count > 0) in each sample")
    gap = (right - left) / len(rows)
    for panel_title, kind, top in panels:
        values_by_sample = []
        for j in range(len(rows)):
            if kind == "cpm":
                values = [math.log2(result["cpm"][i][j]) for i in result["passing"]
                          if result["features"][i]["values"][j] > 0]
            else:
                values = [result["log_normalized"][i][j] for i in result["passing"]
                          if result["features"][i]["values"][j] > 0]
            values_by_sample.append(values)
        flat = [value for values in values_by_sample for value in values]
        minimum = min(flat)
        maximum = max(flat)
        axis_min = min(0, math.floor(minimum / 2) * 2)
        axis_max, ticks = nice_axis(math.ceil(maximum / 2) * 2, intervals=6, minimum=axis_min)
        baseline = top + plot_h
        yscale = plot_h / (axis_max - min(axis_min, 0))
        parts.append(f'<text x="{left}" y="{top-18}" class="panel">{panel_title}</text>')
        for tick in ticks:
            y = baseline - (tick - axis_min) * yscale
            parts.extend([f'<line x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}" class="grid"/>',
                          f'<text x="{left-12}" y="{y+4:.1f}" text-anchor="end" class="axis">{tick:g}</text>'])
        for j, (row, values) in enumerate(zip(rows, values_by_sample)):
            q1, med, q3 = quantile(values, .25), quantile(values, .5), quantile(values, .75)
            iqr = q3 - q1
            low = min((v for v in values if v >= q1 - 1.5 * iqr), default=min(values))
            high = max((v for v in values if v <= q3 + 1.5 * iqr), default=max(values))
            x = left + gap * (j + .5)
            color = COLORS[row["condition"]]
            yq1 = baseline - (q1-axis_min)*yscale
            yq3 = baseline - (q3-axis_min)*yscale
            ymed = baseline - (med-axis_min)*yscale
            ylow = baseline - (low-axis_min)*yscale
            yhigh = baseline - (high-axis_min)*yscale
            bw = min(42, gap*.55)
            parts.extend([f'<line x1="{x}" y1="{yhigh:.1f}" x2="{x}" y2="{ylow:.1f}" stroke="{color}" stroke-width="2"/>',
                          f'<line x1="{x-bw/4}" y1="{yhigh:.1f}" x2="{x+bw/4}" y2="{yhigh:.1f}" stroke="{color}" stroke-width="2"/>',
                          f'<line x1="{x-bw/4}" y1="{ylow:.1f}" x2="{x+bw/4}" y2="{ylow:.1f}" stroke="{color}" stroke-width="2"/>',
                          f'<rect x="{x-bw/2}" y="{yq3:.1f}" width="{bw}" height="{max(1,yq1-yq3):.1f}" fill="{color}" fill-opacity=".72" stroke="{color}"/>',
                          f'<line x1="{x-bw/2}" y1="{ymed:.1f}" x2="{x+bw/2}" y2="{ymed:.1f}" stroke="#17232d" stroke-width="2"/>',
                          f'<text x="{x}" y="{baseline+18}" text-anchor="middle" class="sample">HC{row["donor"]}</text>'])
        caption_y = baseline + 44
        for condition, mid in (("ex_vivo", left + gap*4), ("LPS", left + gap*12)):
            label = "Ex vivo" if condition == "ex_vivo" else "3 h LPS"
            parts.extend([f'<rect x="{mid-38}" y="{caption_y-12}" width="12" height="12" rx="2" fill="{COLORS[condition]}"/>',
                          f'<text x="{mid-20}" y="{caption_y-2}" class="group">{label}</text>'])
        divider_x = left + gap*8
        parts.append(f'<line x1="{divider_x}" y1="{top}" x2="{divider_x}" y2="{baseline+25}" stroke="#9aa8b4" stroke-dasharray="4 4"/>')
    parts.extend(['<text x="120" y="887" class="axis">Whiskers extend to the most extreme values within 1.5 × IQR; outliers are omitted.</text>', '</svg>'])
    return "\n".join(parts) + "\n"


def jacobi_eigen(matrix: list[list[float]]) -> tuple[list[float], list[list[float]]]:
    """Eigenvalues and column eigenvectors of a small symmetric matrix."""
    a = [row[:] for row in matrix]
    n = len(a)
    vectors = [[float(i == j) for j in range(n)] for i in range(n)]
    for _ in range(100 * n * n):
        p, q, largest = 0, 1, 0.0
        for i in range(n):
            for j in range(i + 1, n):
                if abs(a[i][j]) > largest:
                    p, q, largest = i, j, abs(a[i][j])
        if largest < 1e-12:
            break
        app, aqq, apq = a[p][p], a[q][q], a[p][q]
        tau = (aqq - app) / (2 * apq)
        t = (1 if tau >= 0 else -1) / (abs(tau) + math.sqrt(1 + tau*tau))
        c = 1 / math.sqrt(1 + t*t)
        s = t * c
        a[p][p] = app - t * apq
        a[q][q] = aqq + t * apq
        a[p][q] = a[q][p] = 0.0
        for k in range(n):
            if k not in (p, q):
                akp, akq = a[k][p], a[k][q]
                a[k][p] = a[p][k] = c*akp - s*akq
                a[k][q] = a[q][k] = s*akp + c*akq
            vkp, vkq = vectors[k][p], vectors[k][q]
            vectors[k][p] = c*vkp - s*vkq
            vectors[k][q] = s*vkp + c*vkq
    order = sorted(range(n), key=lambda i: a[i][i], reverse=True)
    eigenvalues = [max(0.0, a[i][i]) for i in order]
    eigenvectors = [[vectors[row][i] for i in order] for row in range(n)]
    return eigenvalues, eigenvectors


def pca(result: dict) -> dict:
    passing = result["passing"]
    sample_n = len(result["samples"])
    variances = []
    for i in passing:
        values = [result["log_normalized"][i][j] for j in range(sample_n)]
        mean = sum(values) / sample_n
        variances.append((sum((v-mean)**2 for v in values)/(sample_n-1), i))
    selected = [i for _, i in sorted(variances, reverse=True)[:500]]
    centered = []
    for i in selected:
        values = [result["log_normalized"][i][j] for j in range(sample_n)]
        mean = sum(values) / sample_n
        centered.append([v-mean for v in values])
    covariance = [[sum(feature[i]*feature[j] for feature in centered)/(sample_n-1)
                   for j in range(sample_n)] for i in range(sample_n)]
    eigenvalues, eigenvectors = jacobi_eigen(covariance)
    total = sum(eigenvalues)
    scores = [[eigenvectors[j][pc] * math.sqrt(eigenvalues[pc]*(sample_n-1))
               for pc in range(2)] for j in range(sample_n)]
    exv = [i for i,s in enumerate(result["samples"]) if s["condition"] == "ex_vivo"]
    lps = [i for i,s in enumerate(result["samples"]) if s["condition"] == "LPS"]
    if sum(scores[i][0] for i in exv)/len(exv) > sum(scores[i][0] for i in lps)/len(lps):
        for row in scores:
            row[0] *= -1
    return {"scores": scores, "variance": [value/total for value in eigenvalues],
            "selected": selected}


def svg_pca(results: list[dict], pcas: list[dict]) -> str:
    width, height = 1500, 790
    parts = svg_header(width, height, "PCA separates samples by condition in both sequencing assays",
                       "Counts normalized with median-of-ratios; log2(count + 1); 500 most variable CPM-filtered features")
    for result, coords, x0, panel_width in zip(results, pcas, (65, 770), (665, 665)):
        samples, technology = result["samples"], result["technology"]
        plot = (x0+80, x0+panel_width-20, 140, 650)
        xmin, xmax = min(s[0] for s in coords["scores"]), max(s[0] for s in coords["scores"])
        ymin, ymax = min(s[1] for s in coords["scores"]), max(s[1] for s in coords["scores"])
        xpad=max(1,(xmax-xmin)*.16); ypad=max(1,(ymax-ymin)*.2)
        xmin-=xpad; xmax+=xpad; ymin-=ypad; ymax+=ypad
        _, xticks=nice_axis(xmax,5,xmin); _, yticks=nice_axis(ymax,5,ymin)
        xmin, xmax = min(xmin,xticks[0]), max(xmax,xticks[-1]); ymin, ymax=min(ymin,yticks[0]),max(ymax,yticks[-1])
        xleft,xright,ytop,ybottom=plot
        parts.append(f'<text x="{x0+10}" y="105" class="panel">{technology.replace("_"," ").title()} — {len(samples)} samples</text>')
        for v in xticks:
            x=xleft+(v-xmin)/(xmax-xmin)*(xright-xleft)
            parts.extend([f'<line x1="{x:.1f}" y1="{ytop}" x2="{x:.1f}" y2="{ybottom}" class="grid"/>',f'<text x="{x:.1f}" y="{ybottom+22}" text-anchor="middle" class="axis">{v:.0f}</text>'])
        for v in yticks:
            y=ybottom-(v-ymin)/(ymax-ymin)*(ybottom-ytop)
            parts.extend([f'<line x1="{xleft}" y1="{y:.1f}" x2="{xright}" y2="{y:.1f}" class="grid"/>',f'<text x="{xleft-10}" y="{y+4:.1f}" text-anchor="end" class="axis">{v:.0f}</text>'])
        parts.extend([f'<text x="{(xleft+xright)/2}" y="{ybottom+48}" text-anchor="middle" class="axis">PC1 ({coords["variance"][0]*100:.1f}% variance)</text>',
                      f'<text x="{x0+15}" y="{(ytop+ybottom)/2}" transform="rotate(-90 {x0+15} {(ytop+ybottom)/2})" text-anchor="middle" class="axis">PC2 ({coords["variance"][1]*100:.1f}% variance)</text>'])
        for sample, score in zip(samples,coords["scores"]):
            x=xleft+(score[0]-xmin)/(xmax-xmin)*(xright-xleft); y=ybottom-(score[1]-ymin)/(ymax-ymin)*(ybottom-ytop)
            parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="{COLORS[sample["condition"]]}" stroke="#fff" stroke-width="2"/>')
            parts.append(f'<text x="{x+9:.1f}" y="{y-8:.1f}" class="sample">HC{sample["donor"]}</text>')
    parts.extend(['<circle cx="650" cy="735" r="7" fill="#3977a8"/><text x="665" y="740" class="axis">Ex vivo</text>',
                  '<circle cx="760" cy="735" r="7" fill="#d4773c"/><text x="775" y="740" class="axis">3 h LPS</text>',
                  '</svg>'])
    return "\n".join(parts)+"\n"


def correlations(result: dict) -> list[list[float]]:
    rows = [[result["log_normalized"][i][j] for j in range(len(result["samples"]))]
            for i in result["passing"]]
    n = len(result["samples"])
    means = [sum(row[j] for row in rows)/len(rows) for j in range(n)]
    centered = [[row[j]-means[j] for j in range(n)] for row in rows]
    sumsquares = [sum(row[j]**2 for row in centered) for j in range(n)]
    output=[]
    for i in range(n):
        output.append([sum(row[i]*row[j] for row in centered)/math.sqrt(sumsquares[i]*sumsquares[j])
                       for j in range(n)])
    return output


def heat_color(value: float) -> str:
    low=(247,250,252); high=(38,92,132)
    t=max(0,min(1,(value-.35)/.65))
    rgb=tuple(round(a+(b-a)*t) for a,b in zip(low,high))
    return '#%02x%02x%02x'%rgb


def svg_heatmaps(results: list[dict], matrices: list[list[list[float]]]) -> str:
    width,height=1600,900
    parts=svg_header(width,height,"Samples correlate most strongly within their condition",
                     "Pearson correlation of log2(median-of-ratios normalized counts + 1) across all CPM-filtered features")
    for result,matrix,x0,cell in zip(results,matrices,(55,870),(30,46)):
        samples=sorted(result["samples"],key=lambda s:(s["condition"]!="ex_vivo",s["donor"]))
        original=result["samples"]
        order=[original.index(s) for s in samples]
        matrix=[[matrix[i][j] for j in order] for i in order]
        n=len(samples); xstart=x0+100; ystart=155
        parts.append(f'<text x="{x0+12}" y="112" class="panel">{result["technology"].replace("_"," ").title()} ({n} samples)</text>')
        parts.append(f'<text x="{xstart}" y="140" class="group" fill="{COLORS["ex_vivo"]}">Ex vivo</text>')
        parts.append(f'<text x="{xstart+cell*sum(s["condition"]=="ex_vivo" for s in samples)}" y="140" class="group" fill="{COLORS["LPS"]}">3 h LPS</text>')
        for i,sample in enumerate(samples):
            short=f'HC{sample["donor"]} '+('E' if sample["condition"]=="ex_vivo" else 'L')
            for j in range(n):
                x=xstart+j*cell; y=ystart+i*cell
                value=matrix[i][j]
                parts.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{heat_color(value)}" stroke="#fff" stroke-width="1"/>')
                if cell>=40:
                    parts.append(f'<text x="{x+cell/2}" y="{y+cell/2+4}" text-anchor="middle" font-size="11" fill="{("#fff" if value>.72 else "#243342")}">{value:.2f}</text>')
            y=ystart+i*cell+cell/2+4
            parts.append(f'<text x="{xstart-10}" y="{y}" text-anchor="end" font-size="11" fill="{COLORS[sample["condition"]]}">{short}</text>')
        for j,sample in enumerate(samples):
            x=xstart+j*cell+cell/2
            y=ystart+n*cell+15
            short=f'HC{sample["donor"]} '+('E' if sample["condition"]=="ex_vivo" else 'L')
            parts.append(f'<text x="{x}" y="{y}" transform="rotate(55 {x} {y})" text-anchor="start" font-size="11" fill="{COLORS[sample["condition"]]}">{short}</text>')
        legend_y=ystart+n*cell+95
        for step in range(7):
            value=.4+step*.1; x=xstart+step*44
            parts.append(f'<rect x="{x}" y="{legend_y}" width="44" height="12" fill="{heat_color(value)}"/>')
        parts.append(f'<text x="{xstart}" y="{legend_y+29}" class="axis">Pearson r: 0.4</text><text x="{xstart+264}" y="{legend_y+29}" text-anchor="end" class="axis">1.0</text>')
    parts.extend(['<text x="60" y="860" class="axis">Heatmap labels: E = ex vivo; L = 3 h LPS. Color indicates correlation strength.</text>','</svg>'])
    return "\n".join(parts)+"\n"


def main() -> None:
    metadata=read_metadata()
    results=[]
    for technology,path in TABLES.items():
        samples,features=read_counts(path)
        results.append(analyze(technology,samples,features,metadata))
    pcas=[pca(result) for result in results]
    corrs=[correlations(result) for result in results]
    FIGURES.mkdir(parents=True,exist_ok=True)
    (FIGURES/"short_read_top100_gene_share.svg").write_text(svg_top100_chart(results[0]),encoding="utf-8")
    (FIGURES/"short_read_normalization_boxplots.svg").write_text(svg_boxplots(results[0]),encoding="utf-8")
    (FIGURES/"normalized_counts_pca.svg").write_text(svg_pca(results,pcas),encoding="utf-8")
    (FIGURES/"sample_correlation_heatmaps.svg").write_text(svg_heatmaps(results,corrs),encoding="utf-8")
    fields=["technology","sample_column","geo_accession","donor","condition","features","library_total","cpm_filter_features","median_expressed_cpm","size_factor","median_log2_normalized_expressed","top100_share_percent"]
    with SUMMARY.open("w",encoding="utf-8",newline="") as handle:
        writer=csv.DictWriter(handle,fields,delimiter="\t",lineterminator="\n");writer.writeheader()
        for result in results: writer.writerows(result["metrics"])
    for result,coords in zip(results,pcas):
        print(f'{result["technology"]}: CPM-filtered={len(result["passing"]):,}; PC1={coords["variance"][0]*100:.2f}%; PC2={coords["variance"][1]*100:.2f}%')
    print(f"Wrote four figures to {FIGURES.relative_to(ROOT)} and per-sample values to {SUMMARY.relative_to(ROOT)}")


if __name__=="__main__":
    main()
