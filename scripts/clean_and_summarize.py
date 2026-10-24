#!/usr/bin/env python3
"""Inspect and standardize GSE297994 processed TSV tables (Python 3.9+, stdlib only)."""
from __future__ import annotations
import csv, gzip, json, math, sys
from collections import Counter
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = ROOT / "data" / "processed"
CLEAN_DIR = ROOT / "data" / "cleaned"
RESULTS_DIR = ROOT / "results"
FILES = ("long_read_raw_counts.txt.gz", "long_read_salmon_abundance.txt.gz", "short_read_raw_counts.txt.gz", "short_read_salmon_abundance.txt.gz")

def clean_table(filename: str) -> dict:
    source = SOURCE_DIR / filename
    if not source.is_file():
        raise FileNotFoundError(f"Missing source table: {source}")
    with gzip.open(source, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        raw_header = next(reader, None)
        if not raw_header:
            raise ValueError(f"Empty table: {filename}")
        header = [cell.strip() for cell in raw_header]
        if not header[0]:
            header[0] = "feature_id"
        if any(not cell for cell in header[1:]):
            raise ValueError(f"Blank sample column name in {filename}")
        duplicate_columns = [name for name, n in Counter(header).items() if n > 1]
        if duplicate_columns:
            raise ValueError(f"Duplicate column names in {filename}: {duplicate_columns}")
        # In the native long-read raw-count table, `gene_id` repeats because
        # several transcript rows can belong to one gene. The separate `ids`
        # column carries the row-level transcript identifier.
        feature_id_index = header.index("ids") if "ids" in header else 0
        CLEAN_DIR.mkdir(parents=True, exist_ok=True)
        output = CLEAN_DIR / filename.removesuffix(".gz")
        rows = blank_rows = empty_cells = duplicate_ids = missing_ids = 0
        seen_ids = set()
        seen_first_column = set()
        repeated_first_column_values = 0
        numeric = [dict(n=0, total=0.0, min=None, max=None, zero=0, negative=0, non_numeric=0) for _ in header[1:]]
        annotation_values = [None] * len(header)
        nonempty = [0] * len(header)
        with output.open("w", encoding="utf-8", newline="") as out_handle:
            writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            for line_number, raw_row in enumerate(reader, start=2):
                row = [cell.strip() for cell in raw_row]
                if not any(row):
                    blank_rows += 1
                    continue
                if len(row) != len(header):
                    raise ValueError(f"{filename}:{line_number}: expected {len(header)} fields, found {len(row)}")
                rows += 1
                empty_cells += sum(cell == "" for cell in row)
                for i, cell in enumerate(row):
                    nonempty[i] += bool(cell)
                feature_id = row[feature_id_index]
                if not feature_id:
                    missing_ids += 1
                if feature_id in seen_ids:
                    duplicate_ids += 1
                seen_ids.add(feature_id)
                if row[0] in seen_first_column:
                    repeated_first_column_values += 1
                seen_first_column.add(row[0])
                for i, text in enumerate(row[1:]):
                    if not text:
                        continue
                    try:
                        value = float(text)
                    except ValueError:
                        numeric[i]["non_numeric"] += 1
                        if annotation_values[i + 1] is None:
                            annotation_values[i + 1] = set()
                        annotation_values[i + 1].add(text)
                        continue
                    if not math.isfinite(value):
                        numeric[i]["non_numeric"] += 1
                        if annotation_values[i + 1] is None:
                            annotation_values[i + 1] = set()
                        annotation_values[i + 1].add(text)
                        continue
                    stat = numeric[i]
                    stat["n"] += 1
                    stat["total"] += value
                    stat["min"] = value if stat["min"] is None else min(stat["min"], value)
                    stat["max"] = value if stat["max"] is None else max(stat["max"], value)
                    stat["zero"] += value == 0
                    stat["negative"] += value < 0
                writer.writerow(row)
    column_summaries = []
    annotation_columns = []
    for name, stat, filled in zip(header[1:], numeric, nonempty[1:]):
        if stat["n"] == 0 and stat["non_numeric"] > 0:
            annotation_columns.append({
                "column": name,
                "nonempty_cells": filled,
                "distinct_non_numeric_values": len(annotation_values[header.index(name)] or set()),
                "selected_as_feature_id": header.index(name) == feature_id_index,
            })
            continue
        total = stat.pop("total")
        stat["mean"] = total / stat["n"] if stat["n"] else None
        column_summaries.append({"column": name, "nonempty_cells": filled, **stat})
    return {"source": source.relative_to(ROOT).as_posix(), "cleaned_copy": output.relative_to(ROOT).as_posix(),
            "rows": rows, "columns": len(header), "feature_id_column": header[feature_id_index],
            "duplicate_feature_ids": duplicate_ids, "missing_feature_ids": missing_ids, "blank_rows_removed": blank_rows,
            "repeated_first_column_values": repeated_first_column_values,
            "empty_cells": empty_cells, "columns_detail": [{"name": n, "nonempty_cells": k} for n,k in zip(header,nonempty)],
            "annotation_columns": annotation_columns,
            "numeric_column_summaries": column_summaries,
            "cleaning_policy": "Trim surrounding whitespace and omit fully blank rows. Preserve source files, duplicate IDs, and empty cells. Reject inconsistent row widths."}

def main() -> int:
    try:
        files = [clean_table(name) for name in FILES]
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    report = {"dataset": "NCBI GEO GSE297994", "generated_by": "scripts/clean_and_summarize.py", "files": files}
    path = RESULTS_DIR / "task2_data_summary.json"
    path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {path.relative_to(ROOT)}")
    for item in files:
        print(f"{Path(item['source']).name}: {item['rows']} rows, {item['columns']} columns, {item['empty_cells']} empty cells, {item['duplicate_feature_ids']} duplicate IDs")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
