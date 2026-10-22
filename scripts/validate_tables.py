import csv
import gzip
from collections import Counter
from pathlib import Path

files = sorted(Path("data/processed").glob("*.gz"))

for path in files:
    print("\n" + "=" * 70)
    print(f"FILE: {path}")
    print("=" * 70)

    with gzip.open(path, "rt", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)

        if not header:
            print("ERROR: Empty file")
            continue

        print(f"Header columns: {len(header)}")
        print(f"Header fields: {header}")
        print(f"First data rows:")

        row_widths = Counter()
        first_rows = []
        feature_ids = []
        total_rows = 0
        empty_cells = 0

        for row in reader:
            total_rows += 1
            row_widths[len(row)] += 1
            empty_cells += sum(cell == "" for cell in row)

            if len(first_rows) < 2:
                first_rows.append(row)

            if row:
                feature_ids.append(row[0])

        print(f"Data rows: {total_rows}")
        print(f"Column counts across rows: {dict(row_widths)}")
        print(f"Empty cells: {empty_cells}")
        print(f"Unique first-column IDs: {len(set(feature_ids))}")
        print(f"Duplicate first-column IDs: {len(feature_ids) - len(set(feature_ids))}")

        for i, row in enumerate(first_rows, start=1):
            print(f"Row {i}: {row[:5]}")

        if len(row_widths) == 1 and next(iter(row_widths)) == len(header):
            print("STRUCTURE CHECK: PASS")
        else:
            print("STRUCTURE CHECK: REVIEW REQUIRED")
