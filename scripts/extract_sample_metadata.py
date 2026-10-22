import gzip
import re
import csv
from pathlib import Path

source = Path("metadata/GSE297994_family.soft.gz")
output = Path("metadata/sample_metadata.tsv")

samples = []
current = None

with gzip.open(source, "rt", encoding="utf-8") as handle:
    for line in handle:
        line = line.rstrip("\n")

        if line.startswith("^SAMPLE = "):
            if current:
                samples.append(current)
            current = {
                "geo_accession": line.split("=", 1)[1].strip(),
                "title": "",
                "treatment": "",
                "cell_type": "",
            }

        elif current and line.startswith("!Sample_title = "):
            current["title"] = line.split("=", 1)[1].strip()

        elif current and line.startswith("!Sample_characteristics_ch1 = "):
            value = line.split("=", 1)[1].strip()
            if value.startswith("treatment:"):
                current["treatment"] = value.split(":", 1)[1].strip()
            elif value.startswith("cell type:"):
                current["cell_type"] = value.split(":", 1)[1].strip()

if current:
    samples.append(current)

for sample in samples:
    title = sample["title"].lower()
    sample["technology"] = (
        "long-read" if "long-read" in title else "short-read"
    )
    sample["condition"] = (
        "LPS" if "lps" in title else "ex_vivo"
    )

with output.open("w", newline="", encoding="utf-8") as handle:
    fields = [
        "geo_accession", "title", "technology",
        "condition", "treatment", "cell_type"
    ]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(samples)

print(f"Samples extracted: {len(samples)}")
print(f"Saved to: {output}")
print("\nCounts by technology and condition:")
counts = {}
for sample in samples:
    key = (sample["technology"], sample["condition"])
    counts[key] = counts.get(key, 0) + 1

for key, count in sorted(counts.items()):
    print(f"{key[0]:10s} | {key[1]:8s} | {count}")

print("\nFirst five records:")
with output.open(encoding="utf-8") as handle:
    for _, line in zip(range(6), handle):
        print(line.rstrip())
