#!/bin/bash
# Download the published HuRI interaction table on a Biowulf login node.

set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output="${1:-$repo_root/runtime/reference/huri/HuRI.tsv}"
url="${GBI_HURI_URL:-https://www.interactome-atlas.org/data/HuRI.tsv}"
mkdir -p "$(dirname "$output")"
temporary="${output}.download"
trap 'rm -f -- "$temporary"' EXIT
curl --fail --location --retry 3 --output "$temporary" "$url"

python3 -c 'import re,sys; p=sys.argv[1]; n=sum(1 for line in open(p,encoding="utf-8-sig") if len(line.rstrip("\n").split("\t"))>=2 and re.fullmatch(r"ENSG\d+(?:\.\d+)?",line.split("\t")[0].strip(),re.I)); assert n>0, "download has no HuRI Ensembl pairs"; print(f"Validated {n:,} HuRI rows")' "$temporary"
mv -f "$temporary" "$output"
trap - EXIT
echo "Saved HuRI to $output"
