#!/usr/bin/env bash
# Download the public CosMx SMI dataset GSE234713.
# 9 human colon samples; see the GEO record for the sample table.
# NOTE: the full raw download is large (tens of GB). Per-sample files are fetched individually.
set -euo pipefail
OUT="${1:-./RAW}"
mkdir -p "$OUT"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234713/suppl"

echo "Downloading GSE234713 supplementary files into $OUT"
echo "See: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234713"
echo
echo "Fetch the file list from GEO and download the *_exprMat_file.csv.gz and"
echo "*_metadata_file.csv.gz for each GSM, plus the CellComposite image folders."
echo "Example:"
echo "  curl -O ${BASE}/GSE234713_RAW.tar   # full archive"
