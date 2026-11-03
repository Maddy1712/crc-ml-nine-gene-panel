#!/usr/bin/env bash
# Downloads raw FASTQ files for the three ENA development cohorts.
# TEMPLATE: check accession-to-cohort order and file naming before use.
set -euo pipefail
OUT="${1:-fastq}"
for ACC in PRJNA596460 PRJNA658001 PRJNA218851; do
  mkdir -p "${OUT}/${ACC}"
  wget -q -O "${OUT}/${ACC}/filereport.tsv" \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=run_accession,fastq_ftp,sample_title&format=tsv"
  tail -n +2 "${OUT}/${ACC}/filereport.tsv" | cut -f2 | tr ';' '\n' | grep -v '^$' |
    while read -r url; do wget -c -P "${OUT}/${ACC}" "ftp://${url}"; done
done
