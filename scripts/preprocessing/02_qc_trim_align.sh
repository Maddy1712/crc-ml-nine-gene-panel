#!/usr/bin/env bash
# QC -> trimming -> alignment, run once per cohort directory.  Usage:
#   bash 02_qc_trim_align.sh fastq/PRJNA658001 bam/PRJNA658001
# Parameters follow the Methods (Sections 2.2-2.3). TEMPLATE: set HISAT2_INDEX and
# ADAPTERS, and keep the Trimmomatic step order identical to what you actually ran.
set -euo pipefail
IN="$1"; OUT="$2"
THREADS=16
HISAT2_INDEX="${HISAT2_INDEX:?set path to the GRCh38 HISAT2 index prefix}"
ADAPTERS="${ADAPTERS:?set path to the adapter FASTA (e.g. TruSeq3-PE.fa)}"
mkdir -p "${OUT}/fastqc_raw" "${OUT}/trimmed" "${OUT}/fastqc_trimmed" "${OUT}/bam"

fastqc -t "$THREADS" -o "${OUT}/fastqc_raw" "${IN}"/*.fastq.gz
multiqc -o "${OUT}/fastqc_raw" "${OUT}/fastqc_raw"

for R1 in "${IN}"/*_1.fastq.gz; do
  S=$(basename "$R1" _1.fastq.gz); R2="${IN}/${S}_2.fastq.gz"
  trimmomatic PE -threads "$THREADS" "$R1" "$R2" \
    "${OUT}/trimmed/${S}_1P.fq.gz" "${OUT}/trimmed/${S}_1U.fq.gz" \
    "${OUT}/trimmed/${S}_2P.fq.gz" "${OUT}/trimmed/${S}_2U.fq.gz" \
    ILLUMINACLIP:"${ADAPTERS}":2:30:10 HEADCROP:12 CROP:87 MINLEN:50
done

fastqc -t "$THREADS" -o "${OUT}/fastqc_trimmed" "${OUT}"/trimmed/*P.fq.gz
multiqc -o "${OUT}/fastqc_trimmed" "${OUT}/fastqc_trimmed"

for P1 in "${OUT}"/trimmed/*_1P.fq.gz; do
  S=$(basename "$P1" _1P.fq.gz)
  hisat2 -p "$THREADS" -x "$HISAT2_INDEX" -1 "$P1" -2 "${OUT}/trimmed/${S}_2P.fq.gz" |
    samtools sort -@ "$THREADS" -o "${OUT}/bam/${S}.bam" -
  samtools index "${OUT}/bam/${S}.bam"
done
