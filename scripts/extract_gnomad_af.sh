#!/usr/bin/env bash
set -euo pipefail

# Stream gnomAD v4 genome AF fields from UCSC without downloading full VCFs.
#
# Usage:
#   OUTDIR=data/interim/gnomad_v4_AF bash scripts/extract_gnomad_af.sh
#   sbatch --array=1-24 scripts/extract_gnomad_af.sh

BASE="${BASE:-https://hgdownload.soe.ucsc.edu/gbdb/hg38/gnomAD/v4/genomes}"
OUTDIR="${OUTDIR:-data/interim/gnomad_v4_AF}"
CHROMOSOMES="${CHROMOSOMES:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y}"

mkdir -p "${OUTDIR}/logs"

if command -v module >/dev/null 2>&1; then
    module load bcftools 2>/dev/null || module add bcftools 2>/dev/null || true
fi

if ! command -v bcftools >/dev/null 2>&1; then
    echo "ERROR: bcftools is required but was not found on PATH." >&2
    exit 1
fi

read -r -a CHROMS <<< "${CHROMOSOMES}"

extract_one() {
    local chr="$1"
    local out="${OUTDIR}/gnomad_v4_AF_chr${chr}.txt.gz"
    local url="${BASE}/gnomad.genomes.v4.0.sites.chr${chr}.vcf.bgz"

    echo "$(date) Starting chr${chr}"
    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AF_grpmax\n' \
        "${url}" | gzip > "${out}"
    echo "$(date) Finished chr${chr}: ${out}"
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    idx=$((SLURM_ARRAY_TASK_ID - 1))
    if (( idx < 0 || idx >= ${#CHROMS[@]} )); then
        echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} is outside 1-${#CHROMS[@]}." >&2
        exit 1
    fi
    extract_one "${CHROMS[$idx]}"
else
    for chr in "${CHROMS[@]}"; do
        extract_one "${chr}"
    done
fi

