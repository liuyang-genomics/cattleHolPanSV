#!/bin/bash
# run_holstein_sv_class.sh
# Generate SV classification table for Holstein vs Public based on genotype info

set -euo pipefail

input_vcf="bovHol-2024-08-13.anno_biallelic.filtered.vcf.gz"
output_tsv="bovHol-2024-08-13.anno_biallelic.classified.tsv"
awk_script="classify_holstein_sv.awk"

zcat "$input_vcf" | awk -f "$awk_script" > "$output_tsv"
