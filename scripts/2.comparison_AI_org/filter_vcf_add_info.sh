#!/bin/bash
# filter_vcf_add_info.sh
# Sort, index, annotate, and filter VCFs for SVTYPE/SVLEN > 50bp

set -euo pipefail

# Paths
input_dir="${PROJECT_ROOT}/stat_pan/2.vcf_stats/2.ind"
output_dir="${PROJECT_ROOT}/stat_pan/4.truvari_sv/0.sv_vcfgs"
tmp_dir="${PROJECT_ROOT}/stat_pan/4.truvari_sv/0.tmp"
awk_script="/mnt/data/pangenie_remove50.awk"
chr_list=$(seq 1 29 | tr '\n' ',' && echo X)

mkdir -p "$output_dir" "$tmp_dir"

# Process each sample VCF
for vcf in "$input_dir"/hol-pg2hic-2024-05-22.sample_*.anno_biallelic.filtered.vcf.gz; do
    sample=$(basename "$vcf" .anno_biallelic.filtered.vcf.gz | sed 's/hol-pg2hic-2024-05-22.sample_//')

    # Step 1: Sort and index
    sorted_vcf="$tmp_dir/${sample}.vcf.gz"
    bcftools sort "$vcf" -Oz -o "$sorted_vcf"
    tabix -f -p vcf "$sorted_vcf"

    # Step 2: Extract chr1–29,X and filter + annotate
    bcftools view --regions "$chr_list" "$sorted_vcf" |
        awk -f "$awk_script" |
        bgzip -c > "$output_dir/${sample}.sv.vcf.gz"

    tabix -f -p vcf "$output_dir/${sample}.sv.vcf.gz"
done
