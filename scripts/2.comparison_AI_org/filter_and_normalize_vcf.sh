#!/bin/bash
# filter_and_normalize_vcf.sh
# Apply bcftools filters and normalize VCFs to biallelic

set -euo pipefail

input_vcfs=$(ls ${PROJECT_ROOT}/minigraph-cactus/*/*.vcf.gz | grep -v 'raw')

# Step 1: Apply missingness filter
for vcf in $input_vcfs; do
    out_vcf="${vcf/.vcf.gz/.filtered.vcf.gz}"
    bcftools view --threads 8 -i 'F_MISSING<0.2' "$vcf" -Oz -o "$out_vcf"
done

# Step 2: Normalize to biallelic
for filtered_vcf in *.filtered.vcf.gz; do
    out_biallelic="${filtered_vcf/filtered/filtered-bi}"
    bcftools norm --threads 8 -m- "$filtered_vcf" -Oz -o "$out_biallelic"
done

# Step 3: Subset by sample lists (assuming lists exist)
bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -S hol.sample allbovinePan-2024-07-03.filtered-bi.vcf.gz -Oz -o allbovinePan-2024-07-03_hol.filtered-bi.vcf.gz
bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -S jer.sample allbovinePan-2024-07-03.filtered-bi.vcf.gz -Oz -o allbovinePan-2024-07-03_jer.filtered-bi.vcf.gz
bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -S pub.sample allbovinePan-2024-07-03.filtered-bi.vcf.gz -Oz -o allbovinePan-2024-07-03_pub.filtered-bi.vcf.gz
