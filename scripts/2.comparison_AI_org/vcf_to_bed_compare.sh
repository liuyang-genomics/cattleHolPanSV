#!/bin/bash
# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
# --------------------------

# vcf_to_bed_compare.sh
# Extract BED-like info from VCFs for comparison (first 5 columns)

set -euo pipefail

cd ${PROJECT_ROOT}/comp_pan
mkdir -p 3.pan_vcf_compa

declare -A vcf_map=(
  ["allbovinePan.hol.filter-bi"]="2.vcf_filter/1.grp/allbovinePan-2024-07-03.hol.filter-bi.vcf.gz"
  ["hol.filter-bi"]="2.vcf_filter/0.filterF_M0.2/hol-pg2hic-2024-05-22.filter-bi.vcf.gz"
  ["allbovinePan.jer.filter-bi"]="2.vcf_filter/1.grp/allbovinePan-2024-07-03.jer.filtered.vcf.gz"
  ["jer.filter-bi"]="2.vcf_filter/0.filterF_M0.2/jer-pg-2024-05-15.filter-bi.vcf.gz"
)

for name in "${!vcf_map[@]}"; do
  vcf_path="${vcf_map[$name]}"
  out_path="3.pan_vcf_compa/${name}.bed"

  zcat "$vcf_path" | grep -v "^##" | cut -f1-5 > "$out_path"
done
