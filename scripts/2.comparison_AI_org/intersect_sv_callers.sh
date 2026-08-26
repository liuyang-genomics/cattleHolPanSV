#!/bin/bash
# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PANEL_DIR:?PANEL_DIR is unset - see config.sh.example at the repository root}"
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
# --------------------------

# intersect_sv_callers.sh
# Intersect pangenome SVs with tool-based SV calls, followed by SV count/length/share statistics

set -euo pipefail

base_dir="${PROJECT_ROOT}/stat_pan/3.intersect.ab_rb_sr"
mkdir -p "$base_dir"/{0.panBed,1.ab_rb_srBed,2.sv.intersect_both,6.stats-jer,7.share_stats-jer}

# Function: extract pan-SV BED
generate_pan_bed() {
  zcat "$1" | bcftools query -f "%CHROM %POS %ID %INFO/ID\n" |
  awk '{
    split($4, bc, "-")
    print $1, $2, $2 + bc[5], bc[5], bc[3], $3, bc[4]
  }' | awk '($1 <= 29 || $1 == "X") && $4 > 50 && $4 <= 1000000' |
  sed 's/ /\t/g'
}

# Function: extract tool-SV BED
generate_tool_bed() {
  zcat "$1" | bcftools query -f "%CHROM %POS %ID %INFO/SVTYPE %INFO/SVLEN\n" |
  awk '{
    if ($5 == "." || $5 == "") $5 = 1
    if ($5 < 0) $5 = -$5
    print $1, $2, $2 + $5, $5, $4, "-", $3
  }' | awk '($1 <= 29 || $1 == "X") && $4 > 50 && $4 <= 1000000' |
  sed 's/ /\t/g'
}

# Input VCFs
pan_vcfs=$(ls ${PROJECT_ROOT}/stat_pan/2.vcf_stats/2.ind/hol-*.vcf.gz | grep -v "4439")
tool_vcfs=$(ls ${PANEL_DIR}/sample_*/4.sv_vcf/sample_*.vcf.gz)

# Generate pan BEDs
for vcf in $pan_vcfs; do
  id=$(basename "$vcf" .anno_biallelic.filtered.vcf.gz | sed 's/.*sample_//')
  generate_pan_bed "$vcf" > "$base_dir/0.panBed/${id}.bed"
done

# Generate tool BEDs
for vcf in $tool_vcfs; do
  id=$(basename "$vcf" .vcf.gz)
  generate_tool_bed "$vcf" > "$base_dir/1.ab_rb_srBed/${id}.bed"
done

# Intersect and compute statistics
for pan_bed in "$base_dir"/0.panBed/*.bed; do
  sid=$(basename "$pan_bed" .bed)
  for tool_bed in "$base_dir"/1.ab_rb_srBed/${sid}*.bed; do
    [ -f "$tool_bed" ] || continue
    tid=$(basename "$tool_bed" .bed)

    intersect_out="$base_dir/2.sv.intersect_both/${tid}_pan.bed"
    bedtools intersect -a "$pan_bed" -b "$tool_bed" -f 0.9 -r -wo > "$intersect_out"

    awk -f "$base_dir/pan_stats.awk" "$tool_bed" > "$base_dir/6.stats-jer/${tid}.stats"
    awk -f "$base_dir/pan_share_stats.awk" "$intersect_out" > "$base_dir/7.share_stats-jer/${tid}.sv.share.stats"
  done
done
