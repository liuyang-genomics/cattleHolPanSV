#!/bin/bash
# fasta_stats_gfastats.sh
# Compute assembly statistics using gfastats for all relevant jer_* assemblies

set -euo pipefail

# Directory containing input assemblies
input_dir="${PANEL_DIR}"
output_dir="${PROJECT_ROOT}/stat_pan/0.fa.stats"
mkdir -p "$output_dir"

# Run gfastats
find "$input_dir" -name "jer_*.bp.p_ctg.fa" ! -name "*bak*" | while read fa; do
    id=$(basename "$fa" .p_ctg.fa)
    id=${id/.hifi_hi2c.hic/}
    gfastats "$fa" > "${output_dir}/${id}.gfastats"
done

# Optional: summary views
grep "Total scaffold length" ${output_dir}/*.gfastats
grep "Contig N50:" ${output_dir}/*.gfastats
