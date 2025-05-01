#!/bin/bash
# summarize_rtg_results.sh
# Extract summary statistics from RTG vcfeval results

set -euo pipefail

input_dir="${PROJECT_ROOT}/stat_pan/5.rtg_snv/1.rtg"
output_file="${input_dir}/rtg_summary_all.tsv"

echo -e "Sample\tBase\tCall\tRegion\tPrecision\tRecall\tF1\tTP\tFP\tFN" > "$output_file"

for path in "$input_dir"/*/summary.txt; do
    dir=$(basename "$(dirname "$path")")  # e.g., sample.base.call.region.rtg
    parts=(${dir//./ })
    sample="${parts[0]}"
    base="${parts[1]}"
    call="${parts[2]}"
    region="${parts[3]}"

    data=$(awk '/^Overall/ {getline; print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' "$path")
    echo -e "${sample}\t${base}\t${call}\t${region}\t${data}" >> "$output_file"
done
