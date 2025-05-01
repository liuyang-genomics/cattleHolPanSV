#!/bin/bash
# summarize_truvari_results.sh
# Extract precision, recall, F1-score, and counts from truvari summary.tsv files

set -euo pipefail

result_dir="${PROJECT_ROOT}/stat_pan/4.truvari_sv/1.benchmark"
output_tsv="${result_dir}/truvari_summary_all.tsv"

echo -e "Sample\tRegion\tPrecision\tRecall\tF1\tTP\tFP\tFN" > "$output_tsv"

for summary in "$result_dir"/*/summary.tsv; do
    sample=$(basename "$(dirname "$summary")")
    values=$(awk -F '\t' 'NR==2 {print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' "$summary")
    echo -e "${sample//./\t}\t$values" >> "$output_tsv"
done
