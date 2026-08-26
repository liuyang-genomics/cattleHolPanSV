#!/bin/bash
# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PANEL_DIR:?PANEL_DIR is unset - see config.sh.example at the repository root}"
# --------------------------

# summarize_assembly_metrics.sh
# Summarize BUSCO, QUAST, and gfastats metrics for all assemblies

set -euo pipefail

panel_dir="${PANEL_DIR}"

# Summarize BUSCO completeness
echo -e "Sample\tBUSCO_Complete"
for file in "$panel_dir"/*/2.assembly/*.busco-hifiasm-bp/short_summary*.json; do
    sample=$(basename "$(dirname "$file")")
    value=$(grep "C:" "$file")
    echo -e "${sample}\t${value}"
done > busco_completeness.tsv

# Summarize QUAST N50
echo -e "Sample\tQUAST_Contig_N50"
for file in "$panel_dir"/*/2.assembly/*.busco-hifiasm-*/short_summary*.json; do
    sample=$(basename "$(dirname "$file")")
    value=$(grep "Contigs N50" "$file" | grep -v "of")
    echo -e "${sample}\t${value}"
done > quast_n50.tsv

# Summarize QUAST total length
echo -e "Sample\tQUAST_Total_Length"
for file in "$panel_dir"/*/2.assembly/*.busco-hifiasm-*/short_summary*.json; do
    sample=$(basename "$(dirname "$file")")
    value=$(grep "Total length" "$file" | grep -v "of")
    echo -e "${sample}\t${value}"
done > quast_total_length.tsv

# Summarize gfastats contig N50
echo -e "Sample\tGfastats_Contig_N50"
for file in *hifi_hi2c.hic*.fa.gfastats; do
    sample=$(basename "$file" .p_ctg.fa.gfastats)
    value=$(grep "Contig N50" "$file" | grep -v "of")
    echo -e "${sample}\t${value}"
done > gfastats_n50.tsv

# Summarize gfastats scaffold total length
echo -e "Sample\tGfastats_Total_Scaffold_Length"
for file in *hifi_hi2c.hic*.fa.gfastats; do
    sample=$(basename "$file" .p_ctg.fa.gfastats)
    value=$(grep "Total scaffold length" "$file" | grep -v "of")
    echo -e "${sample}\t${value}"
done > gfastats_total_length.tsv
