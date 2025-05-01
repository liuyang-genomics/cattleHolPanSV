#!/bin/bash
# merge_assembly_summaries.sh
# Merge BUSCO, QUAST, and gfastats summaries into one combined table

set -euo pipefail

# Load individual files
paste busco_completeness.tsv quast_n50.tsv quast_total_length.tsv gfastats_n50.tsv gfastats_total_length.tsv | awk '
BEGIN {OFS="\t"}
NR==1 {
    print "Sample", "BUSCO_Complete", "QUAST_Contig_N50", "QUAST_Total_Length", "Gfastats_Contig_N50", "Gfastats_Total_Scaffold_Length"
}
NR>1 {
    print $1, $2, $4, $6, $8, $10
}
' > combined_assembly_summary.tsv

echo "✅ Combined summary saved as combined_assembly_summary.tsv"
