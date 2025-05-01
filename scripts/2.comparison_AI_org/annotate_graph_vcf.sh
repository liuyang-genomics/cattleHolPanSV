#!/bin/bash
# annotate_graph_vcf.sh
# Annotate graph-based VCFs using their corresponding GFA files

set -euo pipefail

# Settings
annotator="${SOFTWARE_DIR}/genotyping-pipelines/prepare-vcf-MC/workflow/scripts/annotate_vcf.py"
vcf_list=$(find ${PROJECT_ROOT}/minigraph-cactus -name "*.vcf.gz" | grep -vE "raw|filter|bak|bovinePan|hol")

cd ${PROJECT_ROOT}/stat_pan/2.vcf_stats/0.filter

# Annotate
for vcf in $vcf_list; do
    prefix=${vcf%.vcf.gz}
    vcf_unzipped="${prefix}.vcf"
    gfa_file="${prefix}.gfa"
    anno_file="${prefix}.anno.vcf"
    log_file="${prefix}.anno.log"

    gunzip -c "$vcf" > "$vcf_unzipped"
    gunzip -c "$gfa_file.gz" > "$gfa_file"

    python3 "$annotator" -vcf "$vcf_unzipped" -gfa "$gfa_file" -o "$anno_file" &> "$log_file"

    bgzip "$anno_file"
    bcftools sort -Oz -o "${anno_file}.gz" "${anno_file}.gz"
    tabix -f -p vcf "${anno_file}.gz"

    rm "$vcf_unzipped" "$gfa_file"
done
