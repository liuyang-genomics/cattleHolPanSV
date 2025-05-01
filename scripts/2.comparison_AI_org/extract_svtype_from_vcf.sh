#!/bin/bash
# extract_svtype_from_vcf.sh
# Extract and summarize SV types and lengths from a list of VCFs

set -euo pipefail

input_dir="${PROJECT_ROOT}/pangenie_HiFi/3.pan_all"
output_dir="${input_dir}/summary"
mkdir -p "$output_dir"

awk_script="/mnt/data/summarize_svtype.awk"

for vcf in "$input_dir"/*.vcf.gz; do
    filename=$(basename "$vcf")
    sample=${filename%%-*}
    allele=${filename#*.}
    allele=${allele%%.*}
    id="${sample}.${allele}"

    # Extract SV info from VCF
    zcat "$vcf" | bcftools query -f "%ID %INFO/ID\n" | \
        awk '{
            if($2 ~ /,|:/) {
                split($2,a,/,|:/);
                for(i in a) {
                    gsub("-"," ",a[i]);
                    print $1,a[i]
                }
            } else {
                gsub("-"," ",$2);
                print
            }
        }' | awk '{print $1,$2,$3,$5,$4,$6}' | awk -v idd="$id" -f "$awk_script" > "$output_dir/${id}.stats"
done
