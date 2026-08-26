#!/bin/bash
# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
# --------------------------

# run_truvari_benchmark.sh
# Benchmark SV calls using truvari against pangenome reference VCF

set -euo pipefail

ref_vcf="${PROJECT_ROOT}/stat_pan/4.truvari_sv/reference/holstein_pan_merged.vcf.gz"
callset_dir="${PROJECT_ROOT}/stat_pan/4.truvari_sv/0.sv_vcfgs"
out_dir="${PROJECT_ROOT}/stat_pan/4.truvari_sv/1.benchmark"
repeat_dir="${PROJECT_ROOT}/stat_pan/ref/ARS_UCD_v2.0.ref_repeat"
mkdir -p "$out_dir"

# Run truvari bench per sample + repeat region
for vcf in "$callset_dir"/*.sv.vcf.gz; do
    sample=$(basename "$vcf" .sv.vcf.gz)

    # General benchmark
    truvari bench \
        -b "$ref_vcf" \
        -c "$vcf" \
        -o "$out_dir/${sample}.all" \
        --passonly --sizemin 50 --sizemax 1000000 \
        -r 2000 -C 2000 --pctsim 0.7 --pctsize 0.7 --threads 4

    # Repeat-stratified
    for region in LINE SINE LTR Simple_repeat Satellite Low_complexity; do
        bed="${repeat_dir}/rm.${region}.bed"
        [ ! -f "$bed" ] && continue

        truvari bench \
            -b "$ref_vcf" \
            -c "$vcf" \
            -o "$out_dir/${sample}.${region}" \
            --includebed "$bed" \
            --passonly --sizemin 50 --sizemax 1000000 \
            -r 2000 -C 2000 --pctsim 0.7 --pctsize 0.7 --threads 4
    done
done
