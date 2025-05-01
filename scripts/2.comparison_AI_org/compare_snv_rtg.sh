#!/bin/bash
# compare_snv_rtg.sh
# Pairwise comparison of SNV callsets using RTG vcfeval across genome and repeat regions

set -euo pipefail

ref_sdf="${REF_DIR}/ARS_UCD_v2.0.sdf"
vcf_dir="${PROJECT_ROOT}/stat_pan/5.rtg_snv/0.snp_vcfgs"
out_dir="${PROJECT_ROOT}/stat_pan/5.rtg_snv/1.rtg"
repeat_dir="${PROJECT_ROOT}/stat_pan/ref/ARS_UCD_v2.0.ref_repeat"
mkdir -p "$out_dir"

# Methods to compare
tools="snv_hifi snv_deepv snv_hybird pan pangenie-sr pangenie-holForJer snv_rna3 snv_rna4"
regions="ALL LINE SINE LTR Simple_repeat Satellite Low_complexity"

samples="${PROJECT_ROOT}/stat_pan/hol.sample"

# Loop through sample pairs and regions
for sample in $(cat "$samples"); do
  for t1 in $tools; do
    for t2 in $tools; do
      [ "$t1" = "$t2" ] && continue

      base_vcf="${vcf_dir}/${sample}.${t2}.vcf.gz"
      call_vcf="${vcf_dir}/${sample}.${t1}.vcf.gz"
      [ ! -f "$base_vcf" ] || [ ! -f "$call_vcf" ] && continue

      for region in $regions; do
        region_opt=""
        out_suffix="${region}"
        if [ "$region" != "ALL" ]; then
          bed_file="${repeat_dir}/rm.${region}.bed"
          [ -f "$bed_file" ] || continue
          region_opt="-e $bed_file"
        fi

        out_path="${out_dir}/${sample}.${t2}.${t1}.${out_suffix}.rtg"
        rtg vcfeval -b "$base_vcf" -c "$call_vcf" -t "$ref_sdf" \
          -o "$out_path" -T 4 $region_opt
      done
    done
  done
done
