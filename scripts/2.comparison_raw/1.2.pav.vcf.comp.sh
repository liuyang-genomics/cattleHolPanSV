# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
# --------------------------

cd ${PROJECT_ROOT}/comp_pan

mkdir -p 3.pan_vcf_compa

zcat 2.vcf_filter/1.grp/allbovinePan-2024-07-03.hol.filter-bi.vcf.gz | grep -v "##" | cut -f 1-5 > 3.pan_vcf_compa/allbovinePan.hol.filter-bi.bed

zcat 2.vcf_filter/0.filterF_M0.2/hol-pg2hic-2024-05-22.filter-bi.vcf.gz | grep -v "##" | cut -f 1-5 > 3.pan_vcf_compa/hol.filter-bi.bed

zcat 2.vcf_filter/1.grp/allbovinePan-2024-07-03.jer.filtered.vcf.gz | grep -v "##" | cut -f 1-5 > 3.pan_vcf_compa/allbovinePan.jer.filter-bi.bed

zcat 2.vcf_filter/0.filterF_M0.2/jer-pg-2024-05-15.filter-bi.vcf.gz | grep -v "##" | cut -f 1-5 > 3.pan_vcf_compa/jer.filter-bi.bed
