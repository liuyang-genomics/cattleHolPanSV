# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
# --------------------------

cat hol-pg2hic-2024-05-22.pangenome.vcf | grep -v "##" | wc -l
15638895
cat hol-pg2hic-2024-05-22.callset.vcf | grep -v "##" | wc -l
15645048
cat hol-pg2hic-2024-05-22.vcf  | grep -v "##" | wc -l
16234934
cat hol-pg2hic-2024-05-22.filtered_ids.vcf | grep -v "##" | wc -l
15645048

    python3 $annotatePy \
        -vcf $cactusVcf_name.callset.vcf \
        -gfa $cactusVcf_name.gfa \
        -o $cactusVcf_name.callset_annotation


cd ${PROJECT_ROOT}/comp_pan

ls ../minigraph-cactus/*/*vcf.gz | grep -v raw | grep -v filter | xargs -I {} ln -s {} .

ls *gz | xargs -I {} sbatch -A ${SLURM_ACCOUNT} --cpus-per-task=2 --mem-per-cpu=16g --wrap="bcftools view -i 'F_MISSING<0.2' {} -Oz -o {}.filter.gz"


ls $work_path/*/3.pan_data/*_genotyping.vcf.gz |
    while read id; do
    sbatch -A ${SLURM_ACCOUNT} \
    --cpus-per-task=1 \
    --mem-per-cpu=8g \
    --wrap="wc -l ${id} >> 3.pan_all/hol_genotyping.count"
    done

ls $work_path/*/3.pan_data_jer/*_genotyping.vcf.gz |
    while read id; do
    sbatch -A ${SLURM_ACCOUNT} \
    --cpus-per-task=1 \
    --mem-per-cpu=8g \
    --wrap="wc -l ${id} >> 3.pan_all/jer_genotyping.count"
    done


ls *.filter.gz |
    while read id; do
    sbatch -A ${SLURM_ACCOUNT} \
    --cpus-per-task=8 \
    --mem-per-cpu=8g \
    --wrap="bcftools norm --threads 8 -m- $id -Oz -o ${id/filter.gz/filter-bi.gz}"
    done


bcftools view -i 'F_MISSING<0.2' $1 -Oz -o $1.filter.gz



bcftools view -Oz -c1 -s sample1 joined.vcf.gz > sample1.vcf.gz
bcftools view -Oz -c1 -s sample2 joined.vcf.gz > sample2.vcf.gz
and then compare them

bcftools stats sample1.vcf.gz sample2.vcf.gz > joined.stats.txt
	


bcftools norm --threads 8 -m- All_graph_genotyping.merge.vcf.gz -Oz -o All_graph_genotyping.merge-biallelic.vcf.gz



sbatch -A ${SLURM_ACCOUNT} --cpus-per-task=2 --mem-per-cpu=16g --wrap="bcftools merge --threads=16 -R $pangenieVcf -l sa.list -Oz -o jer8_graph_genotyping.merge.vcf.gz"


ls *filtered*gz | 
while read id; do
    sbatch -A ${SLURM_ACCOUNT} \
    --cpus-per-task=4 \
    --mem-per-cpu=8g \
    --wrap="bctools stats --threads 4 "
    done