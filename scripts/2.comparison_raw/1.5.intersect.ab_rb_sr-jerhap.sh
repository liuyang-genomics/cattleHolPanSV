
mkdir -p ${PROJECT_ROOT}/stat_pan/3.intersect.ab_rb_sr
cd ${PROJECT_ROOT}/stat_pan/3.intersect.ab_rb_sr

ls ${PROJECT_ROOT}/stat_pan/2.vcf_stats/2.ind/jerHap-pg-2024-12-18.*.anno_biallelic.filtered.vcf.gz | 
    while read id; do
        idd=$(basename $id | sed -e 's/jerHap-pg-2024-12-18.//' -e 's/\.anno_biallelic.filtered.vcf.gz//')
        sbatch -A ${SLURM_ACCOUNT} -J $idd.panBed \
            -o $idd.panBed.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
zcat $id | bash panBed.sh > 0.panBed/$idd.bed
"
    done

wd="2.vcf_stats/2.ind"
mkdir -p $wd
cat ${PROJECT_ROOT}/stat_pan/jer.sample | 
    while read sample;do 
    id=2.vcf_stats/0.filter/jerPri-pg-2024-12-18.anno_biallelic.filtered.vcf.gz 
            idd=$(basename $id);
            idd_suffix=${idd#*.};
            idd=${idd%%.*};
            sbatch -A ${SLURM_ACCOUNT} -J $idd.ind -c 1 --wrap="
bcftools view --threads 1 -i 'F_MISSING<0.2' -c 1 -s $sample $id -Oz -o $wd/${idd}.${sample}.$idd_suffix
    "
        done

ls ${PROJECT_ROOT}/stat_pan/2.vcf_stats/2.ind/jerPri-pg-2024-12-18.jer_*.anno_biallelic.filtered.vcf.gz | 
    while read id; do
        idd=$(basename $id | sed -e 's/jerPri-pg-2024-12-18.//' -e 's/\.anno_biallelic.filtered.vcf.gz//')
        sbatch -A ${SLURM_ACCOUNT} -J $idd.mcPri \
            -o $idd.mcPri.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
zcat $id | bash panBed.sh > 1.ab_rb_srBed/$idd.mcPri.bed
"
    done


# 1.cute
ls ${PANEL_DIR}/jer_*/4.sv_vcf/jer_*.vcf.gz |
    while read id; do
        idd=$(basename $id | sed -e 's/.vcf.gz//')
        sbatch -A ${SLURM_ACCOUNT} -J $idd \
            -o $idd.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
zcat $id | bash ab_rbBed.sh > 1.ab_rb_srBed/$idd.bed
"
    done
#svision $7 is 0xxx


# 2.srSVtools
svtoolsMerged_vcf=${PROJECT_ROOT}/pangenie_HiFi/5.svtools_all/output.ls.filter.vcf.gz
bcftools view -h $svtoolsMerged_vcf > hdr.txt
sed -i 's|##bcftools_viewVersion|##FORMAT=<ID=CN,Number=A,Type=Float,Description="Copy number">\n##bcftools_viewVersion|' hdr.txt
cat ${PROJECT_ROOT}/stat_pan/jer.sample |
    while read id; do
        sbatch -A ${SLURM_ACCOUNT} -J $idd \
            -o $idd.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
bcftools reheader -h hdr.txt $svtoolsMerged_vcf | bcftools view -i 'F_MISSING<0.2' | 
    bcftools view -s ${id/jer_/} -c 1 --force-samples | bash ab_rbBed.sh > 1.ab_rb_srBed/$id.svtools-sr.bed
"
    done

# lumpy
cat ${PROJECT_ROOT}/stat_pan/jer.sample |
    while read id; do 
        idd=$id
        sbatch -A ${SLURM_ACCOUNT} -J $idd \
            -o $idd.lumpy-sr.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
    zcat ${PROJECT_ROOT}/pangenie_HiFi/*/5.svtools_data/2.lumpy_vcf/${id/jer_/}.part1-smoove.genotyped.vcf.gz | bash ab_rbBed.sh > 1.ab_rb_srBed/$id.lumpy-sr.bed
"
    done


# cnvnator 
ls ${PROJECT_ROOT}/pangenie_HiFi/*/5.svtools_data/5.cnvnator_rawcnv/*.cnvnator.bed |
    while read id; do
        idd=jer_$(basename $id | sed -e 's/\..*//')
        sbatch -A ${SLURM_ACCOUNT} -J $idd \
            -o $idd.cnvnator-sr.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
    cat $id | awk -f cnvnator.bed.covert.awk > 1.ab_rb_srBed/$idd.cnvnator-sr.bed
"
    done



pangenieMerged_vcf=${PROJECT_ROOT}/pangenie_HiFi/3.pan_all/jerHap-pg-2024-12-18_graph_genotyping.merge-biallelic.filter.vcf.gz
cat ${PROJECT_ROOT}/stat_pan/jer.sample |
    while read id; do
        idd=$id
        sbatch -A ${SLURM_ACCOUNT} -J $idd \
            -o $idd.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
bcftools view -s ${id/jer_/} -c 1 $pangenieMerged_vcf  | 
    grep -v '#' | awk -f pangenie-sr_bed.awk > 1.ab_rb_srBed/$idd.pangenie-sr.bed
"
    done


####

pangenieMerged_vcf2=${PROJECT_ROOT}/pangenie_HiFi/3.pan_all/hol-pg2hic-2024-05-22_graph_genotyping.merge-biallelic.filter.vcf.gz
cat ${PROJECT_ROOT}/stat_pan/jer.sample |
    while read id; do
        idd=$id
        sbatch -A ${SLURM_ACCOUNT} -J $idd \
            -o $idd.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
bcftools view -s ${id/jer_/} -c 1 $pangenieMerged_vcf2  | 
    grep -v '#' | awk -f pangenie-sr_bed.awk > 1.ab_rb_srBed/$idd.pangenie-jerForHol.bed
"
    done

#################################################################################
#######################################
#### intersect for pan sv
mkdir -p 5.sv.intersect_both
mkdir -p 6.stats-jer
mkdir -p 7.share_stats-jer

tools="pav_bp pav_diploid sv-cutesv sv-pbsv sv-sniffles sv-svim sv-svim-asm-bp sv-svim-asm-diploid sv-svision svtools-sr lumpy-sr cnvnator-sr pangenie-sr pangenie-jerForHol mcPri"
cat ${PROJECT_ROOT}/stat_pan/jer.sample | 
    while read id; do
        for tool in $tools; do
            idd=$id
            sbatch -A ${SLURM_ACCOUNT} -J $idd.jerBed \
                -o $idd.$tool.out \
                --cpus-per-task=1 \
                --mem-per-cpu=8g \
                --wrap="
bedtools intersect \
    -a 0.panBed/$id.bed \
    -b 1.ab_rb_srBed/$id.$tool.bed \
    -f 0.90 -r -wo > 5.sv.intersect_both/$id.${tool}_pan.bed
cat  0.panBed/$id.bed | awk -f pan_stats.awk > 6.stats-jer/$id.pan.stats
cat 1.ab_rb_srBed/$id.$tool.bed | awk -f pan_stats.awk > 6.stats-jer/$id.$tool.stats
awk -f pan_share_stats.awk 5.sv.intersect_both/$id.${tool}_pan.bed > 7.share_stats-jer/$id.$tool.sv.share.stats
                "
        done
    done


scp -r ${CLUSTER_HOST}:${PROJECT_ROOT}/stat_pan/3.intersect.ab_rb_sr/6.stats-jer
scp -r ${CLUSTER_HOST}:${PROJECT_ROOT}/stat_pan/3.intersect.ab_rb_sr/7.share_stats-jer 


C:\Users\${USER}\OneDrive - University of Maryland\Data\2024-02-07.cattleLR-SR-GWAS\2.analyses\1.5.intersect

tools="pav_bp pav_diploid sv-cutesv sv-pbsv sv-sniffles sv-svim sv-svim-asm-bp sv-svim-asm-diploid sv-svision svtools-sr lumpy-sr cnvnator-sr pangenie-sr pangenie-jerForHol mcPri"
cat ${PROJECT_ROOT}/stat_pan/jer.sample | 
    while read id; do
        for tool in $tools; do
            idd=$id
cat  0.panBed/$id.bed | awk -f pan_stats.awk > 6.stats-jer/$id.pan.stats
cat 1.ab_rb_srBed/$id.$tool.bed | awk -f pan_stats.awk > 6.stats-jer/$id.$tool.stats
awk -f pan_share_stats.awk 5.sv.intersect_both/$id.${tool}_pan.bed > 7.share_stats-jer/$id.$tool.sv.share.stats
        done
    done
