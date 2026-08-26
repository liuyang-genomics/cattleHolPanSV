# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
: "${SOFTWARE_DIR:?SOFTWARE_DIR is unset - see config.sh.example at the repository root}"
# --------------------------


cd ${PROJECT_ROOT}/stat_pan

annotatePy="${SOFTWARE_DIR}/genotyping-pipelines/prepare-vcf-MC/workflow/scripts/annotate_vcf.py"

wd=2.vcf_stats/0.filter

sa cactus
ls ${PROJECT_ROOT}/minigraph-cactus/*/*.vcf.gz | grep -v 'raw' | grep -v 'filter' | grep -v 'bak' | grep -v bovinePan | grep -v hol |
    while read id; do
    sbatch -A ${SLURM_ACCOUNT} \
        --cpus-per-task=8 \
        --mem-per-cpu=8g \
        --wrap="
gunzip -c $id > ${id/.gz/} 
gunzip -c ${id/vcf/gfa} > ${id/vcf.gz/gfa}
python3 $annotatePy \
    -vcf ${id/.gz/} \
    -gfa ${id/vcf.gz/gfa} \
    -o ${id/vcf.gz/anno} &> ${id/vcf.gz/anno}.log
bgzip ${id/vcf.gz/anno}
bcftools sort -Oz -o ${id/vcf.gz/anno}.gz ${id/vcf.gz/anno}.gz
"
    done

ls ${PROJECT_ROOT}/minigraph-cactus/*/*.vcf.gz | grep -v 'raw' | grep -v 'filter' | grep -v 'bak' | grep -v bovinePan | grep -v hol |
    while read id; do
        idd=$(basename $id);
        idd=${idd/vcf/filtered.vcf};
        sbatch -A ${SLURM_ACCOUNT} \
        --cpus-per-task=2 \
        --mem-per-cpu=8g \
        --wrap="bcftools view --threads 2 -i 'F_MISSING<0.2' $id -Oz -o $wd/${idd}"
    done

####

wd=x.anno.stats # ========================================
mkdir -p $wd

ls ${PROJECT_ROOT}/minigraph-cactus/*/*.anno.log |  while read id; do 
    idd=$(basename $id);
    idd=${idd%%-*};
    tail -n2 $id | awk '{print idd,$3,$2}' idd=$idd; 
done > $wd/1.pan.raw.anno.stats


### vcf count
cat > vcf.svtype.sh <<'EOF'
bcftools query -f "%ID %INFO/ID\n" | 
    awk '{if($2 ~ /,|:/){split($2,a,/,|:/);for(i in a){gsub("-"," ",a[i]);print $1,a[i]}} else {gsub("-"," ",$2);print}}' | 
    awk '{print $1,$2,$3,$5,$4,$6}'
EOF

cat > count.svtype.sh <<'EOF'
awk '{gc[$5]++;gl[$5]+=$6;if($2 <= 29){ac[$5]++;al[$5]+=$6}} END {for(i in gc) print idd,i,gc[i],gl[i],ac[i],al[i]}' idd=$1
EOF


rm $wd/2.pan.raw.vcf.stats
rm $wd/2.pan.raw.vcf.stats.auto
ls ${PROJECT_ROOT}/minigraph-cactus/*/*.anno*vcf.gz |
    while read id; do
        idd=$(basename $id);
        ale=${idd#*.};
        ale=${ale%.*};
        idd=${idd%%-*};
#zcat $id | bash vcf.svtype.sh > ${id/vcf.gz/svtype}
cat ${id/vcf.gz/svtype} | bash count.svtype.sh $idd.$ale >> $wd/2.pan.raw.vcf.stats
 sbatch -A ${SLURM_ACCOUNT} \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap=""
    done


cat > list.svid.sh <<'EOF'
bcftools query -f "%ID %CHROM %POS %INFO/AT\n" | 
    awk '{if($4 ~ /,|:/){split($4,a,/,|:/);for(i in a){gsub("-"," ",a[i]);print a[i]}} else {print $4}}'
EOF


echo "${PROJECT_ROOT}/comp_pan/2.vcf_filter/0.filterF_M0.2/hol-pg2hic-2024-05-22.filtered.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/1.grp/allbovinePan-2024-07-03.hol.filtered.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/0.filterF_M0.2/jer-pg-2024-05-15.filtered.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/1.grp/allbovinePan-2024-07-03.jer.filtered.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/0.filterF_M0.2/bovinePan-2024-06-24.filtered.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/1.grp/allbovinePan-2024-07-03.pub.filtered.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/1.grp/bovinePanPhase-2024-07-03.pub.filtered.vcf.gz" |
    while read id; do
    idd=$(basename $id);
    iid=${idd/.filtered.vcf.gz/};
    idd=${idd%%.*}; 
    sbatch -A ${SLURM_ACCOUNT} \
        --cpus-per-task=1 \
        --mem-per-cpu=8g \
        --wrap="
    zcat $id | bash list.svid.sh > ${id/vcf.gz/svid}
    awk 'FNR==NR{a[\$1]++;next} {if(a[\$4]){print \$0}}' \
    ${id/vcf.gz/svid} \
    ${PROJECT_ROOT}/minigraph-cactus/$idd/$idd.anno.svtype \
    > ${id/vcf.gz/svtype}
    awk '{a[\$5]++} END {for(i in a) print idd,i,a[i]}' idd=$iid ${id/vcf.gz/svtype} >> $wd/grp.stats
    "
    done


echo "${PROJECT_ROOT}/comp_pan/2.vcf_filter/0.filterF_M0.2/hol-pg2hic-2024-05-22.filter-bi.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/1.grp/allbovinePan-2024-07-03.hol.filter-bi.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/0.filterF_M0.2/jer-pg-2024-05-15.filter-bi.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/1.grp/allbovinePan-2024-07-03.jer.filter-bi.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/0.filterF_M0.2/bovinePan-2024-06-24.filter-bi.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/1.grp/allbovinePan-2024-07-03.pub.filter-bi.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/1.grp/bovinePanPhase-2024-07-03.pub.filter-bi.vcf.gz" |
    while read id; do
    idd=$(basename $id);
    iid=${idd/.filter-bi.vcf.gz/};
    idd=${idd%%.*}; 
    sbatch -A ${SLURM_ACCOUNT} \
        --cpus-per-task=1 \
        --mem-per-cpu=8g \
        --wrap="
    zcat $id | bash tab.svtype.sh > ${id/vcf.gz/svid}
    awk 'FNR==NR{a[\$1]++;next} {if(a[\$4]){print \$0}}' \
    ${id/vcf.gz/svid} \
    ${PROJECT_ROOT}/minigraph-cactus/$idd/$idd.anno_biallelic.svtype \
    > ${id/vcf.gz/svtype}
    awk '{a[\$5]++} END {for(i in a) print idd,i,a[i]}' idd=$iid ${id/vcf.gz/svtype} >> $wd/grp-bi.stats
    "
    done



#### anno vcf count ind ####

wd=2.annvcf.stats
mkdir -p $wd

echo "${PROJECT_ROOT}/comp_pan/2.vcf_filter/2.ind/bovinePanPhase-2024-07-03_sample_4232.filtered.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/2.ind/bovinePanPhase-2024-07-03_Jersey_441.filtered.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/2.ind/hol-pg2hic-2024-05-22_sample_4232.filtered.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/2.ind/jer-pg-2024-05-15_Jersey_441.filtered.vcf.gz" |
    while read id; do
    idd=$(basename $id);
    iid=${idd/.filtered.vcf.gz/};
    idd=${idd%%_*}; 
    sbatch -A ${SLURM_ACCOUNT} \
        --cpus-per-task=1 \
        --mem-per-cpu=8g \
        --wrap="
    zcat $id | bash tab.svtype.sh > ${id/vcf.gz/svid}
    awk 'FNR==NR{a[\$1]++;next} {if(a[\$4]){print \$0}}' \
    ${id/.vcf.gz/svid} \
    ${PROJECT_ROOT}/minigraph-cactus/$idd/$idd.anno.svtype \
    > ${id/vcf.gz/svtype}
    awk '{a[\$5]++} END {for(i in a) print idd,i,a[i]}' idd=$iid ${id/vcf.gz/svtype} >> $wd/ind.stats
    "
    done


echo "${PROJECT_ROOT}/comp_pan/2.vcf_filter/2.ind/bovinePanPhase-2024-07-03_sample_4232.filter-bi.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/2.ind/bovinePanPhase-2024-07-03_Jersey_441.filter-bi.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/2.ind/hol-pg2hic-2024-05-22_sample_4232.filter-bi.vcf.gz
${PROJECT_ROOT}/comp_pan/2.vcf_filter/2.ind/jer-pg-2024-05-15_Jersey_441.filter-bi.vcf.gz" |
    while read id; do
    idd=$(basename $id);
    iid=${idd/.filter-bi.vcf.gz/};
    idd=${idd%%_*}; 
    sbatch -A ${SLURM_ACCOUNT} \
        --cpus-per-task=1 \
        --mem-per-cpu=8g \
        --wrap="
    zcat $id | bash tab.svtype.sh > ${id/vcf.gz/svid}
    awk 'FNR==NR{a[\$1]++;next} {if(a[\$4]){print \$0}}' \
    ${id/.vcf.gz/svid} \
    ${PROJECT_ROOT}/minigraph-cactus/$idd/$idd.anno_biallelic.svtype \
    > ${id/vcf.gz/svtype}
    awk '{a[\$5]++} END {for(i in a) print idd,i,a[i]}' idd=$iid ${id/vcf.gz/svtype} >> $wd/ind-bi.stats
    "
    done