# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${DATA_ROOT:?DATA_ROOT is unset - see config.sh.example at the repository root}"
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
: "${REF_DIR:?REF_DIR is unset - see config.sh.example at the repository root}"
# --------------------------

mkdir -p ${PROJECT_ROOT}/stat_pan/5.rtg_snv
cd ${PROJECT_ROOT}/stat_pan/5.rtg_snv
mkdir logs

chrs="$(echo X | awk '{for(i=1;i<=29;i++)printf i","; print}')"

#snv_hifi=${DATA_ROOT}/longReadCons/$sample/
#snv_hybird_mer=${PROJECT_ROOT}/pangenie_HiFi/4.deepv_all/deepv.merge7_hifi.vcf.gz
#snv_deepv_mer=${PROJECT_ROOT}/pangenie_HiFi/4.deepv_all/deepv.merge2.vcf.gz

snv_hifi=${DATA_ROOT}/longReadCons/${id}/merge_output.vcf.gz
snv_hybird=${PROJECT_ROOT}/pangenie_HiFi/${id/sample_/}/7.deepv_hybri/${id/sample_/}.vcf.gz
snv_deepv=${PROJECT_ROOT}/pangenie_HiFi/${id/sample_/}/4.deepv_data/${id/sample_/}.vcf.gz


${PROJECT_ROOT}/pangenie_HiFi/snps
${DATA_ROOT}/longReadCons/submit/*.merge_output.gvcf.gz

wget http://www.bio8.cs.hku.hk/clair3_trio/config/clair3.yml

ls ${DATA_ROOT}/longReadCons/submit/*.merge_output.gvcf.gz > input.list

glnexus_cli \
    --config clair3.yml \
    --dir GLnexus.DB \
    --mem-gbytes $((($nthreads-2)*8)) \
    --threads $(($nthreads-2)) \
    --list input.list > deepv.merge2.vcf.gz.bcf \
    2> deepv.merge2.vcf.gz.bcf.log
bcftools view --threads $nthreads deepv.merge2.vcf.gz.bcf -Oz > deepv.merge2.vcf.gz


mkdir -p 0.tmp 0.snp_vcfgs 1.rtg 2.regions
ls ${PROJECT_ROOT}/stat_pan/2.vcf_stats/2.ind/hol-pg2hic-2024-05-22.sample_*.anno_biallelic.filtered.vcf.gz | 
    grep -v "4439" | 
    while read id; do
        idd=$(basename $id | sed -e 's/hol-pg2hic-2024-05-22.//' -e 's/\.anno_biallelic.filtered.vcf.gz//')
        sbatch -A ${SLURM_ACCOUNT} -J $idd.panBed \
            -o logs/$idd.panBed.out \
            -c 1 --time "02-00:00:00" \
            --wrap="
bcftools sort $id -Oz -o 0.tmp/$idd.vcf.gz && tabix 0.tmp/$idd.vcf.gz
bcftools view --regions $chrs -v snps -e 'GT=\".\"' 0.tmp/$idd.vcf.gz | 
    bcftools sort -Oz -o 0.snp_vcfgs/$idd.pan.vcf.gz && tabix 0.snp_vcfgs/$idd.pan.vcf.gz
"
    done
    
pangenieMerged_vcf=${PROJECT_ROOT}/pangenie_HiFi/3.pan_all/hol-pg2hic-2024-05-22_graph_genotyping.merge-biallelic.filter.vcf.gz
cat ${PROJECT_ROOT}/stat_pan/hol.sample |
    while read id; do
        sbatch -A ${SLURM_ACCOUNT} -J $id.pangenie \
            -o logs/$id.pangenie.out \
            -c 1 --time "02-00:00:00" \
            --wrap="
bcftools view -s ${id/sample_/} --regions $chrs -v snps -e 'GT=\".\"' $pangenieMerged_vcf | 
    bcftools sort -Oz -o 0.snp_vcfgs/$id.pangenie-sr.vcf.gz && 
    tabix 0.snp_vcfgs/$id.pangenie-sr.vcf.gz
"
    done


###
vcf=${PROJECT_ROOT}/pangenie_HiFi/8.rna_all/holPub.filter_ind.g-rnasnps.vcf.gz
cat ${PROJECT_ROOT}/stat_pan/hol.sample |
    while read id; do
        sbatch -A ${SLURM_ACCOUNT} -J $id.rna \
            -o logs/$id.rna.out \
            -c 1 --time "02-00:00:00" \
            --wrap="
bcftools view -s ${id/sample_/} --regions $chrs -v snps $vcf  -e 'GT=\".\"' | 
    bcftools sort -Oz -o 0.snp_vcfgs/$id.snv_rna3.vcf.gz && 
    tabix 0.snp_vcfgs/$id.snv_rna3.vcf.gz
"
    done

vcf=${PROJECT_ROOT}/pangenie_HiFi/8.rna_all/rna.hifi19.filter_ind.rna-filter.vcf.gz
cat ${PROJECT_ROOT}/stat_pan/hol.sample |
    while read id; do
        sbatch -A ${SLURM_ACCOUNT} -J $id.rna \
            -o logs/$id.rna.out \
            -c 1 --time "01-00:00:00" \
            --wrap="
bcftools view -s ${id/sample_/} --regions $chrs -v snps $vcf -e 'GT=\".\"' | 
    bcftools sort -Oz -o 0.snp_vcfgs/$id.snv_rna4.vcf.gz && 
    tabix 0.snp_vcfgs/$id.snv_rna4.vcf.gz
"
    done

pangenieMerged_vcf2=${PROJECT_ROOT}/pangenie_HiFi/3.pan_all/jerHap-pg-2024-12-18_graph_genotyping.merge-biallelic.filter.vcf.gz
cat ${PROJECT_ROOT}/stat_pan/hol.sample |
    while read id; do
        idd=$id
        sbatch -A ${SLURM_ACCOUNT} -J $idd \
            -o $idd.out \
            -c 1 --time "01-00:00:00" \
            --wrap="
bcftools view -s ${id/sample_/} --regions $chrs -v snps $pangenieMerged_vcf2 -e 'GT=\".\"' | 
    bcftools sort -Oz -o 0.snp_vcfgs/$id.pangenie-holForJer.vcf.gz && 
    tabix 0.snp_vcfgs/$id.pangenie-holForJer.vcf.gz
"
    done

####################################################
# 1.6.sv.com.sh

ref_path=${REF_DIR}
ref_fa=$ref_path/ARS_UCD_v2.0.fa

tools1="pan pangenie-sr pangenie-holForJer snv_deepv snv_hybird snv_hifi snv_rna3 snv_rna4"
tools2="pan pangenie-sr pangenie-holForJer snv_deepv snv_hybird snv_hifi snv_rna3 snv_rna4"
cat ${PROJECT_ROOT}/stat_pan/hol.sample | 
    while read id; do
        for tool1 in $tools1; do
        for tool2 in $tools2; do
        if [ $tool1 == $tool2 ]; then continue; fi
            idd=$id
            i=1
            sbatch -A ${SLURM_ACCOUNT} -J $idd.$tool1.$tool2 \
                -o $idd.$tool1.$tool2.out \
                -c 4 --time "01-00:00:00" \
                --mem 200g \
                --wrap="
if [ -d 1.rtg/$id.$tool2.$tool1.rtg ]; then rm -rf 1.rtg/$id.$tool2.$tool1.rtg ; fi
rtg vcfeval -b 0.snp_vcfgs/$id.$tool2.vcf.gz -c 0.snp_vcfgs/$id.$tool1.vcf.gz -t ${ref_fa/fa/sdf} -o 1.rtg/$id.$tool2.$tool1.All.rtg -T 4
                "
            let i++
        done
    done
done

##


ref_path=${REF_DIR}
ref_fa=$ref_path/ARS_UCD_v2.0.fa
ref_rm=${REF_DIR}/ARS_UCD_v2.0.ref_repeat

#rms="sv ins del cpx None RM DNA LTR Low_complexity LINE srpRNA rRNA Unknown RNA RC scRNA SINE Satellite Simple_repeat tRNA snRNA"
rms="sv ins del cpx LTR Low_complexity LINE  SINE Satellite Simple_repeat"


cat ${PROJECT_ROOT}/stat_pan/hol.sample | 
    while read id; do
        for tool1 in $tools1; do
        for tool2 in $tools2; do
        if [ $tool1 == $tool2 ]; then continue; fi
            idd=$id
            i=1
            for rm in $rms; do
            sbatch -A ${SLURM_ACCOUNT} -J $idd.$tool1.$tool2 \
                -o $idd.$tool1.$tool2.out \
                -c 4 --time "02-00:00:00" \
                --mem 128g \
                --wrap="
if [ -d 1.rtg/$id.$tool2.$tool1.rtg ]; then rm -rf 1.rtg/$id.$tool2.$tool1.rtg ; fi
rtg vcfeval -b 0.snp_vcfgs/$id.$tool2.vcf.gz -c 0.snp_vcfgs/$id.$tool1.vcf.gz -t ${ref_fa/fa/sdf} -o 1.rtg/$id.$tool2.$tool1.$rm.rtg -T 4 -e $ref_rm/rm.$rm.bed 
                "
            let i++
        done
    done
done
done

echo "Sample Comp Base Threshold True-pos-baseline True-pos-call False-pos False-neg Precision Sensitivity F-measure" > all.summary.txt
i=0
ls 1.rtg/*/summary.txt | while read id ; do 
cat $id | awk '$1 == "None"{printf ID" ";for(i=1;i<=NF;i++){printf $i" "};printf "\n"}' \
    ID="$(echo $id |sed -e 's|1.rtg/||' -e 's|.rtg/summary.txt||' )"  >> all.summary.txt &
    if [ $i -eq 30 ]; then wait; i=0; fi
    let i++
    #echo $i $id
    done 



cat all.summary.txt | awk 'NR==1{print "Name Sample Comp Base Region Threshold True-pos-baseline True-pos-call False-pos False-neg Precision Sensitivity F-measure";next} {split($1,a,/\./);print $1,a[1],a[2],a[3],a[4],$2,$3,$4,$5,$6,$7,$8,$9}' > all.summary.table

scp -r ${CLUSTER_HOST}:${PROJECT_ROOT}/stat_pan/5.rtg_snv/all.summary.table .

# <local workstation path redacted - results were copied to a local workstation here>



##

zcat ${PROJECT_ROOT}/imputation/1.holPub_imp/holPub.pangenie-sv.filter.vcf.gz | awk '$1 !~ /#/{split($3,a,/[:_-]/);print a[1]"\t"a[2]"\t"a[3]"\t"a[4]}'  > $ref_rm/rm.sv.bed
cat $ref_rm/rm.sv.bed | awk '$4 == "INS"' > $ref_rm/rm.ins.bed
cat $ref_rm/rm.sv.bed | awk '$4 == "DEL"' > $ref_rm/rm.del.bed
cat $ref_rm/rm.sv.bed | awk '$4 == "COMPLEX"' > $ref_rm/rm.cpx.bed
