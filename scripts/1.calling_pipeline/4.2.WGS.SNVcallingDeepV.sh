#!/bin/bash

set -o nounset
set -o errexit

# Submit to sbatch ###################################################################################
sample=$1
nthreads=$2
work_path=$3

source $work_path/env.config

# ========= functions ==

module load apptainer
export SINGULARITY_CACHEDIR="${SCRATCH_DIR}/.singularity"
export APPTAINER_CACHEDIR="${SCRATCH_DIR}/.singularity"

work_path="${PROJECT_ROOT}/pangenie_holPub"

cactusVcf_name="jerHap-pg-2024-12-18"

sample_good="$work_path/fsample.list"
raw_path="$work_path/$sample/0.raw_data/"
clean_path="$work_path/$sample/1.clean_data/"
bam_path="$work_path/$sample/2.bam_data/"
pan_path="$work_path/$sample/3.pan_${cactusVcf_name%%-*}/"

deepv_data="$work_path/$sample/4.deepv_data/"
svtools_path="$work_path/$sample/5.svtools_data/"

bam_all="$work_path/2.bam_all/"
if [ ! -f $bam_all/bam.coverage.txt ]; then
    echo "auto_coverage auto_meandepth X_coverage X_meandepth Y_coverage Y_meandepth" > $bam_all/bam.coverage.txt
fi
pan_all="$work_path/3.pan_all/"
deepv_all="$work_path/4.deepv_all/"
svtools_all="$work_path/5.svtools_all/"

conda_svtools="${CONDA_BASE}/envs/svtools"
conda_pangenie="${CONDA_BASE}/envs/pangenie"

ref_path="${REF_DIR}/"
fa="ARS_UCD_v2.0.chr.fa"
ref_fa=$ref_path/$fa
ref_gap=$ref_path/Bos_taurus.ARS-UCD1.2.dna.toplevel.genomic_gaps.txt
ref_repeat="${EXT_PROJECT_DIR}/ruminant_t2t/existing_NCBI_references/Cattle/RM_GCF_002263795/GCF_002263795.3_ARS-UCD2.0_genomic.fna.out.gz"


#cactusVcf_name="hol-pg2hic-2024-05-22"
#cactusVcf_name="bovinePan-2024-08-13"
#cactusVcf_name="bovHol-2024-08-13"


pangenie_data="${PROJECT_ROOT}/pangenie_indexs/"
pangenieIndex="$pangenie_data/pangenieIndex_${fa}_${cactusVcf_name}/${fa}_pangenieIndex"
pangenieVcf="$pangenie_data/$cactusVcf_name.pangenome.vcf"
pangenieCallsetVcf="$pangenie_data/$cactusVcf_name.callset-biallelic.vcf"

pangenie_scripts="${SOFTWARE_DIR}/pangenie/pipelines/run-from-callset/scripts"
convertAllelicPy="$pangenie_scripts/convert-to-biallelic.py"
annotatePy="${SOFTWARE_DIR}/genotyping-pipelines/prepare-vcf-MC/workflow/scripts/annotate_vcf.py"

sratoolkit_bin=${SOFTWARE_DIR}/sratoolkit.3.1.0-centos_linux64/bin/fastq-dump

wags_sif="${SOFTWARE_DIR}/wags/wags.sif"

cnvnatorbin=100
cnvnator_ref="${REF_DIR}/0.refChr"

BIN_VERSION="1.6.1"
gpu_partition="${SLURM_GPU_PARTITION:-gpu}"


main(){
    #module load miniconda3

    mkdir -p ${work_path}/tmp
    export TMPDIR=$work_path/tmp

    cd $work_path

    bwa_alignment

    #module load miniconda
    mkdir -p ${deepv_data}/tmp
    export TMPDIR="${deepv_data}/tmp"

    module load apptainer
    export SINGULARITY_CACHEDIR="${SCRATCH_DIR}/.singularity"
    export APPTAINER_CACHEDIR="${SCRATCH_DIR}/.singularity"
    export SINGULARITY_TMPDIR=$TMPDIR
    export APPTAINER_TMPDIR=$TMPDIR

    # if [ -d "${deepv_data}/intermediate_results_dir" ]; then
    #     mv ${deepv_data}/intermediate_results_dir ${deepv_data}/${sample}.intermediate_results_dir
    # fi
    run_deepvariant_cpu
    postprocess_variants
    #rm_temp
    vcf_stats_report


    mkdir -p GLnexus.DB vcf_chunk sh_chunk

    #sumbi
    glnexus_cli_db


} 

function bwa_alignment () {
    # samtools faidx $ref_fa
    # bwa-mem2 index $ref_fa
    
    mkdir -p ${bam_path}
    #3 map to the genome
    bwa-mem2 mem -t $nthreads \
        -R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina\tLB:${sample}\tPU:${sample}" \
        ${ref_fa} \
        $clean_path/${sample}_1.clean.fq.gz \
        $clean_path/${sample}_2.clean.fq.gz |
        samtools sort -@ $nthreads -o ${bam_path}/${sample}.sorted.bam -
        # samtools sort -@ $nthreads -m 8G -o ${bam_path}/${sample}.sorted.bam - 
        
    bamtools stats -in ${bam_path}/${sample}.sorted.bam > ${bam_path}/${sample}.algin.stats
    samtools index ${bam_path}/${sample}.sorted.bam
    samtools coverage ${bam_path}/${sample}.sorted.bam > ${bam_path}/${sample}.algin.coverage
    
    awk 'NR > 1 {if($1 <= 29){a+=$6;b+=$7} else {ao[$1]=$6;bo[$1]=$7}} 
    END {printf sample" "a/29" "b/29" "; printf ao["X"]" "bo["X"]" "; printf ao["Y"]" "bo["Y"]"\n";}' sample=${sample} ${sample}/2.bam_data/${sample}.algin.coverage >> $bam_all/bam.coverage.txt
    # coverage meandepth X Y_coverage
    
    #samtools depth -a -q 20 -Q 20 -d 1000 -b $ref_gap $bam_path/${sample}.sorted.bam > $bam_path/${sample}.depth.txt
}

function bamMD () {
    gatk MarkDuplicates \
        --TMP_DIR $TMPDIR \
        -I ${bam_path}/${sample}.sorted.bam \
        -O ${bam_path}/${sample}.md.bam \
        -M ${bam_path}/${sample}.md.metrics.txt
    
    rm ${bam_path}/${sample}.sorted.bam &
    bamtools stats -in ${bam_path}/${sample}.md.bam > ${bam_path}/${sample}.algin.stats
    samtools index ${bam_path}/${sample}.md.bam
    samtools idxstats ${bam_path}/${sample}.md.bam | 
        awk '$1 == "1" || $1 == "Y"{a[$1]=$3} END{print ID,a["Y"]/a["1"]}' ID=$sample \
        > ${bam_path}/${sample}.chrY1_ratio.txt
    samtools coverage ${bam_path}/${sample}.md.bam > ${bam_path}/${sample}.algin.coverage

}

function run_deepvariant_cpu () {
    mkdir -p ${deepv_data}
    singularity run \
        -B /usr/lib/locale/:/usr/lib/locale/ \
        -B "${ref_path}":"${ref_path}" \
        -B "${bam_path}":"${bam_path}" \
        -B "${deepv_data}":"${deepv_data}" \
        -B "${deepv_data}/tmp":"$TMPDIR" \
        docker://google/deepvariant:1.6.1 \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref="$ref_fa" \
        --reads="${bam_path}/${sample}.sorted.bam" \
        --output_vcf="${deepv_data}"/${sample}.vcf.gz \
        --output_gvcf="${deepv_data}"/${sample}.g.vcf.gz \
        --intermediate_results_dir "${deepv_data}/${sample}.intermediate_results_dir" \
        --num_shards=$nthreads
#        
}

function postprocess_variants () {
    cd ${deepv_data}
    apptainer run \
    -B /usr/lib/locale/:/usr/lib/locale/ \
    -B "${ref_path}":"${ref_path}" \
    -B "${deepv_data}":"${deepv_data}" \
    -B "${deepv_data}/tmp":"$TMPDIR" \
    docker://google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/postprocess_variants \
    --ref "$ref_fa" \
    --infile "${sample}.intermediate_results_dir/call_variants_output.tfrecord.gz" \
    --outfile "${sample}.vcf.gz" \
    --cpus "$nthreads" \
    --gvcf_outfile "${sample}.g.vcf.gz" \
    --nonvariant_site_tfrecord_path "${sample}.intermediate_results_dir/gvcf.tfrecord@$nthreads.gz"
}

function rm_temp () {
    rm -rf ${deepv_data}/${sample}.intermediate_results_dir
}

function vcf_stats_report () {
    cd ${deepv_data}
    apptainer run \
        -B /usr/lib/locale/:/usr/lib/locale/ \
        -B "${ref_path}":"${ref_path}" \
        -B "${deepv_data}":"${deepv_data}" \
        -B "${deepv_data}/tmp":"$TMPDIR" \
        docker://google/deepvariant:1.6.1 \
        /opt/deepvariant/bin/vcf_stats_report \
        --input_vcf "${sample}.vcf.gz" \
        --outfile_base "${sample}"
}

glnexus_cli_db () {
    ls $work_path/*/4.deepv_data/*.g.vcf.gz | grep -f $sample_good > input.list

    rm -rf GLnexus.DB
    glnexus_cli \
        --config DeepVariant_unfiltered \
        --dir $deepv_all/GLnexus.DB \
        --mem-gbytes $((($nthreads-2)*8)) \
        --threads $(($nthreads-2)) \
        --list input.list > $deepv_all/deepv.merge2.vcf.gz.bcf \
        2> $deepv_all/deepv.merge2.vcf.gz.bcf.log
    bcftools view --threads $nthreads $deepv_all/deepv.merge2.vcf.gz.bcf -Oz > $deepv_all/deepv.merge2.vcf.gz
}


main


