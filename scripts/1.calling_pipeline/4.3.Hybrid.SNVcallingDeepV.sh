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
    #module load miniconda
    cd $work_path
    mkdir -p ${deepv_hybri}/tmp
    export TMPDIR="${deepv_hybri}/tmp"

    module load apptainer
    export SINGULARITY_CACHEDIR="${SCRATCH_DIR}/.singularity"
    export APPTAINER_CACHEDIR="${SCRATCH_DIR}/.singularity"
    export SINGULARITY_TMPDIR=$TMPDIR
    export APPTAINER_TMPDIR=$TMPDIR

    if [[ $sample -gt 1000 ]]; then
        hifi_sample="sample_"$sample
    else
       hifi_sample="jer_"$sample
    fi

    minimap2_align
    merge_LR_SR
    run_deepvariant_cpu
    #submit_next

    postprocess_variants
    #rm_temp
    vcf_stats_report

}

function minima2_index() {
    minimap2 -d $ref_fa.mmi $ref_fa
}
function minimap2_align() {
    mkdir -p $deepv_hybri
    minimap2 -ax map-hifi $ref_fa.mmi \
        -t $nthreads $hifi_loc/${hifi_sample}.hifi_reads.fastq.gz --MD |
        samtools sort -@ $nthreads - -o $deepv_hybri/${hifi_sample}.sorted.bam
    samtools index -@ $nthreads $deepv_hybri/${hifi_sample}.sorted.bam
}

function merge_LR_SR () {
    mkdir -p $deepv_hybri
    samtools merge -@ $nthreads \
        $deepv_hybri/$sample.hybri.bam \
        $deepv_hybri/${hifi_sample}.sorted.bam \
        ${bam_path}/${sample}.sorted.bam

    samtools index $deepv_hybri/$sample.hybri.bam
}

function run_deepvariant_cpu () {
    mkdir -p ${deepv_hybri}
    singularity run \
        -B /usr/lib/locale/:/usr/lib/locale/ \
        -B "${ref_path}":"${ref_path}" \
        -B "${bam_path}":"${bam_path}" \
        -B "${deepv_hybri}":"${deepv_hybri}" \
        -B "${deepv_hybri}/tmp":"$TMPDIR" \
        docker://google/deepvariant:1.6.1 \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type="HYBRID_PACBIO_ILLUMINA" \
        --ref="$ref_fa" \
        --reads="${deepv_hybri}"/${sample}.hybri.bam \
        --output_vcf="${deepv_hybri}"/${sample}.vcf.gz \
        --output_gvcf="${deepv_hybri}"/${sample}.g.vcf.gz \
        --intermediate_results_dir "${deepv_hybri}"/${sample}.intermediate_results_dir \
        --num_shards=$nthreads        
}


function postprocess_variants () {
    cd ${deepv_hybri}
    apptainer run \
    -B /usr/lib/locale/:/usr/lib/locale/ \
    -B "${ref_path}":"${ref_path}" \
    -B "${deepv_hybri}":"${deepv_hybri}" \
    -B "${deepv_hybri}/tmp":"$TMPDIR" \
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
    rm -rf ${deepv_hybri}/${sample}.intermediate_results_dir
}

function vcf_stats_report () {
    cd ${deepv_hybri}
    apptainer run \
        -B /usr/lib/locale/:/usr/lib/locale/ \
        -B "${ref_path}":"${ref_path}" \
        -B "${deepv_hybri}":"${deepv_hybri}" \
        -B "${deepv_hybri}/tmp":"$TMPDIR" \
        docker://google/deepvariant:1.6.1 \
        /opt/deepvariant/bin/vcf_stats_report \
        --input_vcf "${sample}.vcf.gz" \
        --outfile_base "${sample}"
}


glnexus_cli_db () {
    ls $work_path/*/7.deepv_hybri/*.g.vcf.gz | grep -f $sample_good > input7.list

    rm -rf GLnexus.DB
    glnexus_cli \
        --config DeepVariant_unfiltered \
        --dir $deepv_all/GLnexus.DB \
        --mem-gbytes $((($nthreads-2)*8)) \
        --threads $(($nthreads-2)) \
        --list input7.list |
        bcftools view - | bgzip -c \
        > $deepv_all/deepv.mergeh7.vcf.gz
}

main


