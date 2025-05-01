#cat > sv-panel.config << 'EOF'

#!/bin/bash

#set -o nounset
#set -o errexit

sample=$sample
nthreads=$SLURM_CPUS_PER_TASK

#481
#sample_4439

#hic="hic"
hic=""

work_dir=${PANEL_DIR}
conda_envs=${EXT_PROJECT_DIR}/uvm_mckay/software/miniconda3/envs

ref_path=${REF_DIR}
ref_fa=$ref_path/ARS_UCD_v2.0.fa
ref_gff3=$ref_path/ARS_UCD_v2.0.gff
ref_tdr=$ref_path/ARS-UCD2.0_trf2.bed
# https://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/bigZips/bosTau9.chromAlias.txt

hic_raw="${DATA_DIR}/clean_Hi-C_8samples/"
hicSuffix1="_1.clean.fq.gz"
hicSuffix2="_2.clean.fq.gz"

hifi_raw='${PROJECT_ROOT}/0.reads_HiFi/'
hifiSuffix=".hifi_reads.fastq.gz"
#assembly_path=${PROJECT_ROOT}/0.hifiasm

export SINGULARITY_CACHEDIR="${SCRATCH_DIR}/.singularity"
export APPTAINER_CACHEDIR="${SCRATCH_DIR}/.singularity"
export SINGULARITY_TMPDIR=$TMPDIR
export APPTAINER_TMPDIR=$TMPDIR

odb_path=${EXT_PROJECT_DIR}/uvm_mckay/bovine_genome
odb_name=mammalia_odb10
model_path='${CONDA_BASE}/envs/clair3/bin/models/hifi'
platform='hifi' # 'ont' or 'pacbio'
#### hifi_reads

singularity_cache="${SCRATCH_DIR}"

export SINGULARITY_CACHEDIR="$singularity_cache/.singularity"
#export SINGULARITY_TMPDIR=""
export APPTAINER_CACHEDIR="$singularity_cache/.singularity"

# singularity pull pav_latest.sif library://becklab/pav/pav:latest

# module load miniconda3


source ${CONDA_BASE}/etc/profile.d/conda.sh

module load samtools

mkdir -p $work_dir/${sample}/
cd $work_dir/${sample}/

mkdir -p 0.raw_data 1.clean_data 2.assembly 3.mapping x.sv_vcf

function main() {

    # module load miniconda3
    source ${CONDA_BASE}/etc/profile.d/conda.sh

    module load samtools

    mkdir -p $work_dir/${sample}/
    cd $work_dir/${sample}/

    mkdir -p 0.raw_data 1.clean_data 2.assembly 3.mapping 

    run_clair3
   
}


function run_clair3() {
    conda activate ${CONDA_BASE}/envs/clair3
    run_clair3.sh \
        --bam_fn=3.mapping/${sample}.map-pbmm2_ccs_mq30.bam \
        --ref_fn=$ref_fa \
        --threads=$nthreads \
        --platform=${platform} \
        --sample_name=${sample} \
        --model_path=${model_path} \
        --remove_intermediate_dir \
        --include_all_ctgs \
        --gvcf \
        --output=3.mapping/${sample}_clair3

    # --ctg_name=${chr} \
    # --include_all_ctgs
}


main