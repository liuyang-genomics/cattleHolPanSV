# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
: "${DATA_DIR:?DATA_DIR is unset - see config.sh.example at the repository root}"
: "${PANEL_DIR:?PANEL_DIR is unset - see config.sh.example at the repository root}"
: "${SCRATCH_DIR:?SCRATCH_DIR is unset - see config.sh.example at the repository root}"
: "${REF_DIR:?REF_DIR is unset - see config.sh.example at the repository root}"
: "${REF_FA:?REF_FA is unset - see config.sh.example at the repository root}"
: "${REF_GFF:?REF_GFF is unset - see config.sh.example at the repository root}"
: "${REF_TRF_BED:?REF_TRF_BED is unset - see config.sh.example at the repository root}"
: "${EXT_PROJECT_DIR:?EXT_PROJECT_DIR is unset - see config.sh.example at the repository root}"
: "${CONDA_BASE:?CONDA_BASE is unset - see config.sh.example at the repository root}"
: "${BUSCO_ODB_DIR:?BUSCO_ODB_DIR is unset - see config.sh.example at the repository root}"
# --------------------------

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
ref_fa=${REF_FA}
ref_gff3=${REF_GFF}
ref_tdr=${REF_TRF_BED}
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

odb_path=${BUSCO_ODB_DIR}
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

    #hap="hic"
    hap="bp"
    svim_asm_haploid
    pav_haploid

    svim_asm_diploid
    pav_diploid
    #jasmine_merge_svim-asm
}

function svim_asm_haploid () {
    conda activate $conda_envs/svim-asm
    svim-asm haploid --min_sv_size 50 \
        --max_sv_size 200000 \
        --tandem_duplications_as_insertions \
        --interspersed_duplications_as_insertions \
        2.assembly/${sample}.svim-asm_haploid/ \
        2.assembly/${sample}.hifiasm-${hic}.bam \
        $ref_fa
    # svim-asm diploid
}

function svim_asm_diploid () {
    conda activate $conda_envs/svim-asm
    svim-asm diploid --min_sv_size 50 \
        --max_sv_size 200000 \
        --tandem_duplications_as_insertions \
        --interspersed_duplications_as_insertions \
        2.assembly/${sample}.svim-asm_diploid/ \
        2.assembly/${sample}.hifiasm-${hic}.hap1.bam \
        2.assembly/${sample}.hifiasm-${hic}.hap1.bam \
        $ref_fa
}

function pav_haploid() {
    mkdir -p 2.assembly/${sample}.pav_haploid
    cd 2.assembly/${sample}.pav_haploid

    echo "{
    \"reference\": $ref_fa.gz
}" > config.json 

    cat >assemblies.tsv << EOF2 
NAME HAP1 HAP2
$sample $work_dir/${sample}/2.assembly/${sample}.${hic}.p_ctg.fa
EOF2

    sed -i 's/ /\t/g' assemblies.tsv
    conda activate base
    singularity run \
        --bind "$singularity_cache:$singularity_cache" \
        --bind "$(pwd):$(pwd)" \
        --bind "$ref_path:$ref_path" \
        --bind "$work_dir:$work_dir" \
        $work_dir/pav_latest.sif -c $nthreads

    cd $work_dir/${sample}/
}

function pav_diploid() {
    mkdir -p 2.assembly/${sample}.pav_diploid
    cd 2.assembly/${sample}.pav_diploid

    echo "{
    \"reference\": $ref_fa.gz
}" > config.json 

    cat >assemblies.tsv << EOF2 
NAME HAP1 HAP2
$sample $work_dir/${sample}/2.assembly/${sample}.${hic}.hap1.p_ctg.fa $work_dir/${sample}/2.assembly/${sample}.${hic}.hap2.p_ctg.fa
EOF2

    sed -i 's/ /\t/g' assemblies.tsv
    conda activate base
    singularity run \
        --bind "$singularity_cache:$singularity_cache" \
        --bind "$(pwd):$(pwd)" \
        --bind "$ref_path:$ref_path" \
        --bind "$work_dir:$work_dir" \
        $work_dir/pav_latest.sif -c $nthreads

    cd $work_dir/${sample}/
}

main
