# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${CONDA_BASE:?CONDA_BASE is unset - see config.sh.example at the repository root}"
: "${DATA_DIR:?DATA_DIR is unset - see config.sh.example at the repository root}"
: "${PANEL_DIR:?PANEL_DIR is unset - see config.sh.example at the repository root}"
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
: "${REF_DIR:?REF_DIR is unset - see config.sh.example at the repository root}"
: "${SCRATCH_DIR:?SCRATCH_DIR is unset - see config.sh.example at the repository root}"
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

    # merge_fq or link_hifi_rawdata
    if [[ ! -f $hifi_raw/${sample}${hifiSuffix} ]]; then
        merge_fq
    fi

    link_hifi_rawdata

    fastp_hifi
    # markdup
    bioawk_stat
    seqkit_stat
    # rm_intermediate_reads

    ## assembly
    if [[ "$hic" == "hic" ]]; then
        haps="hic hic.hap1 hic.hap2"
        fastp_hic 
        # link_hic_clean
        asm_hic
    else
        haps="bp bp.hap1 bp.hap2"
        asm
    fi

    echo "line,length,type,coverage" > 2.assembly/${sample}.coverage-hifiasm.csv

    # minima2_index
    
    for hap in $haps; do
        echo ${sample}.${hap}.p_ctg.fa
        gfa2fa
        # link_fa
        mapFa_minimap2
        busco_stat
        quast_stat
        inspector_stat
        coverage_hifiasm
        btk_snail
    done

}


function merge_fq () {
    ## provide a list of samples to merge
    zcat $hifi_raw/${sample}/$hifiSuffix | gzip -c > 0.raw_data/${sample}${hifiSuffix}
}

function link_hifi_rawdata () {
    ln -s $hifi_raw/${sample}${hifiSuffix} 0.raw_data/${sample}${hifiSuffix}
    ln -s $hifi_raw/${sample}${hifiSuffix} 1.clean_data/${sample}${hifiSuffix}
}

function fastp_hifi () {
    conda activate base
    fastp -i 0.raw_data/${sample}${hifiSuffix} \
        -l 1000 \
        -w $nthreads \
        -j 1.clean_data/${sample}.hifi_reads_qc.json \
        -h 1.clean_data/${sample}.hifi_reads_qc.html \
        -o 1.clean_data/${sample}.hifi_reads_qc.fastq.gz \
        2> 1.clean_data/${sample}.hifi_reads_qc.log
}

function markdup () {
    conda activate $conda_envs/pbmarkdup

    pbmarkdup \
        1.clean_data/${sample}.hifi_reads_qc.fastq.gz \
        1.clean_data/${sample}.hifi_reads_qc_rmdup.fastq.gz \
        -j $nthreads -r \
        2> 1.clean_data/${sample}.hifi_reads_qc_rmdup.log
}

function bioawk_stat () {
    conda activate $conda_envs/bioawk
    bioawk -c fastx '{print "PacBio_HiFi," length($seq)}' \
        1.clean_data/${sample}${hifiSuffix} \
        > 1.clean_data/${sample}.hifi_reads.csv
    # bioawk -c fastx '{print "PacBio_HiFi," length($seq)}' \
    #     1.clean_data/${sample}.hifi_reads_qc.fastq.gz \
    #     > 1.clean_data/${sample}.hifi_reads_qc.csv
    # bioawk -c fastx '{print "PacBio_HiFi," length($seq)}' \
    #     1.clean_data/${sample}.hifi_reads_qc_rmdup.fastq.gz \
    #     > 1.clean_data/${sample}.hifi_reads_qc_rmdup.csv        
}

function seqkit_stat () {
    module load seqkit
    seqkit stats 1.clean_data/${sample}${hifiSuffix} -a > 1.clean_data/${sample}.hifi_reads.txt
    seqkit stats 1.clean_data/${sample}.hifi_reads_qc.fastq.gz -a > 1.clean_data/${sample}.hifi_reads_qc.txt
   # seqkit stats 1.clean_data/${sample}.hifi_reads_qc_rmdup.fastq.gz -a > 1.clean_data/${sample}.hifi_reads_qc_rmdup.txt
}

function rm_intermediate_reads () {
    rm 0.raw_data/${sample}${hifiSuffix}\
        1.clean_data/${sample}${hifiSuffix} \
        1.clean_data/${sample}.hifi_reads_qc_rmdup.fastq.gz
}

#### hifi_reads assembly based sv calling

function fastp_hic () {
    fastp -w $nthreads \
        -i ${hic_raw}/${sample}${hicSuffix1} \
        -I ${hic_raw}/${sample}${hicSuffix2} \
        -o 1.clean_data/${sample}.hic.R1.fq.gz \
        -O 1.clean_data/${sample}.hic.R2.fq.gz \
        -h 1.clean_data/${sample}.hic_reads.fastq.html \
        -j 1.clean_data/${sample}.hic_reads.fastq.json
}

function link_hic_clean () {
    ln -s ${hic_raw}/${sample}${hicSuffix1} 1.clean_data/${sample}.hic.R1.fq.gz
    ln -s ${hic_raw}/${sample}${hicSuffix2} 1.clean_data/${sample}.hic.R2.fq.gz
}

function asm() {
    conda activate $conda_envs/hifiasm
    hifiasm -t $nthreads -r 3 1.clean_data/${sample}.hifi_reads.fastq.gz -o 2.assembly/${sample}
}

function asm_hic() {
    conda activate $conda_envs/hifiasm
    hifiasm -o 2.assembly/${sample} -t $nthreads -r 3 \
    1.clean_data/${sample}.hifi_reads.fastq.gz \
    --h1 1.clean_data/${sample}.hic.R1.fq.gz \
    --h2 1.clean_data/${sample}.hic.R2.fq.gz
}

function gfa2fa() {
    conda activate $conda_envs/gfatools
    gfatools gfa2fa 2.assembly/${sample}.${hap}.p_ctg.gfa >2.assembly/${sample}.${hap}.p_ctg.fa
}

function minima2_index() {
    minimap2 -d ${ref_fa/fa/mmi} $ref_fa -t $nthreads
}

function mapFa_minimap2() {
    conda activate $conda_envs/minimap2
    minimap2 -t $nthreads -x asm20 -m 10000 \
        -z 10000,50 -r 50000 --end-bonus=100 \
        --secondary=no -O 5,56 -E 4,1 -B 5 -a --eqx \
        -Y $ref_fa 2.assembly/${sample}.${hap}.p_ctg.fa |
            samtools sort -@$nthreads -o  2.assembly/${sample}.hifiasm-${hap}.bam
    samtools index 2.assembly/${sample}.hifiasm-${hap}.bam
}

function busco_stat() {
    conda activate $conda_envs/busco
    busco -i 2.assembly/${sample}.${hap}.p_ctg.fa -f \
        -c $nthreads -o 2.assembly/${sample}.busco-hifiasm-${hap} -m genome -l $odb_path/$odb_name --offline
}

function quast_stat() {
    conda activate $conda_envs/quast
    quast -t $nthreads 2.assembly/${sample}.${hap}.p_ctg.fa -r $ref_fa \
        -g $ref_gff3 -o 2.assembly/${sample}.quast-hifiasm-${hap}
}

function inspector_stat() {
    conda activate $conda_envs/inspector
    inspector.py -t $nthreads \
        -c 2.assembly/${sample}.${hap}.p_ctg.fa \
        -r 0.raw_data/${sample}${hifiSuffix} \
        -o 2.assembly/${sample}.inspector-hifiasm-${hap} --datatype hifi
}

function coverage_hifiasm() {
    conda activate $conda_envs/bioawk

    REF1=2.assembly/${sample}.${hap}.p_ctg.fa
    TYPE1=contig                                                      # Specify your genome assembly type, such as contig, scaffold, chromosome, etc.
    LEN1=$(bioawk -c fastx '{sum+=length($seq)}END{print sum}' $REF1) # Size of assembled genome
    #echo "line,length,type,coverage" > 1.4-coverage-hifiasm.csv
    STRAIN1=Hifi
    cat $REF1 | bioawk -c fastx -v line="$STRAIN1" '{print line","length($seq)","length($seq)}' | 
        sort -k3rV -t "," | 
        awk -F "," -v len="$LEN1" -v type="$TYPE1" 'OFS=","{print $1,$2,type,(sum+0)/len; sum+=$3 }' \
        >> 2.assembly/${sample}.coverage-hifiasm.csv
}

function btk_snail (){
    conda activate $conda_envs/btk

    blobtools create \
        --fasta 2.assembly/${sample}.${hap}.p_ctg.fa \
        --busco 2.assembly/${sample}.busco-hifiasm-${hap}/run_mammalia_odb10/full_table.tsv \
        2.assembly/${sample}.${hap}.snail
    
    blobtools view \
    --plot \
    --view snail \
    2.assembly/${sample}.${hap}.snail
}


main
