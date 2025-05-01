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

    ## mapping
    align_pbmm2 || minimap2_align
    bam_fn
    bam_calmd && sniffles_sv
    cuteSV_sv
    pbsv_sv
    svim_sv
    SVision_sv

    #module load miniconda3
    #source ${CONDA_BASE}/etc/profile.d/conda.sh

}

### mapping-based SV calling
index_pbmm2 () {
    if [[ ! -f "/path/to/file" ]]; then
        pbmm2 index $ref_fa ${ref_fa/fa/mmi} --preset "HiFi" -j $nthreads
    fi
}

function align_pbmm2() {
    conda activate $conda_envs/pbmm2

    pbmm2 align $ref_fa \
        1.clean_data/${sample}.hifi_reads.fastq.gz \
        --preset CCS \
        --rg "@RG\tID:$sample\tSM:$sample" \
        --log-level INFO \
        -j $nthreads \
        --sort | 
            samtools view -@ $nthreads -h -b -q 30 \
                -o 3.mapping/${sample}.map-pbmm2_ccs_mq30.bam

    samtools index -@ $nthreads 3.mapping/${sample}.map-pbmm2_ccs_mq30.bam
}

function minima2_index() {
    minimap2 -d ${ref_fa/fa/mmi} $ref_fa -t $nthreads
}

function minimap2_align() {
    minimap2 -ax map-hifi ${ref_fa/fa/mmi} -t $nthreads \
        1.clean_data/${sample}.hifi_reads.fastq.gz \
        -L --cs -R "@RG\tID:$sample\tSM:$sample" -H -Y --MD |
        samtools sort -@ $nthreads |
        samtools view -@ $nthreads -h -b -q 30 \
        -o 3.mapping/${sample}.map-mm2_ccs_mq30.bam
    samtools index -@ $nthreads 3.mapping/${sample}.map-mm2_ccs_mq30.bam
    ln -s 3.mapping/${sample}.map-mm2_ccs_mq30.bam 3.mapping/${sample}.map-pbmm2_ccs_mq30.bam
    samtools coverage -@ $nthreads 3.mapping/${sample}.map-pbmm2_ccs_mq30.bam > 3.mapping/${sample}.map-pbmm2_ccs_mq30.algin.coverage
}

function bam_fn() {
    conda activate ${CONDA_BASE}/envs/clair3
    mosdepth --threads $nthreads 3.mapping/${sample} 3.mapping/${sample}.map-pbmm2_ccs_mq30.bam
}

function bam_calmd() {
    samtools calmd -@ $nthreads -bS 3.mapping/${sample}.map-pbmm2_ccs_mq30.bam $ref_fa |
        samtools sort -@$nthreads > 3.mapping/${sample}.sv-sniffles.md.sort.bam
    samtools index -@ $nthreads 3.mapping/${sample}.sv-sniffles.md.sort.bam
}

function sniffles_sv() {
    conda activate $conda_envs/sniffles

    sniffles -t $nthreads \
        --reference $ref_fa \
        --sample-id ${sample} \
        -i 3.mapping/${sample}.sv-sniffles.md.sort.bam \
        -v 3.mapping/${sample}.sv-sniffles-v2.vcf \
        --minsupport 3 \
        --minsvlen 50 \
        --mapq 20 \
        --tandem-repeats $ref_tdr \
        --allow-overwrite
}

function cuteSV_sv() {
    conda activate $conda_envs/cuteSV
    cuteSV -t $nthreads \
        3.mapping/${sample}.map-pbmm2_ccs_mq30.bam \
        $ref_fa \
        3.mapping/${sample}.sv-cutesv.vcf \
        ./ \
        --max_cluster_bias_INS 1000 \
        --diff_ratio_merging_INS 0.9 \
        --max_cluster_bias_DEL 1000 \
        --diff_ratio_merging_DEL 0.5 \
        --min_size 50 \
        --min_support 3 \
        --genotype
}

function pbsv_sv() {
    conda activate $conda_envs/pbsv

    pbsv discover 3.mapping/${sample}.map-pbmm2_ccs_mq30.bam \
        3.mapping/${sample}.sv-pbsv.svsig.gz \
        --ccs --tandem-repeats $ref_tdr
    pbsv call -j $nthreads $ref_fa 3.mapping/${sample}.sv-pbsv.svsig.gz \
        3.mapping/${sample}.sv-pbsv-all.vcf --ccs --min-sv-length 50
    
    # mkdir -p 3.mapping/2-sv-pbsv
    # for i in $(samtools view -H 3.mapping/${sample}.map-pbmm2_ccs_mq30.bam | grep '^@SQ' | cut -f2 | cut -d':' -f2) ; do
    #     pbsv discover --region $i 3.mapping/${sample}.map-pbmm2_ccs_mq30.bam \
    #         3.mapping/2-sv-pbsv/${sample}.sv-pbsv.$i.svsig.gz \
    #         --hifi --tandem-repeats $ref_tdr -s ${sample}
    # done
    #   pbsv call -j $nthreads $ref_fa 3.mapping/ 3.mapping/2-sv-pbsv/${sample}.sv-pbsv.*.svsig.gz \
    #     3.mapping/${sample}.sv-pbsv-all.vcf --ccs --min-sv-length 50
}

function svim_sv() {
    conda activate $conda_envs/svim
    svim alignment --segment_gap_tolerance 10 --segment_overlap_tolerance 5 \
        --interspersed_duplications_as_insertions --tandem_duplications_as_insertions \
        --read_names --max_sv_size 1000000 \
        3.mapping/2-sv-svim/ 3.mapping/${sample}.map-pbmm2_ccs_mq30.bam $ref_fa
}

function SVision_sv() {
    conda activate $conda_envs/svisionenv-1.3.8
    model="${EXT_PROJECT_DIR}/uvm_mckay/software/SVision-1.3.8/svision_model"
    SVision -o 3.mapping/2-sv-svision/ \
        -b 3.mapping/${sample}.map-pbmm2_ccs_mq30.bam \
        -m ${model}/svision-cnn-model.ckpt \
        -g $ref_fa -n ${sample} -s 5 --graph --qname
}

function copy_sv_vcf () {
    dir=x.sv_vcf
    mkdir -p $dir
    bcftools reheader --sample <(echo Sample $sample) 2.assembly/$sample.svim-asm-bp/variants.vcf | bcftools sort  | bgzip -c > $dir/$sample.sv-svim-asm-bp.vcf.gz && tabix -p vcf $dir/$sample.sv-svim-asm-bp.vcf.gz
    bcftools reheader --sample <(echo Sample $sample) 2.assembly/$sample.svim-asm-diploid/variants.vcf | bcftools sort | bgzip -c > $dir/$sample.sv-svim-asm-diploid.vcf.gz && tabix -p vcf $dir/$sample.sv-svim-asm-diploid.vcf.gz
    bcftools reheader --sample <(echo Sample $sample) 3.mapping/2-sv-svim/variants.vcf | bcftools sort | bgzip -c > $dir/$sample.sv-svim.vcf.gz && tabix -p vcf $dir/$sample.sv-svim.vcf.gz
    bcftools reheader --sample <(echo 426 $sample) 3.mapping/2-sv-svision/*.svision.s5.graph.vcf | awk 'BEGIN{FS=OFS="\t"} {if($1 ~ "##INFO=<ID=GFA_L")$1=$1"\n##INFO=<ID=GFA_ID,Number=1,Type=String,Description=\"GFA_ID\">\n"; print}' | bcftools sort | bgzip -c > $dir/$sample.sv-svision.vcf.gz && tabix -p vcf $dir/$sample.sv-svision.vcf.gz

    bcftools sort 3.mapping/$sample.sv-pbsv-all.vcf | bgzip -c > $dir/$sample.sv-pbsv.vcf.gz && tabix -p vcf $dir/$sample.sv-pbsv.vcf.gz
    bcftools sort 3.mapping/$sample.sv-sniffles-v2.vcf | bgzip -c > $dir/$sample.sv-sniffles.vcf.gz && tabix -p vcf $dir/$sample.sv-sniffles.vcf.gz
    bcftools sort 3.mapping/$sample.sv-cutesv.vcf | bgzip -c > $dir/$sample.sv-cutesv.vcf.gz && tabix -p vcf $dir/$sample.sv-cutesv.vcf.gz

    cp 2.assembly/$sample.pav-bp/pav_$sample.vcf.gz $dir/$sample.pav_bp.vcf.gz && tabix -p vcf $dir/$sample.pav_bp.vcf.gz ||
    cp 2.assembly/$sample.pav-bp/$sample.vcf.gz $dir/$sample.pav_bp.vcf.gz && tabix -p vcf $dir/$sample.pav_bp.vcf.gz

    cp 2.assembly/$sample.pav_diploid/pav_$sample.vcf.gz $dir/$sample.pav_diploid.vcf.gz && tabix -p vcf $dir/$sample.pav_diploid.vcf.gz ||
    cp 2.assembly/$sample.pav_diploid/$sample.vcf.gz $dir/$sample.pav_diploid.vcf.gz && tabix -p vcf $dir/$sample.pav_diploid.vcf.gz
}

main
