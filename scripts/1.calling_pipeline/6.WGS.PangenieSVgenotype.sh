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

    mkdir -p ${work_path}

   # module load miniconda
   # source activate $conda_pangenie
    export PATH=${SOFTWARE_DIR}/pangenie/build/src:$PATH

   # cactusVcf2PanGenie
    PanGenie_index # run for once
    
    PanGenie_genotyping # run for each sample
    PanGenie_genotyping_biallelic # run for each sample

    PanGenie_genotyping_merge3 # run for once for all samples
    
    #submit_next
    # https://github.com/eblerjana/pangenie/blob/master/pipelines/run-from-callset/Snakefile


}

function cactusVcf2PanGenie () {
    # vcfbub -l 0 -r 100000 --input <your-vcf-file> > pangenie-ready.vcf
    # already done by minigraph-cactus with output $cactusVcf_name.vcf.gz no raw tagged

    gunzip -c $panVcf > $cactusVcf_name.vcf

    cat $cactusVcf_name.vcf |
        python3 $pangenie_scripts/prepare-vcf.py --missing 0.2 | 
        awk '$1 !~ /^N/ && $1 !~ "##contig=<ID=N" {print}' \
        > $cactusVcf_name.filtered_ids.vcf \
        2> $cactusVcf_name.filtered.vcf.log

    cat $cactusVcf_name.filtered_ids.vcf | 
        python3 $pangenie_scripts/add-ids.py \
        > $cactusVcf_name.callset.vcf \
        2> $cactusVcf_name.callset.vcf.log

    bcftools norm --threads $nthreads -m- $cactusVcf_name.callset.vcf > $cactusVcf_name.callset-biallelic.vcf \
    2> $cactusVcf_name.callset-biallelic.vcf.log

    python3 $pangenie_scripts/merge_vcfs.py merge \
        -vcf $cactusVcf_name.callset-biallelic.vcf \
        -r $ref_fa -ploidy 2 > $cactusVcf_name.pangenome.vcf \
        2> $cactusVcf_name.pangenome.vcf.log
}

function PanGenie_index () {
    # cd $ref_path
    # if [[] 0 ]]; then
    #     chr=$(seq 1 29)" X Y MT"
    #     for i in $chr ; do 
    #     samtools faidx ARS_UCD_v2.0.fa ${i} > chr${i}.fa
    #     done
    #     cat chr*fa > ARS_UCD_v2.0.chr.fa
    #     rm chr*fa
    # fi
     
    mkdir -p pangenieIndex_${fa}_${cactusVcf_name}/

    
    PanGenie-index -v $pangenieVcf \
       -r $ref_fa \
       -t $nthreads \
       -o $pangenieIndex
}

function PanGenie_genotyping () {
    PanGenie -f $pangenieIndex \
        -i <(zcat $clean_path/${sample}_1.clean.fq.gz $clean_path/${sample}_2.clean.fq.gz) \
        -s $sample -j $nthreads -t $nthreads -o ${sample}_graph
}

function PanGenie_genotyping_biallelic () { # do not run in this case
    cat ${sample}_graph_genotyping.vcf | 
        python3 $pangenie_scripts/convert-to-biallelic.py \
        $pangenieVcf > ${sample}_genotypes-biallelic.vcf

    bgzip ${sample}_genotypes-biallelic.vcf
    tabix -p vcf ${sample}_genotypes-biallelic.vcf.gz

    bcftools norm -m+ ${sample}_genotypes-biallelic.vcf.gz -Oz -o ${sample}_genotypes.vcf.gz
}

function PanGenie_genotyping_merge3 () {

    bcftools merge --threads $nthreads \
        $(ls $work_path/*/3.pan_data_jer/*_graph_genotyping.vcf.gz) \
        -Oz -o ${cactusVcf_name}_graph_genotyping.merge.vcf.gz
    tabix -p vcf ${cactusVcf_name}_graph_genotyping.merge.vcf.gz
   
    bcftools norm --threads $nthreads -m- ${cactusVcf_name}_graph_genotyping.merge.vcf.gz -Oz -o ${cactusVcf_name}_graph_genotyping.merge-biallelic.vcf.gz
    tabix -p vcf ${cactusVcf_name}_graph_genotyping.merge-biallelic.vcf.gz
}

main


