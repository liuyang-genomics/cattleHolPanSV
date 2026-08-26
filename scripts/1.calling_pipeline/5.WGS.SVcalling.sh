#!/bin/bash
# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
: "${SCRATCH_DIR:?SCRATCH_DIR is unset - see config.sh.example at the repository root}"
: "${REF_DIR:?REF_DIR is unset - see config.sh.example at the repository root}"
: "${EXT_PROJECT_DIR:?EXT_PROJECT_DIR is unset - see config.sh.example at the repository root}"
: "${SOFTWARE_DIR:?SOFTWARE_DIR is unset - see config.sh.example at the repository root}"
: "${CONDA_BASE:?CONDA_BASE is unset - see config.sh.example at the repository root}"
# --------------------------

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

    mkdir -p $svtools_all
    cd $svtools_all
    
    #module load miniconda3
    #source activate $conda_svtools

    ## run for individual

    smoove_call

    ## run for all
    svtools_merge(){
        svtools lsort \
            $(ls $work_path/*/5.svtools_data/2.lumpy_vcf/*.part1-smoove.genotyped.vcf.gz | 
                grep -f $sample_good) | 
            svtools lmerge -i /dev/stdin -f 20 > $svtools_all/merged.total.sites.vcf
        create_coordinates -i $svtools_all/merged.total.sites.vcf -o $svtools_all/coordinates
    }


    ## run for individual
	svtools_genotype
    check1_svcount
    run_cnvnator
    call_cnvnator
    cnvnator_genotype
    cnvnatorsex
    svtools_copynumber
    check2_svcount


    # mkdir -p 2.lumpy_vcf 3.genotype_vcf 0.refChr 4.cnvnator_root 5.cnvnator_rawcnv 6.cnvnator_copynumber 7.copynumber_vcf
    ## run for all
    svtools_vcfpaste
	svtools_prune


}


install_packages () {
    mamba install -c bioconda -c conda-forge -c default speedseq root cnvnator svtools smoove samtools=1.9 bgzip vawk 
    pip install git+https://github.com/hall-lab/svtyper.git

    conda create -n svtools -c bioconda -c conda-forge -c default speedseq root cnvnator svtools smoove samtools bgzip vawk
    
}

smoove_call(){
	mkdir -p 2.lumpy_vcf
	
    smoove call \
        -p $nthreads \
        --outdir 2.lumpy_vcf \
        --exclude $ref_gap \
        --name ${sample}.part1 \
        --fasta $ref_fa \
        --genotype $bam_path/${sample}.sorted.bam 
}


svtools_merge(){
    svtools lsort \
        $(ls $work_path/*/5.svtools_data/2.lumpy_vcf/*.part1-smoove.genotyped.vcf.gz | 
            grep -f $sample_good) | 
        svtools lmerge -i /dev/stdin -f 20 > $svtools_all/merged.total.sites.vcf
    create_coordinates -i $svtools_all/merged.total.sites.vcf -o $svtools_all/coordinates
}


svtools_genotype(){
    mkdir -p 3.genotype_vcf
    
    cat $svtools_all/merged.total.sites.vcf \
        | vawk --header '{  $6="."; print }' \
        | svtools genotype \
          -B $bam_path/${sample}.sorted.bam \
          -l 3.genotype_vcf/${sample}.bam.json \
        | sed 's/PR...=[0-9\.e,-]*\(;\)\{0,1\}\(\t\)\{0,1\}/\2/g' - \
        > 3.genotype_vcf/${sample}.genotype.vcf 
}

check1_svcount(){
    wc -l 3.genotype_vcf/${sample}.genotype.vcf | 
        sed "s/^/$sample\t/" >> $svtools_all/part1.genotype.vcf.count
}

run_cnvnator(){
    mkdir -p 4.cnvnator_root
    if [[ -f 4.cnvnator_root/${sample}.root ]]; then
        rm 4.cnvnator_root/${sample}.root
    fi
    cnvnator -root 4.cnvnator_root/${sample}.root \
        -tree $bam_path/${sample}.sorted.bam \
        -chrom $(seq 1 29) X Y

    cnvnator -root 4.cnvnator_root/${sample}.root \
        -his $cnvnatorbin \
        -d $cnvnator_ref
            
    cnvnator -root 4.cnvnator_root/${sample}.root \
        -stat $cnvnatorbin
            
    cnvnator -root 4.cnvnator_root/${sample}.root \
        -partition $cnvnatorbin
}

call_cnvnator () {
    mkdir -p 5.cnvnator_rawcnv
    cnvnator -root 4.cnvnator_root/${sample}.root \
        -call $cnvnatorbin | awk -v OFS="\t" -v sample=$sample '{if($5 <= 0.05 && $9 <= 0.5){gsub("[:-]","\t",$2);print($2,$3,$1,$4,sample)}}' > 5.cnvnator_rawcnv/${sample}.cnvnator.bed

    # tittle of CNVnator output: CNV_type coordinates CNV_size normalized_RD e-val1 e-val2 e-val3 e-val4 q0 
    # normalized_RD -- read depth normalized to 1.
    # e-val1 -- is calculated using t-test statistics.
    # e-val2 -- is from the probability of RD values within the region to be in the tails of a gaussian distribution describing frequencies of RD values in bins.
    # e-val3 -- same as e-val1 but for the msampledle of CNV
    # e-val4 -- same as e-val2 but for the msampledle of CNV
    # q0 -- fraction of reads mapped with q0 quality
}

cnvnator_genotype(){
    mkdir -p 6.cnvnator_copynumber
    # zcat merged.total.sites.vcf.gz | create_coordinates -i - -o coordinates 
    ## end exit 
    
    cat $svtools_all/coordinates | \
        cnvnator -root 4.cnvnator_root/${sample}.root \
        -genotype $cnvnatorbin > 6.cnvnator_copynumber/$sample.cnvnator.copynumber 2>/dev/null
}

cnvnatorsex(){
    if [[ -f "6.cnvnator_copynumber/$sample.cnvnator.copynumber" ]] ; then
        ch_sex=$(cat 6.cnvnator_copynumber/$sample.cnvnator.copynumber | grep "Assuming male indivsampleual" | head -n1)
        if [[ -n "$ch_sex" ]] ; then 
            echo "$sample 1" | sed 's/ /\t/g' >> $svtools_all/ceph.sex.txt
        else 
            echo "$sample 2" | sed 's/ /\t/g' >> $svtools_all/ceph.sex.txt
        fi
    else 
        echo "$sample none" | sed 's/ /\t/g' >> $svtools_all/ceph.sex.txt
    fi
}

svtools_copynumber(){
    # change /public/agis/likui_group/shuli/0.ly/miniconda3/envs/smoove/lib/python2.7/site-packages/svtools/copynumber.py for cnvnator genotyped input
    #l=$(wc -l 3.genotype_vcf/${sample}.genotype.vcf)
    #if [[ $l != 478103 ]]; then
    #    cat 3.genotype_vcf/${sample}.genotype.vcf | awk '$1 ~ "^#" || $1 !~ "M|A|F"' > 3.genotype_vcf/${sample}.genotype.rmUnmap.vcf
    #    mv 3.genotype_vcf/${sample}.genotype.rmUnmap.vcf 3.genotype_vcf/${sample}.genotype.vcf
    #fi
    
    mkdir -p 7.copynumber_vcf
    
    svtools copynumber \
        --cnvnator cnvnator \
        -s $sample \
        -w $cnvnatorbin \
        -c $svtools_all/coordinates \
        -r 4.cnvnator_root/${sample}.root \
        -i 3.genotype_vcf/${sample}.genotype.vcf \
        > 7.copynumber_vcf/$sample.copynumber.vcf

        # -g 6.cnvnator_copynumber/$sample.cnvnator.copynumber \
}

check2_svcount(){
    wc -l 7.copynumber_vcf/${sample}.copynumber.vcf | sed "s/^/$sample\t/" >> $svtools_all/part7.copynumber.vcf.count
}

svtools_vcfpaste(){
    #wget https://github.com/hall-lab/svtools/blob/master/resources/training_vars.bedpe.gz
    ls -1 $work_path/*/5.svtools_data/7.copynumber_vcf/*copynumber.vcf | 
        grep -f $sample_good > cn.list
    svtools vcfpaste \
        -m merged.total.sites.vcf \
        -f cn.list \
        -q \
        | bgzip -c -@ $nthreads \
        > merged.sv.gt.cn.vcf.gz
}

svtools_prune(){
    zcat merged.sv.gt.cn.vcf.gz \
        | svtools afreq \
        | svtools vcftobedpe \
        | svtools bedpesort \
        | svtools prune -s -d 100 -e "AF" \
        | svtools bedpetovcf \
        | bgzip -c -@ $nthreads > merged.sv.pruned.vcf.gz
    #239192
}


main


