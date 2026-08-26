# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PANEL_DIR:?PANEL_DIR is unset - see config.sh.example at the repository root}"
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
: "${REF_DIR:?REF_DIR is unset - see config.sh.example at the repository root}"
: "${SOFTWARE_DIR:?SOFTWARE_DIR is unset - see config.sh.example at the repository root}"
# --------------------------

cd ${SOFTWARE_DIR}
mamba create -n cactus  -c conda-forge python=3
conda activate cactus
git clone https://github.com/ComparativeGenomicsToolkit/cactus.git --recursive
cd cactus
source python
pip install -U setuptools pip wheel
pip install -U .
pip install -U -r ./toil-requirement.txt

####

## 2. Run cactus-pangenome
cactusVcf_name="hol-pg2hic-2024-05-22"
cd ${PROJECT_ROOT}/minigraph-cactus
###  Create a sequence file  $cactusVcf_name-seqfile.txt
ls ${PROJECT_ROOT}/0.hifiasm/s*.bp.hap?.p_ctg.fa | egrep -v "4494|4611" |
    awk '{a=$1;gsub(/.*\//,"",a);gsub(/\..*/,"",a);b[a]+=1;print a"."b[a]"\t"$1}' >  $cactusVcf_name-seqfile.txt
echo -e "bosTau9\t${REF_DIR}/ARS_UCD_v2.0.fa" >>  $cactusVcf_name-seqfile.txt
ls  ${PROJECT_ROOT}/hifi-hic-assembly/s*_hi2c.hic.hap?.p_ctg.fa | 
    awk '{a=$1;gsub(/.*\//,"",a);gsub(/\..*/,"",a);b[a]+=1;print a"."b[a]"\t"$1}' >>  $cactusVcf_name-seqfile.txt


## 2. Run cactus-pangenome
cactusVcf_name="holPri-2024-12-03"
echo -e "bosTau9\t${REF_DIR}/ARS_UCD_v2.0.fa" >  $cactusVcf_name-seqfile.txt

ls ${PANEL_DIR}/sample_*/2.assembly/sample_*.bp.p_ctg.fa | egrep -v "4494|4611" |
    awk '{a=$1;gsub(/.*\//,"",a);gsub(/\..*/,"",a);b[a]+=1;print a"\t"$1}' >> $cactusVcf_name-seqfile.txt

ls ${PROJECT_ROOT}/hifi-hic-assembly/s*_hi2c.hic.p_ctg.fa | 
    awk '{a=$1;gsub(/.*\//,"",a);gsub(/\..*/,"",a);b[a]+=1;print a"\t"$1}' >>  $cactusVcf_name-seqfile.txt

## 2. Run cactus-pangenome
cactusVcf_name="bovinePan-2024-08-13"
cd ${PROJECT_ROOT}/minigraph-cactus
###  Create a sequence file  $cactusVcf_name-seqfile.txt
ls ${REF_DIR}/bostauFaDownload/*.fna.gz |
     grep -v "ARS-UCD2.0
GCA_000003055
GCA_000003205" | 
    awk '{a=$1;gsub(/.*\//,"",a);gsub(/\..*/,"",a);print a"\t"$1}' >  $cactusVcf_name-seqfile.txt
ls /${REF_DIR}/bostauPanDownload/*_hifiasm.fa |  grep -v "ARS-UCD2.0
G_hifiasm" | awk '{a=$1;gsub(/.*\//,"",a);gsub(/\..*/,"",a);print a"\t"$1}' >> $cactusVcf_name-seqfile.txt
echo -e "bosTau9\t${REF_DIR}/ARS_UCD_v2.0.fa" >>  $cactusVcf_name-seqfile.txt


## 2. Run cactus-pangenome
cactusVcf_name="jerHap-pg-2024-12-18"
cd ${PROJECT_ROOT}/minigraph-cactus
###  Create a sequence file  $cactusVcf_name-seqfile.txt
ls ${PANEL_DIR}/jer_*/2.assembly/jer_*.bp.hap?.p_ctg.fa | grep -v "bak" |
    awk '{a=$1;gsub(/.*\//,"",a);gsub(/\..*/,"",a);b[a]+=1;print a"."b[a]"\t"$1}' \
    >  $cactusVcf_name-seqfile.txt
echo -e "bosTau9\t${REF_DIR}/ARS_UCD_v2.0.fa" >>  $cactusVcf_name-seqfile.txt

cactusVcf_name="jerPri-pg-2024-12-18"
cd ${PROJECT_ROOT}/minigraph-cactus
###  Create a sequence file  $cactusVcf_name-seqfile.txt
ls ${PANEL_DIR}/jer_*/2.assembly/jer_*.bp.p_ctg.fa | grep -v "bak" |
    awk '{a=$1;gsub(/.*\//,"",a);gsub(/\..*/,"",a);b[a]+=1;print a"\t"$1}' \
    >  $cactusVcf_name-seqfile.txt
echo -e "bosTau9\t${REF_DIR}/ARS_UCD_v2.0.fa" >>  $cactusVcf_name-seqfile.txt


### clean toil 
toil clean $cactusVcf_name-seqfile.js
export PATH=${SOFTWARE_DIR}/cactus-bin-v2.7.1/bin:$PATH
export SLURM_ACCOUNT=${SLURM_ACCOUNT:?SLURM_ACCOUNT is unset - see config.sh.example}
export SALLOC_ACCOUNT=$SLURM_ACCOUNT
export SBATCH_ACCOUNT=$SLURM_ACCOUNT
export SALLOC_PARTITION=${SLURM_PARTITION}
export SBATCH_PARTITION=$SALLOC_PARTITION
export SALLOC_NODES=2
export SBATCH_NODES=$SALLOC_NODES
export SALLOC_CPUS_PER_TASK=32
export SBATCH_CPUS_PER_TASK=$SALLOC_CPUS_PER_TASK
export SALLOC_MEM=256G
export SBATCH_MEM=$SALLOC_MEM

sa cactus
work_path="${PROJECT_ROOT}/minigraph-cactus"
cd $work_path
mkdir -p ${work_path}/tmp
export TMPDIR=$work_path/tmp
export PATH=${SOFTWARE_DIR}/cactus-bin-v2.7.1/bin:$PATH

### run cactus-pangenome 
####  $cactusVcf_name

cactus-pangenome --binariesMode local \
    $cactusVcf_name.js $cactusVcf_name-seqfile.txt \
    --outDir  $cactusVcf_name \
    --outName  $cactusVcf_name \
    --reference bosTau9 \
    --giraffe clip filter \
    --gbz clip filter full \
    --gfa clip filter full \
    --vcf --permissiveContigFilter \
    --haplo --chrom-vg clip filter \
    --chrom-og full --viz \
    --consCores 48 --consMemory 380G \
    --indexCores 48 --indexMemory 380G \
    --mgCores 48 --mgMemory 380G \
    --mapCores 24 \
    --batchSystem slurm \
    --logFile  $cactusVcf_name.log \
    2>  $cactusVcf_name.stderr &
# --restart
# --refContigs $(for i in $(seq 29); do printf "chr$i "; done ; echo "chrX chrY chrM") 
# --permissiveContigFilter
# // --odgi

