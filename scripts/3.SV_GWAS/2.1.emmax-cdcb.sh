# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
# --------------------------



cd ${PROJECT_ROOT}/pangenie_cdcb

### +=== EMMAX Kin by SNP thined one in 10kb

plink --vcf cohort.chr_1.snp.vcf.gz --maf 0.05 --bp-space 10000 --chr 1 --make-bed --out cohort.chr_1.snp.maf5.thin10000
# 15122 variants remained
plink --bfile cohort.chr_1.snp.maf5.thin10000 --recode12 --output-missing-genotype 0 --transpose --out cohort.chr_1.snp.maf5.thin10000.t --chr-set 29 
./emmax-kin-intel64 -v -s -d 10 cohort.chr_1.snp.maf5.thin10000.t

# wget https://csg.sph.umich.edu//kang/emmax/download/emmax-intel-binary-20120210.tar.gz
# tar xf emmax-intel-binary-20120210.tar.gz

##### SV dataset

code=${PROJECT_ROOT}/imputation/1.holPub_imp/00.codes
cd ${PROJECT_ROOT}/pangenie_cdcb/3.pan_all/

bcftools annotate --threads $nthreads --remove QUAL,FORMAT \
    hol-pg2hic-2024-05-22_graph_genotyping.merge-biallelic.vcf.gz \
    -Oz -o All2_BTAN_WGS.vcf.gz && 
    tabix -f -p vcf All2_BTAN_WGS.vcf.gz

bcftools view --threads $nthreads -v indels \
    -r $(seq 1 29 | tr '\n' ',') -i 'AC > 0' \
    All2_BTAN_WGS.vcf.gz | 
    awk -f  $code/pan.indel.awk |
    bgzip -c > All2_BTAN_WGS-indels.vcf.gz &&
    tabix -f -p vcf All2_BTAN_WGS-indels.vcf.gz &

bcftools view --threads $nthreads -v other,mnp \
    -r $(seq 1 29 | tr '\n' ',') -i 'AC > 0' \
    All2_BTAN_WGS.vcf.gz  | 
    awk -f  $code/pan.complex.awk |
    bgzip -c > All2_BTAN_WGS-complex.vcf.gz &&
    tabix -f -p vcf All2_BTAN_WGS-complex.vcf.gz &
wait

bcftools concat --threads $nthreads \
    --allow-overlaps \
    All2_BTAN_WGS-complex.vcf.gz \
    All2_BTAN_WGS-indels.vcf.gz \
    -Oz -o All2_BTAN_WGS-nonsnps.vcf.gz &&
    tabix -f -p vcf All2_BTAN_WGS-nonsnps.vcf.gz

zcat All2_BTAN_WGS-nonsnps.vcf.gz | 
    awk '$1 ~ /#/{print; next}{split($3,a,/:|-|_/);if(a[5] > 50) {print}}' |
    bgzip -c > All2_BTAN_WGS-sv.vcf.gz &&
    tabix -f -p vcf All2_BTAN_WGS-sv.vcf.gz

bcftools view --threads $nthreads -v snps \
    -r $(seq 1 29 | tr '\n' ',') -i 'AC > 0' \
    All2_BTAN_WGS.vcf.gz  |
    awk -f  $code/pan.snp.awk |
    bgzip -c > All2_BTAN_WGS-snps.vcf.gz && 
    tabix -f -p vcf All2_BTAN_WGS-snps.vcf.gz
bcftools sort All2_BTAN_WGS-snps2.vcf.gz -Oz -o All2_BTAN_WGS-snps.vcf.gz
 tabix -f -p vcf All2_BTAN_WGS-snps.vcf.gz

bcftools sort All2_BTAN_WGS-snps.vcf.gz -Oz -o All2_BTAN_WGS-snps2.vcf.gz &&
    mv All2_BTAN_WGS-snps2.vcf.gz All2_BTAN_WGS-snps.vcf.gz &&
    tabix -f -p vcf All2_BTAN_WGS-snps.vcf.gz

bcftools concat --threads 4 \
    --allow-overlaps \
    All2_BTAN_WGS-nonsnps.vcf.gz \
    All2_BTAN_WGS-snps.vcf.gz \
    -Oz -o All2_BTAN_WGS-var.vcf.gz &&
    tabix -f -p vcf All2_BTAN_WGS-var.vcf.gz

bcftools annotate --threads 4 --remove QUAL,INFO,FORMAT All2_BTAN_WGS-sv.vcf.gz | 
    grep -v "##" | cut -f 3,10- | gzip -c > All2_BTAN_WGS.sv.bed.gz

cd ${PROJECT_ROOT}/pangenie_cdcb/9.emmax

plink2  --threads $nthreads --vcf ../3.pan_all/All2_BTAN_WGS-sv.vcf.gz --const-fid 0 --chr-set 29 --chr 1-29 \
    --keep-allele-order --mind 0.1 --maf 0.05 --make-bed --out All2_BTAN_WGS.sv.bfile.maf5

plink2  --threads $nthreads --vcf ../3.pan_all/All2_BTAN_WGS-var.vcf.gz --const-fid 0 --chr-set 29 --chr 1-29 \
    --keep-allele-order --mind 0.1 --maf 0.05 --make-bed --out All2_BTAN_WGS.var.bfile.maf5


plink2 --threads $nthreads --bfile All2_BTAN_WGS.sv.bfile.maf5 \
    --chr-set 29 --chr 1-29 --freq --out All2_BTAN_WGS.sv.bfile.maf5.allele_freq

plink2 --threads $nthreads --bfile All2_BTAN_WGS.var.bfile.maf5 \
    --chr-set 29 --chr 1-29 --freq --out All2_BTAN_WGS.var.bfile.maf5.allele_freq

# SV 61249 8.163399e-07
# var 12656969 p 3.950393e-09
# gsnp 11870248 p 4.212212e-09


plink --bfile All2_BTAN_WGS.sv.bfile.maf5 --recode12 --output-missing-genotype 0 --transpose --out All2_BTAN_WGS.sv.bfile.maf5.t --chr-set 29 --threads $nthreads
plink --bfile All2_BTAN_WGS.var.bfile.maf5 --recode12 --output-missing-genotype 0 --transpose --out All2_BTAN_WGS.var.bfile.maf5.t --chr-set 29 --threads $nthreads

#plink --bfile All2_BTAN_WGS.sv.bfile.maf5 --output-missing-genotype 0 --recode transpose --out All2_BTAN_WGS.sv.bfile.maf5.t --chr-set 29 

trait="milk"
nt=8

ls ind_154animals_dPTA.*.txt | cut -d. -f2 | while read trait; do
    echo "./emmax-intel64 -v -d 10 -t All2_BTAN_WGS.sv.bfile.maf5.t -p ind_154animals_dPTA.$trait.txt -k cohort.chr_1.snp.maf5.thin10000.t.aIBS.kinf -o GWAS173_PG.${trait}_sv.emx_k 
./emmax-intel64 -v -d 10 -t All2_BTAN_WGS.var.bfile.maf5.t -p ind_154animals_dPTA.$trait.txt -k cohort.chr_1.snp.maf5.thin10000.t.aIBS.kinf -o GWAS173_PG.${trait}_var.emx_k"
done > run_line3.txt

cat run_line3.txt | while read line; do
let "n_run+=1"
    sbatch -A ${SLURM_ACCOUNT} \
            -D $PWD \
        --export=ALL \
        -J gwas_$n_run \
        -c $nt \
        -o logs/gwas_$n_run.out \
        --wrap="
$line
        "
done
