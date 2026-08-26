# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
: "${SLURM_ACCOUNT:?SLURM_ACCOUNT is unset - see config.sh.example at the repository root}"
# --------------------------

cd ${PROJECT_ROOT}/pangenie_cdcb/9.emmax


plink2  --threads $nthreads --vcf ../3.pan_all/All2_BTAN_WGS-snps.vcf.gz --const-fid 0 --chr-set 29 --chr 1-29 \
    --keep-allele-order --mind 0.1 --maf 0.05 --make-bed --out All2_BTAN_WGS.var-snp.bfile.maf5


shuf -n 61249 All2_BTAN_WGS.var-snp.bfile.maf.bim | cut -f 2 > All2_BTAN_WGS.var-snp.bfile.maf.snplist

plink2 --threads $nthreads --bfile All2_BTAN_WGS.var-snp.bfile.maf --extract All2_BTAN_WGS.var-snp.bfile.maf.snplist --const-fid 0 --chr-set 29 --chr 1-29  --make-bed --out All2_BTAN_WGS.var-snp-dsv.bfile.maf

for typ in var-snp-dsv var-snp; do
plink --threads $nthreads --bfile All2_BTAN_WGS.$typ.bfile.maf --recode12 --output-missing-genotype 0 --transpose --out All2_BTAN_WGS.$typ.bfile.maf5.t --chr-set 29 
done

trait="fat"
nt=8
for typ in var-snp-dsv var-snp; do
    ls ind_154animals_dPTA.*.txt | cut -d. -f2 | while read trait; do
    echo "./emmax-intel64 -v -d 10 -t All2_BTAN_WGS.$typ.bfile.maf5.t -p ind_154animals_dPTA.$trait.txt -k cohort.chr_1.snp.maf5.thin10000.t.aIBS.kinf -o GWAS173_PG.${trait}_$typ.emx_k"
done 
done > run_line1.txt

cat run_line1.txt | while read line; do
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

#10597560
#SV 61249 8.163399e-07

thr=4.71807e-09
cat GWAS173_PG.fat_var-snp.emx_k.ps | awk '($4 < Thr*1e1 && NR % 2 == 1) || ($4 < Thr*1e2 && NR % 5 == 1) || ($4 < Thr*1e3 && NR % 10 == 1) || ($4 > Thr*1e3 && NR % 200 == 1)' Thr=$thr > sig.fat_var-snp.emx_k.base.ps

cp GWAS173_PG.fat_var-snp-dsv.emx_k.ps sig.fat_var-snp-dsv.emx_k.base.ps
thr=3.950393e-09
