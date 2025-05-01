ls ${PROJECT_ROOT}/minigraph-cactus/*/*.vcf.gz | grep -v 'raw' | 
    while read id; do
        sbatch -A ${SLURM_ACCOUNT} \
        --cpus-per-task=8 \
        --mem-per-cpu=8g \
        --wrap="bcftools view --threads 8 -i 'F_MISSING<0.2' $id -Oz -o ${id/vcf.gz/.filtered.vcf.gz}"
    done
    

ls *.filter.gz |
    while read id; do
        sbatch -A ${SLURM_ACCOUNT} \
        --cpus-per-task=8 \
        --mem-per-cpu=8g \
        --wrap="bcftools norm --threads 8 -m- $id -Oz -o ${id/filtered/filtered-bi}"
    done

cat ../minigraph-cactus/allbovinePan-2024-07-03/allbovinePan-2024-07-03-seqfile.txt | awk '$1 ~ /sample/ {gsub(/\..*/,"",$1);print $1}' | sort -u > hol.sample

cat ../minigraph-cactus/allbovinePan-2024-07-03/allbovinePan-2024-07-03-seqfile.txt | awk '$1 ~ /Jersey/ {gsub(/\..*/,"",$1);print $1}' | sort -u > jer.sample

cat ../minigraph-cactus/allbovinePan-2024-07-03/allbovinePan-2024-07-03-seqfile.txt | awk '$1 ~ /hifiasm/ || $1 ~ /GCA/ {gsub(/\..*/,"",$1);print $1}' | sort -u > pub.sample

cat hol.sample jer.sample pub.sample | sort -u > all.sample

sbatch -A ${SLURM_ACCOUNT} --cpus-per-task=8 --wrap="bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -S hol.sample allbovinePan-2024-07-03.filtered-bi.vcf.gz -Oz -o allbovinePan-2024-07-03_hol.filtered-bi.vcf.gz"
sbatch -A ${SLURM_ACCOUNT} --cpus-per-task=8 --wrap="bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -S jer.sample allbovinePan-2024-07-03.filtered-bi.vcf.gz -Oz -o allbovinePan-2024-07-03_jer.filtered-bi.vcf.gz"
sbatch -A ${SLURM_ACCOUNT} --cpus-per-task=8 --wrap="bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -S pub.sample allbovinePan-2024-07-03.filtered-bi.vcf.gz -Oz -o allbovinePan-2024-07-03_pub.filtered-bi.vcf.gz"

sbatch -A ${SLURM_ACCOUNT} --cpus-per-task=8 --wrap="bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -S hol.sample allbovinePan-2024-07-03.filtered.vcf.gz -Oz -o allbovinePan-2024-07-03_hol.filtered.vcf.gz"
sbatch -A ${SLURM_ACCOUNT} --cpus-per-task=8 --wrap="bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -S jer.sample allbovinePan-2024-07-03.filtered.vcf.gz -Oz -o allbovinePan-2024-07-03_jer.filtered.vcf.gz"
sbatch -A ${SLURM_ACCOUNT} --cpus-per-task=8 --wrap="bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -S pub.sample allbovinePan-2024-07-03.filtered.vcf.gz -Oz -o allbovinePan-2024-07-03_pub.filtered.vcf.gz"

sbatch -A ${SLURM_ACCOUNT} --cpus-per-task=8 --wrap="bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -s sample_4232 bovinePanPhase-2024-07-03.filtered-bi.vcf.gz -Oz -o bovinePanPhase-2024-07-03_hol.filtered-bi.vcf.gz"
sbatch -A ${SLURM_ACCOUNT} --cpus-per-task=8 --wrap="bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -s Jersey_441 bovinePanPhase-2024-07-03.filtered-bi.vcf.gz -Oz -o bovinePanPhase-2024-07-03_jer.filtered-bi.vcf.gz"
sbatch -A ${SLURM_ACCOUNT} --cpus-per-task=8 --wrap="bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -S pub.sample bovinePanPhase-2024-07-03.filtered-bi.vcf.gz -Oz -o bovinePanPhase-2024-07-03_pub.filtered-bi.vcf.gz"

sbatch -A ${SLURM_ACCOUNT} --cpus-per-task=8 --wrap="bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -s sample_4232 bovinePanPhase-2024-07-03.filtered.vcf.gz -Oz -o bovinePanPhase-2024-07-03_hol.filtered.vcf.gz"
sbatch -A ${SLURM_ACCOUNT} --cpus-per-task=8 --wrap="bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -s Jersey_441 bovinePanPhase-2024-07-03.filtered.vcf.gz -Oz -o bovinePanPhase-2024-07-03_jer.filtered.vcf.gz"
sbatch -A ${SLURM_ACCOUNT} --cpus-per-task=8 --wrap="bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -S pub.sample bovinePanPhase-2024-07-03.filtered.vcf.gz -Oz -o bovinePanPhase-2024-07-03_pub.filtered.vcf.gz"

ls *.filter*vcf.gz | grep -v "$(ls *.filter*vcf.gz.tbi | sed 's/.tbi//')" |
    while read id; do
        sbatch -A ${SLURM_ACCOUNT} \
        --cpus-per-task=2 \
        --mem-per-cpu=8g \
        --wrap="tabix -p vcf $id"
    done

ls *.filter*vcf.gz | 
    while read id; do
        sbatch -A ${SLURM_ACCOUNT} \
        --cpus-per-task=4 \
        --mem-per-cpu=8g \
        --wrap="bcftools stats --threads $id > ${id/gz/stats}"
    done

######### 2.1.sec #####

cat all.sample | while read id;
do 
ls 2.0.filtered/*.filtered.vcf.gz | while read file;
do 
file_name=$(basename $file)
file_suffix=${file_name#*.}
file_name=${file_name%%.*}
    sbatch -A ${SLURM_ACCOUNT} \
    --cpus-per-task=8 \
    --wrap="bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -s $id $file -Oz -o 2.2.filtered.ind/${file_name}_${id}.$file_suffix"
done
done


cat all.sample | while read id;
do 
ls 3.0.bi-filtered/*.filtered-bi.vcf.gz | while read file;
do 
file_name=$(basename $file)
file_suffix=${file_name#*.}
file_name=${file_name%%.*}
    sbatch -A ${SLURM_ACCOUNT} \
    --cpus-per-task=8 \
    --wrap="bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -s $id $file -Oz -o 3.2.bi-filtered.ind/${file_name}_${id}.$file_suffix"
done
done


mkdir -p 4.vcf.stats

ls */*.filter*vcf.gz | 
    while read id; do
        idd=$(basename $id);
        idd=${idd/gz/stats}
        sbatch -A ${SLURM_ACCOUNT} \
        -p ${SLURM_PARTITION} \
        --cpus-per-task=4 \
        --mem-per-cpu=8g \
        --wrap="bcftools stats --threads 4 $id > 4.vcf.stats/$idd"
    done

cat > 4.vcf.stats/vcf.filtered.stats.tab <<EOF
ID samples records no-ALTs SNPs MNPs indels others multiallelic-sites multiallelic-SNP
EOF

ls -1 4.vcf.stats/*.filtered.vcf.stats |
    while read id; do
        idd=$(basename $id);
        idd=${idd%.*};
        awk '$1 == "SN"{a[NR]=$NF} 
        END {
            printf "%s ",ID
            for(i in a)printf "%s ",a[i]
            print ""
        }' ID=$idd $id >> 4.vcf.stats/vcf.filtered.stats.tab 
    done

cat > 4.vcf.stats/vcf.filtered-bi.stats.tab <<EOF
ID samples records no-ALTs SNPs MNPs indels others multiallelic-sites multiallelic-SNP
EOF

ls -1 4.vcf.stats/*.filtered-bi.vcf.stats |
    while read id; do
        idd=$(basename $id);
        idd=${idd%%.*};
        awk '$1 == "SN"{a[NR]=$NF} 
        END {
            printf "%s ",ID
            for(i in a)printf "%s ",a[i]
            print ""
        }' ID=$idd $id >> 4.vcf.stats/vcf.filtered-bi.stats.tab 
    done
