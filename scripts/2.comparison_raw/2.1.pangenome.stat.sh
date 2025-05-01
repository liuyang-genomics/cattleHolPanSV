




cd ${PROJECT_ROOT}/pangenie_HiFi/3.pan_all

cat > vcf.svtype.sh <<'EOF'
bcftools query -f "%ID %INFO/ID\n" | 
    awk '{if($2 ~ /,|:/){split($2,a,/,|:/);for(i in a){gsub("-"," ",a[i]);print $1,a[i]}} else {gsub("-"," ",$2);print}}' | 
    awk '{print $1,$2,$3,$5,$4,$6}'
EOF

cat > count.svtype.sh <<'EOF'
awk '{gc[$5]++;gl[$5]+=$6;if($2 <= 29){ac[$5]++;al[$5]+=$6}} END {for(i in gc) print idd,i,gc[i],gl[i],ac[i],al[i]}' idd=$1
EOF


ls *vcf.gz |
    while read id; do
        idd=$(basename $id);
        ale=${idd#*.};
        ale=${ale%.*};
        idd=${idd%%-*};
 sbatch -A ${SLURM_ACCOUNT} \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --time=02:00:00 \
            --wrap="
#zcat $id | bash vcf.svtype.sh > ${id/vcf.gz/svtype}
cat ${id/vcf.gz/svtype} | bash count.svtype.sh $idd.$ale > ${id/vcf.gz/svtype}.stats
"
    done




ls *vcf.gz |
    while read id; do
        idd=$(basename $id);
        ale=${idd#*.};
        ale=${ale%.*};
        idd=${idd%%-*};

#zcat $id | bash vcf.svtype.sh > ${id/vcf.gz/svtype}
cat ${id/vcf.gz/svtype} | bash count.svtype.sh $idd.$ale > ${id/vcf.gz/svtype}.stats &

    done