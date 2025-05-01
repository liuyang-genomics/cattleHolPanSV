
alias cdc='cd ${PROJECT_ROOT}/stat_pan'
cdc
####################################################################################################################################################
####################################################################################################################################################

##### 0. fa stats #####
wd='0.fa.stats' # ========================================
mkdir -p $wd

ls ${PROJECT_ROOT}/hifi-hic-assembly/s*_hi2c.hic*p_ctg.fa |
    while read id; do
        idd=$(basename $id);
        idd=${idd%.p_ctg.fa};
        idd=${idd/.hifi_hi2c.hic/}
        sbatch -A ${SLURM_ACCOUNT} \
            --wrap="gfastats $id > 0.fa.stats/$idd.gfastats"
    done

cat ${PROJECT_ROOT}/minigraph-cactus/*-seqfile.txt  | 
    awk '{print $2}' | sed -e  's|//|/|' | sort -u 
while read id ; do
    idd=$(basename $id);
    idd=${idd%%.*};
    sbatch -A ${SLURM_ACCOUNT} \
        --wrap="gfastats $id > $wd/$idd.gfastats "
done

l 0.fa.stats/
cat 0.fa.stats/*.gfastats | grep "Total scaffold length"
cat 0.fa.stats/*.gfastats | grep "Contig N50:"


ls ${PANEL_DIR}/*/2.assembly/*.busco-hifiasm-*/short_summary.specific.mammalia_odb10.*.busco-hifiasm*.txt | grep -v bak |
while read id; do
cat $id | awk '$0 ~ "Contigs N50" || $0 ~ "Total length" {a[NR]=$1} END{for (i in a){printf a[i]"\t"};print idd}' idd=$(echo $id | cut -d/ -f8 | sed 's/.busco-hifiasm-bp//'); 
done > 0.fa.stats/busco-hifiasm.stats

### fast sta 0.fa.stats/busco-hifiasm.stats

####################################################################################################################################################
####################################################################################################################################################

##### 1. pan stats #####
cdc

wd=1.pan.stats # ========================================
mkdir -p $wd

ls ${PROJECT_ROOT}/minigraph-cactus/*/*.stats.tgz | grep -v bak |
    while read id; do
        idd=$(basename $id);
        idd=${idd%.stats.tgz};
        mkdir -p mstat/$idd;
        #mv ${idd}* $idd;
        tar -xzf $id -C mstat/$idd;
        ls mstat/$idd/local/scratch/${USER}/*/*/*/*/tmp*/* | 
            while read file; do
                mv $file $wd/${idd}.$(basename $file)
            done
    done
##node stats
### 

sa cactus
wd=1.pan.stats # ========================================
mkdir -p $wd
ls ${PROJECT_ROOT}/minigraph-cactus/*/*.gfa | grep holPri-2024-12-26 |
    while read id; do
        idd=$(basename $id);
        idd=$wd/${idd%.gfa};
    cat > run_hist.sh << EOF
cat $id | grep -e '^W' | cut -f2-6 | awk '$1 !~ "bosTau9{ print $1 "#" $2 "#" $3 ":" $4 "-" $5 }' > $idd.paths.txt ###
cat $idd.paths.txt | awk '$1 !~ "bosTau9"' | awk '{
    if($1 ~ "Jersey") {j[$1]=$0}else if($1 ~ "sample") {h[$1]=$0}else if($1 ~ "B_hifiasm") {b[$1]=$0} else {print}} END {for(i in b){print b[i]}; n=asorti(j, sortedj);for(i in sortedj){print j[sortedj[i]]}; n=asorti(h, sortedh);for(i in sortedh){print h[sortedh[i]]}}' > $idd.paths.haplotypes.txt

    awk -F# '{p=$1"_hap"$2;a[p]++;if(a[p] == 1) {print p}}' $idd.paths.haplotypes.txt > $idd.paths.histgrowth.order.txt

RUST_LOG=info panacus histgrowth -c all -t16 -l 1,2,1,1 -q 0,0,0.05,0.1,0.5,0.9,1 -S -a -s $idd.paths.haplotypes.txt $id > $idd.histgrowth.all.tsv
RUST_LOG=info panacus histgrowth -c all -t16 -l 1,2,1,1 -q 0,0,0.05,0.1,0.5,0.9,1 -H -a -s $idd.paths.haplotypes.txt $id > $idd.histgrowth.all-hap.tsv
EOF
    sbatch -A ${SLURM_ACCOUNT} \
    --cpus-per-task=16 \
    --mem-per-cpu=8g \
    --wrap="
run_hist.sh
        "
    done


################################
# ========================================
cat > hist.order.txt <<EOF
bosTau9
sample_4232
sample_4245
sample_4294
sample_4298
sample_4420
sample_4443
sample_4450
sample_4494
sample_4611
sample_4681
sample_4801
sample_4816
sample_4828
sample_4837
sample_4852
sample_4873
sample_4879
sample_4889
sample_4890
sample_4896
GCA_021347905
GCA_021234555
GCA_028973685
GCA_039881175
GCA_905123515
GCA_947034695
H_hifiasm
O_hifiasm
P_hifiasm
B_hifiasm
GCA_034097375
GCA_905123885
N_hifiasm
jer_426
jer_441
jer_454
jer_456
jer_467
jer_471
jer_477
jer_481
jer_485
jer_486
EOF

cat > hist-hap.order.txt <<EOF
bosTau9#0
sample_4232#1
sample_4232#2
sample_4245#1
sample_4245#2
sample_4294#1
sample_4294#2
sample_4298#1
sample_4298#2
sample_4420#1
sample_4420#2
sample_4443#1
sample_4443#2
sample_4450#1
sample_4450#2
sample_4494#1
sample_4494#2
sample_4611#1
sample_4611#2
sample_4681#1
sample_4681#2
sample_4801#1
sample_4801#2
sample_4816#1
sample_4816#2
sample_4828#1
sample_4828#2
sample_4837#1
sample_4837#2
sample_4852#1
sample_4852#2
sample_4873#1
sample_4873#2
sample_4879#1
sample_4879#2
sample_4889#1
sample_4889#2
sample_4890#1
sample_4890#2
sample_4896#1
sample_4896#2
GCA_021347905#0
GCA_021234555#0
GCA_028973685#0
GCA_039881175#0
GCA_905123515#0
GCA_947034695#0
H_hifiasm#0
O_hifiasm#0
P_hifiasm#0
B_hifiasm#0
GCA_034097375#0
GCA_905123885#0
N_hifiasm#0
jer_426#1
jer_426#2
jer_441#1
jer_441#2
jer_454#1
jer_454#2
jer_456#1
jer_456#2
jer_467#1
jer_467#2
jer_471#1
jer_471#2
jer_477#1
jer_477#2
jer_481#1
jer_481#2
jer_485#1
jer_485#2
jer_486#1
jer_486#2
EOF

sa cactus
wd=1.pan.stats_bp
mkdir -p $wd
countytpe=bp
ls ${PROJECT_ROOT}/minigraph-cactus/*/*.gfa | grep holPri-2024-12-26 |
    while read id; do
        idd=$(basename $id);
        idd=$wd/${idd%.gfa};
        
        #cat $id | grep -e '^W' | cut -f2-6 | awk '$1 !~ "bosTau9{ print $1 "#" $2 "#" $3 ":" $4 "-" $5 }' > $idd.paths.txt ###
        cat $idd.paths.txt | awk '$1 !~ "bosTau9"'> $idd.bosTau9.paths.txt

echo "RUST_LOG=info panacus ordered-histgrowth -c $countytpe -t16 -l 1,2,1 -q 0.0,1 -H -O hist-hap.order.txt -s $idd.bosTau9.paths.txt $id > $idd.ordered-histgrowth.$countytpe-hap.noref.tsv
RUST_LOG=info panacus ordered-histgrowth -c $countytpe -t16 -l 1,2,1 -q 0.0,1 -H -O hist-hap.order.txt  $id > $idd.ordered-histgrowth.$countytpe-hap.adref.tsv
RUST_LOG=info panacus ordered-histgrowth -c $countytpe -t16 -l 1,2,1 -q 0.0,1 -S -O hist.order.txt -s $idd.bosTau9.paths.txt $id > $idd.ordered-histgrowth.$countytpe-pri.noref.tsv
RUST_LOG=info panacus ordered-histgrowth -c $countytpe -t16 -l 1,2,1 -q 0.0,1 -S -O hist.order.txt  $id > $idd.ordered-histgrowth.$countytpe-pri.adref.tsv"
    done > stats_$countytpe.sh

cat stats_$countytpe.sh | while read line; do
    sbatch -A ${SLURM_ACCOUNT} \
    -J hist \
    --cpus-per-task=4 \
    --mem-per-cpu=8g \
    --wrap="
$line
        "
    done

C:\Users\${USER}\OneDrive - University of Maryland\Data\2024-02-07.cattleLR-SR-GWAS\2.analyses\1.1.pan.plot
scp ${CLUSTER_HOST}:${PROJECT_ROOT}/stat_pan/1.pan.stats_bp/*ref.tsv .
# ========================================
################################


##### 1.pan.raw.anno.stats #####
annotatePy="${SOFTWARE_DIR}/genotyping-pipelines/prepare-vcf-MC/workflow/scripts/annotate_vcf.py"

ls ${PROJECT_ROOT}/minigraph-cactus/*/*.vcf.gz | grep -v 'raw' | grep -v 'filter' | grep -v 'bak' | grep -v "anno" |
    while read id; do
    sbatch -A ${SLURM_ACCOUNT} \
        --cpus-per-task=8 \
        --mem-per-cpu=8g \
        --wrap="
gunzip -c $id > ${id/.gz/} 
gunzip -c ${id/vcf/gfa} > ${id/vcf.gz/gfa}
python3 $annotatePy \
    -vcf ${id/.gz/} \
    -gfa ${id/vcf.gz/gfa} \
    -o ${id/vcf.gz/anno} &> ${id/vcf.gz/anno}.log
    bgzip -@ 8 ${id/vcf.gz/anno.vcf}
    bgzip -@ 8 ${id/vcf.gz/anno_biallelic.vcf}
    rm ${id/.gz/}
"
    done

##
# ========================================
mkdir -p x.anno.stats

ls ${PROJECT_ROOT}/minigraph-cactus/*/*.anno.log |  while read id; do 
    idd=$(basename $id);
    idd=${idd%%-*};
    tail -n2 $id | awk '{print idd,$3,$2}' idd=$idd; 
done > x.anno.stats/1.pan.raw.anno.stats

####################################################################################################################################################
####################################################################################################################################################

#### 1.pan.raw.vcf.stats ####

cat > vcf.svtype.sh <<'EOF'
bcftools query -f "%ID %INFO/ID\n" | 
    awk '{if($2 ~ /,|:/){split($2,a,/,|:/);for(i in a){gsub("-"," ",a[i]);print $1,a[i]}} else {gsub("-"," ",$2);print}}' | 
    awk '{print $1,$2,$3,$5,$4,$6}'
EOF


cat > count.svtype.sh <<'EOF'
awk '{gc[$5]++;gl[$5]+=$6;if($2 <= 29){ac[$5]++;al[$5]+=$6}} END {for(i in gc) print idd,i,gc[i],gl[i],ac[i],al[i]}' idd=$1
EOF

cat > max.svtype.sh <<'EOF'
awk '{if(a[$1] < $6){a[$1]=$6;c[$1]=$2}} END {for(i in a){gc++;gl+=a[i];if(c[i] <= 29){ac++;al+=a[i]}}; print idd,"MAX",gc,gl,ac,al}' idd=$1
EOF

ls ${PROJECT_ROOT}/minigraph-cactus/*/*.anno*vcf.gz | grep -v 'raw' | grep -v 'filter' | grep -v 'bak' | grep -v "full" | grep -v "d2" | grep -v "sv" |  grep -v "node" | grep -v Mul | grep -v bov| grep -v hol-pg | 

rm x.anno.stats/2.pan.raw.vcf.stats
ls ${PROJECT_ROOT}/minigraph-cactus/*/*.anno*vcf.gz |
    while read id; do
        idd=$(basename $id);
        ale=${idd#*.};
        ale=${ale%.*};
        idd=${idd%%-*};
        sbatch -A ${SLURM_ACCOUNT} -J $idd \
            --cpus-per-task=2 \
            --mem-per-cpu=8g \
            --wrap="
zcat $id | bash vcf.svtype.sh > ${id/vcf.gz/svtype}
cat ${id/vcf.gz/svtype} | bash count.svtype.sh $idd.$ale >> x.anno.stats/2.pan.raw.vcf.stats
cat ${id/vcf.gz/svtype} | bash max.svtype.sh $idd.$ale >> x.anno.stats/2.pan.raw.vcf.stats
"
    done

####################################################################################################################################################
####################################################################################################################################################

#### 2. vcf stats ####
cdc
wd=2.vcf_stats/0.filter
mkdir -p $wd

ls ${PROJECT_ROOT}/minigraph-cactus/*/*vcf.gz | grep -v 'raw' | grep -v 'filter' | grep -v 'bak' |
    while read id; do
        idd=$(basename $id);
        idd=${idd/vcf/filtered.vcf};
        sbatch -A ${SLURM_ACCOUNT} -J $idd.filter -c 8 --wrap="
bcftools view --threads 8 -i 'F_MISSING<0.2' $id -Oz -o $wd/${idd}
"
    done


## 2.1 vcf grp stats

cat ../minigraph-cactus/bovHol-2024-08-13/bovHol-2024-08-13-seqfile.txt | awk '$1 ~ /sample/ {gsub(/\..*/,"",$1);print $1}' | sort -u > hol.sample

cat ../minigraph-cactus/bovHol-2024-08-13/bovHol-2024-08-13-seqfile.txt | awk '$1 ~ /hifiasm/ || $1 ~ /GCA/ {gsub(/\..*/,"",$1);print $1}' | sort -u > pub.sample

cat hol.sample pub.sample | sort -u > all.sample

wd="2.vcf_stats/1.grp"
mkdir -p $wd

for sample in hol.sample pub.sample; do
    ls 2.vcf_stats/0.filter/bovHol-2024-08-13.*filtered.vcf.gz  |
    while read id; do
        idd=$(basename $id);
        idd_suffix=${idd#*.};
        idd=${idd%%.*};
        sbatch -A ${SLURM_ACCOUNT} -J $idd.grp -c 8 --wrap="
bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -S $sample $id -Oz -o $wd/$idd.${sample%.*}.$idd_suffix
"
    done
done

wd="2.vcf_stats/1.grp"
mkdir -p $wd

for sample in hol.sample pub.sample jer.sample; do
    ls 2.vcf_stats/0.filter/jerMulHol-2024-08-13*filtered.vcf.gz  |
    while read id; do
        idd=$(basename $id);
        idd_suffix=${idd#*.};
        idd=${idd%%.*};
        sbatch -A ${SLURM_ACCOUNT} -J $idd.grp -c 8 --wrap="
bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -S $sample $id -Oz -o $wd/$idd.${sample%.*}.$idd_suffix
"
    done
done

## 2.2 vcf ind stats

wd="2.vcf_stats/2.ind"
mkdir -p $wd

cat all.sample | while read sample;do 
    ls 2.vcf_stats/0.filter/*.*.filtered.vcf.gz | 
        while read id; do
            idd=$(basename $id);
            idd_suffix=${idd#*.};
            idd=${idd%%.*};
            sbatch -A ${SLURM_ACCOUNT} -J $idd.ind -c 8 --wrap="
bcftools view --threads 8 -i 'F_MISSING<0.2' -c 1 -s $sample $id -Oz -o $wd/${idd}.${sample}.$idd_suffix
    "
        done
    done


### stat svtype

rm x.anno.stats/3.pan.filter.vcf.stats
ls 2.vcf_stats/0.filter/*.anno_biallelic.filtered.vcf.gz | grep .filtered.vcf.gz |
    while read id; do
        idd=$(basename $id);
        ale=${idd#*.};
        ale=${ale%.filtered.vcf.gz};
        idd=${idd%%-*};
        sbatch -A ${SLURM_ACCOUNT} -J f.sta.$idd -c 8 --wrap="
zcat $id | bash vcf.svtype.sh > ${id/vcf.gz/svtype}
cat ${id/vcf.gz/svtype} | bash count.svtype.sh $idd.$ale >> x.anno.stats/3.pan.filter.vcf.stats
cat ${id/vcf.gz/svtype} | bash max.svtype.sh $idd.$ale >> x.anno.stats/3.pan.filter.vcf.stats
"
    done

rm x.anno.stats/4.grp.filter.vcf.stats
ls 2.vcf_stats/1.grp/*.anno*vcf.gz |
        while read id; do
        idd=$(basename $id);
        ale=${idd#*.};
        ale=${ale%.filtered.vcf.gz};
        idd=${idd%%-*};
        sbatch -A ${SLURM_ACCOUNT} -J f.sta.$idd -c 8 --wrap="
zcat $id | bash vcf.svtype.sh > ${id/vcf.gz/svtype}
cat ${id/vcf.gz/svtype} | bash count.svtype.sh $idd.$ale >> x.anno.stats/4.grp.filter.vcf.stats
cat ${id/vcf.gz/svtype} | bash max.svtype.sh $idd.$ale >> x.anno.stats/4.grp.filter.vcf.stats
"
    done

rm x.anno.stats/5.ind.filter.vcf.stats
ls 2.vcf_stats/2.ind/*.anno*vcf.gz |
        while read id; do
        idd=$(basename $id);
        ale=${idd#*.};
        ale=${ale%.filtered.vcf.gz};
        idd=${idd%%-*};
        sbatch -A ${SLURM_ACCOUNT} -J f.sta.$idd -c 8 --wrap="
zcat $id | bash vcf.svtype.sh > ${id/vcf.gz/svtype}
cat ${id/vcf.gz/svtype} | bash count.svtype.sh $idd.$ale >> x.anno.stats/5.ind.filter.vcf.stats
cat ${id/vcf.gz/svtype} | bash max.svtype.sh $idd.$ale >> x.anno.stats/5.ind.filter.vcf.stats
"
    done

### stat allele 
ls 2.vcf_stats/*/*.*.filtered.vcf.gz |
    while read id; do
        sbatch -A ${SLURM_ACCOUNT} -c 2 --wrap="
bcftools stats --threads 2 $id > ${id/vcf.gz/vcfstats}
"
    done

cat > x.anno.stats/6.allele.vcf.stats <<EOF
ID samples records no-ALTs SNPs MNPs indels others multiallelic-sites multiallelic-SNP
EOF

ls 2.vcf_stats/*/*.*.filtered.vcfstats |
    while read id; do
        idd=$(basename $id);
        idd=${idd%.filtered.vcfstats};
        awk '$1 == "SN"{a[NR]=$NF} 
        END {
            printf "%s ",ID
            for(i in a)printf "%s ",a[i]
            print ""
        }' ID=$idd $id >> x.anno.stats/6.allele.vcf.stats
    done

#### stat allele auto.

ls 2.vcf_stats/*/*.*.filtered.vcf.gz | 
    while read id; do
        sbatch -A ${SLURM_ACCOUNT} -c 2 --wrap="
zcat $id | awk '\$1 ~ /#/ || \$1 <= 29 ' | bcftools stats > ${id/vcf.gz/vcfautostats}
"
    done

cat > x.anno.stats/7.allele.autovcf.stats <<EOF
ID samples records no-ALTs SNPs MNPs indels others multiallelic-sites multiallelic-SNP
EOF

ls 2.vcf_stats/*/*.*.filtered.vcfautostats |
    while read id; do
        idd=$(basename $id);
        idd=${idd%.filtered.vcfautostats};
        awk '$1 == "SN"{a[NR]=$NF} 
        END {
            printf "%s ",ID
            for(i in a)printf "%s ",a[i]
            print ""
        }' ID=$idd $id >> x.anno.stats/7.allele.autovcf.stats 
    done


#### stat allele X.

ls 2.vcf_stats/*/*.*.filtered.vcf.gz | 
    while read id; do
        sbatch -A ${SLURM_ACCOUNT} -c 2 --wrap="
zcat $id | awk '\$1 ~ /#/ || \$1 <= 29 ' | bcftools stats > ${id/vcf.gz/vcfautostats}
"
    done

cat > x.anno.stats/7.allele.autovcf.stats <<EOF
ID samples records no-ALTs SNPs MNPs indels others multiallelic-sites multiallelic-SNP
EOF

ls 2.vcf_stats/*/*.*.filtered.vcfautostats |
    while read id; do
        idd=$(basename $id);
        idd=${idd%.filtered.vcfautostats};
        awk '$1 == "SN"{a[NR]=$NF} 
        END {
            printf "%s ",ID
            for(i in a)printf "%s ",a[i]
            print ""
        }' ID=$idd $id >> x.anno.stats/7.allele.autovcf.stats 
    done

###
# ID samples records no-ALTs SNPs MNPs indels others multiallelic-sites multiallelic-SNP
cat 2.vcf_stats/vcf.filtered.stats.tab | grep sample_4232
bovHol-2024-08-13_sample_4232.filtered.vcf 1 7393872 0 6279662 406367 1331751 221751 1313106 220873
bovinePanPhase-2024-07-03_sample_4232.filtered.vcf 1 7351652 0 6241777 380540 1306328 203047 1251295 208799
hol-pg2hic-2024-05-22_sample_4232.filtered.vcf 1 7747253 0 6482229 242282 1128333 80129 534066 57537

cat 2.vcf_stats/vcf.filtered.stats.tab | grep Jersey_441
bovHol-2024-08-13_Jersey_441.filtered.vcf 1 7097704 0 5988366 386535 1323668 215390 1268052 209956
bovinePanPhase-2024-07-03_Jersey_441.filtered.vcf 1 7064850 0 5958040 363901 1299790 197979 1210532 198459
jer-pg-2024-05-15_Jersey_441.filtered.vcf 1 7281863 0 6050232 218376 1068610 59733 382568 35721

cat 2.vcf_stats/vcf.filter-bi.stats.tab | grep sample_4232
bovHol-2024-08-13_sample_4232 1 7543642 0 6119415 199980 1138730 85517 0 0
bovinePanPhase-2024-07-03_sample_4232 1 7498346 0 6096699 195150 1125286 81211 0 0
hol-pg2hic-2024-05-22_sample_4232 1 7885014 0 6444363 198724 1185062 56865 0 0

cat 2.vcf_stats/vcf.filter-bi.stats.tab | grep Jersey_441
bovHol-2024-08-13_Jersey_441 1 7230833 0 5830686 187519 1130989 81639 0 0
bovinePanPhase-2024-07-03_Jersey_441 1 7194905 0 5814587 183503 1118955 77860 0 0
jer-pg-2024-05-15_Jersey_441 1 7398756 0 6028346 196044 1126220 48146 0 0


























######### others #####
cd ${REF_DIR}/bostauFaDownload




!!! GCA_030378505