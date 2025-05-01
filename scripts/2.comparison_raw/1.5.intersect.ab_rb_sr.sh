
mkdir -p ${PROJECT_ROOT}/stat_pan/3.intersect.ab_rb_sr
cd ${PROJECT_ROOT}/stat_pan/3.intersect.ab_rb_sr
mkdir -p 0.panBed 1.ab_rb_srBed

# 0.panBed
cat > panBed.sh <<'EOF'
bcftools query -f "%CHROM %POS %ID %INFO/ID\n" |
    awk '{split($4,bc,"-");print $1,$2,$2+bc[5],bc[5],bc[3],$3,bc[4]}' |
    awk '($1 <= 29 || $1 == "X") && $4 > 50 && $4 <= 1000000' |
    sed 's/ /\t/g'
EOF

ls ${PROJECT_ROOT}/stat_pan/2.vcf_stats/2.ind/hol-pg2hic-2024-05-22.sample_*.anno_biallelic.filtered.vcf.gz | 
    grep -v "4439" | 
    while read id; do
        idd=$(basename $id | sed -e 's/hol-pg2hic-2024-05-22.//' -e 's/\.anno_biallelic.filtered.vcf.gz//')
        sbatch -A ${SLURM_ACCOUNT} -J $idd.panBed \
            -o $idd.panBed.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
zcat $id | bash panBed.sh > 0.panBed/$idd.bed
"
    done


ls ${PROJECT_ROOT}/stat_pan/2.vcf_stats/2.ind/holPri-2024-12-03.sample_*.anno_biallelic.filtered.vcf.gz | 
    grep -v "4439" | 
    while read id; do
        idd=$(basename $id | sed -e 's/holPri-2024-12-03.//' -e 's/\.anno_biallelic.filtered.vcf.gz//')
        sbatch -A ${SLURM_ACCOUNT} -J $idd.mcPri \
            -o $idd.mcPri.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
zcat $id | bash panBed.sh > 1.ab_rb_srBed/$idd.mcPri.bed
"
    done

# 1.cutesv
cat > ab_rbBed.sh <<'EOF'
bcftools query -f "%CHROM %POS %ID %INFO/SVTYPE %INFO/SVLEN\n" | 
    awk '{if($5 == "."){$5=1}else if($5 < 0){$5=-$5};print $1,$2,$2+$5,$5,$4,"-",$3}' |
    awk '($1 <= 29 || $1 == "X") && $4 > 50 && $4 <= 1000000' |
    sed 's/ /\t/g'
EOF

ls ${PANEL_DIR}/sample_*/4.sv_vcf/sample_*.vcf.gz |
    while read id; do
        idd=$(basename $id | sed -e 's/.vcf.gz//')
        sbatch -A ${SLURM_ACCOUNT} -J $idd \
            -o $idd.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
zcat $id | bash ab_rbBed.sh > 1.ab_rb_srBed/$idd.bed
"
    done
#svision $7 is 0xxx


# 2.srSVtools
svtoolsMerged_vcf=${PROJECT_ROOT}/pangenie_HiFi/5.svtools_all/output.ls.filter.vcf.gz
bcftools view -h $svtoolsMerged_vcf > hdr.txt
sed -i 's|##bcftools_viewVersion|##FORMAT=<ID=CN,Number=A,Type=Float,Description="Copy number">\n##bcftools_viewVersion|' hdr.txt
cat ${PROJECT_ROOT}/stat_pan/hol.sample |
    while read id; do
        idd=$id
        sbatch -A ${SLURM_ACCOUNT} -J $idd \
            -o $idd.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
bcftools reheader -h hdr.txt $svtoolsMerged_vcf | bcftools view -i 'F_MISSING<0.2' | 
    bcftools view -s ${id/sample_/} -c 1 | bash ab_rbBed.sh > 1.ab_rb_srBed/$idd.svtools-sr.bed
"
    done

# lumpy
ls ${PROJECT_ROOT}/pangenie_HiFi/*/5.svtools_data/2.lumpy_vcf/*.part1-smoove.genotyped.vcf.gz |
    while read id; do
        idd=sample_$(basename $id | sed -e 's/\..*//')
        sbatch -A ${SLURM_ACCOUNT} -J $idd \
            -o $idd.lumpy-sr.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
    zcat $id | bash ab_rbBed.sh > 1.ab_rb_srBed/$idd.lumpy-sr.bed
"
    done

# cnvnator 
cat > cnvnator.bed.covert.awk <<'EOF'
BEGIN{OFS="\t"}
{
    print $1,$2,$3,$4,toupper(substr($5,1,3)),"-",NR
}
EOF

ls ${PROJECT_ROOT}/pangenie_HiFi/*/5.svtools_data/5.cnvnator_rawcnv/*.cnvnator.bed |
    while read id; do
        idd=sample_$(basename $id | sed -e 's/\..*//')
        sbatch -A ${SLURM_ACCOUNT} -J $idd \
            -o $idd.cnvnator-sr.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
    zcat $id | awk -f cnvnator.bed.covert.awk > 1.ab_rb_srBed/$idd.cnvnator-sr.bed
"
    done


# pangenie

cat > pangenie-sr_bed.awk <<'EOF'
BEGIN{
FS=OFS="\t"
}
{
    altlen=length($5);
    reflen=length($4);
    leng=0;
    if(altlen > reflen && reflen == 1){
        typ="INS"
        leng=altlen-1
    } else if(altlen < reflen && altlen == 1){
        typ="DEL"
        leng=reflen-1
    } else if(altlen == reflen && reflen == 1){
        typ="SNP"
        leng=1
    } else {
        typ="COMPLEX"
        leng=altlen-1
    } 
    if(($1 <= 29 || $1 == "X") && leng > 50 && leng <= 1000000){
        print $1,$2,$2+leng,leng,typ,"-",typ"-"NR
    }
}
EOF

pangenieMerged_vcf=${PROJECT_ROOT}/pangenie_HiFi/3.pan_all/hol-pg2hic-2024-05-22_graph_genotyping.merge-biallelic.filter.vcf.gz
cat ${PROJECT_ROOT}/stat_pan/hol.sample |
    while read id; do
        idd=$id
        sbatch -A ${SLURM_ACCOUNT} -J $idd \
            -o $idd.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
bcftools view -s ${id/sample_/} -c 1 $pangenieMerged_vcf  | 
    grep -v '#' | awk -f pangenie-sr_bed.awk > 1.ab_rb_srBed/$idd.pangenie-sr.bed
"
    done


#### pangenie2

pangenieMerged_vcf2=${PROJECT_ROOT}/pangenie_HiFi/3.pan_all/jerHap-pg-2024-12-18_graph_genotyping.merge-biallelic.filter.vcf.gz
cat ${PROJECT_ROOT}/stat_pan/hol.sample |
    while read id; do
        idd=$id
        sbatch -A ${SLURM_ACCOUNT} -J $idd \
            -o $idd.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
bcftools view -s ${id/sample_/} -c 1 $pangenieMerged_vcf2  | 
    grep -v '#' | awk -f pangenie-sr_bed.awk > 1.ab_rb_srBed/$idd.pangenie-holForJer.bed
"
    done

#### intersect for pan sv
mkdir -p 2.sv.intersect_both

#### stats

mkdir -p 3.stats
mkdir -p 4.share_stats
cat > pan_stats.awk <<'EOF'
#!/bin/awk -f
    function count_num(bed,chr,type,lenCate,extp)
    {
        array["count",bed,"ALL","ALL",lenCate,extp]++;
        array["count",bed,chr,"ALL",lenCate,extp]++;
        array["count",bed,"ALL",type,lenCate,extp]++;
        array["count",bed,chr,type,lenCate,extp]++;
    }

    function count_len(bed,chr,type,lenCate,extp,shareLen)
    {
        array["length",bed,"ALL","ALL",lenCate,extp]+=shareLen;
        array["length",bed,chr,"ALL",lenCate,extp]+=shareLen;
        array["length",bed,"ALL",type,lenCate,extp]+=shareLen;
        array["length",bed,chr,type,lenCate,extp]+=shareLen;
    }

    function count_point(bed,chr,type,lenCate,extp,shareLen)
    {
        count_num(bed,chr,type,lenCate,"no");
        count_len(bed,chr,type,lenCate,"no",shareLen);
    }

    function each_count_point(bed,chr,type,lenCate,extp,shareLen)
    {
        count_point(bed,chr,type,"ALL",extp,shareLen);
        split(lenCate,alenCate_each,",")
        for(i=1;i<=length(alenCate_each);i++){
            if  (shareLen <= alenCate_each[i]){
                count_point(bed,chr,type,alenCate_each[i],extp,shareLen);
            }
        } 
    }

    BEGIN {
        lenCate="200,1000,10000,100000,1000000"
    }

    {   
        bed=""
        extp=""
        chr=$1
        type=$5
        shareLen=$4
        each_count_point(bed,chr,type,lenCate,extp,shareLen)
    }
        
    END {
        for (comb in array) {
            #split(comb,sep,SUBSEP);
            split(comb,a,SUBSEP)
            for(i=1;i<=length(a);i++){
                printf a[i]"\t"
            }
            print array[comb]
        }
    }
EOF

cat > pan_share_stats.awk <<'EOF'
#~/bin/awk -f
    function count_num(bed,chr,type,lenCate,extp)
    {
        array["count",bed,"ALL","ALL",lenCate,extp]++;
        array["count",bed,chr,"ALL",lenCate,extp]++;
        array["count",bed,"ALL",type,lenCate,extp]++;
        array["count",bed,chr,type,lenCate,extp]++;
    }

    function count_len(bed,chr,type,lenCate,extp,shareLen)
    {
        array["length",bed,"ALL","ALL",lenCate,extp]+=shareLen;
        array["length",bed,chr,"ALL",lenCate,extp]+=shareLen;
        array["length",bed,"ALL",type,lenCate,extp]+=shareLen;
        array["length",bed,chr,type,lenCate,extp]+=shareLen;
    }

    function count_point(bed,chr,type,lenCate,extp,shareLen)
    {
        count_num(bed,chr,type,lenCate,"no");
        count_len(bed,chr,type,lenCate,"no",shareLen);
        if(extp){
            count_num(bed,chr,type,lenCate,"yes");
            count_len(bed,chr,type,lenCate,"yes",shareLen);
        }
    }

    function each_count_point(bed,chr,type,lenCate,extp,shareLen)
    {
        count_point(bed,chr,type,"ALL",extp,shareLen);
        split(lenCate,alenCate_each,",")
        for(i=1;i<=length(alenCate_each);i++){
            if  (shareLen <= alenCate_each[i]){
                count_point(bed,chr,type,alenCate_each[i],extp,shareLen);
            }
        } 
    }
    
    BEGIN {
        lenCate="200,1000,10000,100000,1000000"
    }

    {
        if($5 != $12){$5=$5"-"$12}
        if($15 < 0){$15 = 1}
        uniq_pan[$7]++;
        uniq_tool[$14]++;
        
        extp=$2 == $9 && $3 == $10
        chr=$1
        type=$5
        shareLen=$15
        if(uniq_pan[$7] == 1){
            bed="pan"
            each_count_point(bed,chr,type,lenCate,extp,shareLen)
        }

        if(uniq_tool[$14] == 1){
            bed="tool"
            each_count_point(bed,chr,type,lenCate,extp,shareLen)
        }
    }

    END {
        for (comb in array) {
            #split(comb,sep,SUBSEP);
            split(comb,a,SUBSEP)
            for(i=1;i<=length(a);i++){
                printf a[i]"\t"
            }
            print array[comb]
        }
    }

EOF


tools="pav_bp pav_diploid sv-cutesv sv-pbsv sv-sniffles sv-svim sv-svim-asm-bp sv-svim-asm-diploid sv-svision svtools-sr lumpy-sr cnvnator-sr pangenie-sr pangenie-holForJer mcPri"
cat ${PROJECT_ROOT}/stat_pan/hol.sample | 
    while read id; do
        for tool in $tools; do
            idd=$id
            sbatch -A ${SLURM_ACCOUNT} -J $idd.$tool \
                -o $idd.$tool.out \
                --cpus-per-task=1 \
                --mem-per-cpu=8g \
                --wrap="
bedtools intersect \
    -a 0.panBed/$id.bed \
    -b 1.ab_rb_srBed/$id.$tool.bed \
    -f 0.90 -r -wo > 2.sv.intersect_both/$id.${tool}_pan.bed
                
cat  0.panBed/$id.bed | awk -f pan_stats.awk > 3.stats/$id.pan.stats
cat 1.ab_rb_srBed/$id.$tool.bed | awk -f pan_stats.awk > 3.stats/$id.$tool.stats
awk -f pan_share_stats.awk 2.sv.intersect_both/$id.${tool}_pan.bed > 4.share_stats/$id.$tool.sv.share.stats
                "
        done
    done

scp -r ${CLUSTER_HOST}:${PROJECT_ROOT}/stat_pan/3.intersect.ab_rb_sr/3.stats .
scp -r ${CLUSTER_HOST}:${PROJECT_ROOT}/stat_pan/3.intersect.ab_rb_sr/4.share_stats .


C:\Users\${USER}\OneDrive - University of Maryland\Data\2024-02-07.cattleLR-SR-GWAS\2.analyses\1.5.intersect




