mkdir -p ${PROJECT_ROOT}/stat_pan/4.truvari_jer
cd ${PROJECT_ROOT}/stat_pan/4.truvari_jer
mkdir -p logs

chrs="$(echo X | awk '{for(i=1;i<=29;i++)printf i","; print}')"


cat > pangenie_remove50.awk <<'EOF'
BEGIN{
FS=OFS="\t"
}
{
    if($1 ~ "##"){
    print
    } else if($1 ~ "#" && $1 !~ "##"){
        print "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
        print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"
        print "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"
        print;
    } else {
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
            typ="MNP"
            leng=altlen-1
        } 
        #print $1,$2,$2+leng,leng,typ,"-",typ"-"NR 
        
        if(leng > 50){
        $8="SVTYPE="typ";END="$2+leng";SVLEN="leng";"$8
            print
        }
    }
}
EOF

mkdir -p 0.tmp 0.sv_vcfgs

ls ${PROJECT_ROOT}/stat_pan/2.vcf_stats/2.ind/jerHap-pg-2024-12-18.jer_*.anno_biallelic.filtered.vcf.gz | 
    while read id; do
        idd=$(basename $id | sed -e 's/jerHap-pg-2024-12-18.//' -e 's/\.anno_biallelic.filtered.vcf.gz//')
        sbatch -A ${SLURM_ACCOUNT} -J $idd.panBed \
            -o logs/$idd.panBed.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
zcat $id | bcftools sort -Oz -o 0.tmp/$idd.vcf.gz && tabix 0.tmp/$idd.vcf.gz
chrs="$(echo X | awk '{for(i=1;i<=29;i++)printf i","; print}')"
bcftools view --regions $chrs 0.tmp/$idd.vcf.gz | awk -f pangenie_remove50.awk | bcftools sort -Oz -o 0.sv_vcfgs/$idd.pan.vcf.gz && tabix 0.sv_vcfgs/$idd.pan.vcf.gz
"
    done

## mc pri
ls ${PROJECT_ROOT}/stat_pan/2.vcf_stats/2.ind/jerPri-pg-2024-12-18.j*.anno_biallelic.filtered.vcf.gz |  
    while read id; do
        idd=$(basename $id | sed -e 's/jerPri-pg-2024-12-18.//' -e 's/\.anno_biallelic.filtered.vcf.gz//')
        sbatch -A ${SLURM_ACCOUNT} -J $idd.mcPri \
            -o logs/$idd.mcPri.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
zcat $id | bcftools sort -Oz -o 0.tmp/$idd.vcf.gz && tabix 0.tmp/$idd.vcf.gz
chrs="$(echo X | awk '{for(i=1;i<=29;i++)printf i","; print}')"
bcftools view --regions $chrs 0.tmp/$idd.vcf.gz | awk -f pangenie_remove50.awk | bcftools sort -Oz -o 0.sv_vcfgs/$idd.mcPri.vcf.gz && tabix 0.sv_vcfgs/$idd.mcPri.vcf.gz
"
    done

# long sSVs
chrs="$(echo X | awk '{for(i=1;i<=29;i++)printf i","; print}')"
ls ${PANEL_DIR}/jer_*/4.sv_vcf/jer_*.vcf.gz |
    while read id; do
        idd=$(basename $id | sed -e 's/.vcf.gz//')
        sbatch -A ${SLURM_ACCOUNT} -J $idd \
            -o logs/$idd.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
bcftools view --regions $chrs -i '(INFO/SVLEN >50 && INFO/SVLEN < 1e6) || (INFO/SVLEN < -50 && INFO/SVLEN >- 1e6)' $id | 
    bcftools sort -Oz -o 0.sv_vcfgs/$idd.vcf.gz && tabix 0.sv_vcfgs/$idd.vcf.gz
"
    done


svtoolsMerged_vcf=${PROJECT_ROOT}/pangenie_HiFi/5.svtools_all/output.ls.filter.vcf.gz
bcftools view -h $svtoolsMerged_vcf > hdr.txt
sed -i 's|##bcftools_viewVersion|##FORMAT=<ID=CN,Number=A,Type=Float,Description="Copy number">\n##bcftools_viewVersion|' hdr.txt
cat ${PROJECT_ROOT}/stat_pan/jer.sample |
    while read id; do
        idd=$id
        sbatch -A ${SLURM_ACCOUNT} -J $idd.svtools \
            -o logs/$idd.svtools.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
bcftools reheader -h hdr.txt $svtoolsMerged_vcf | bcftools view -i 'F_MISSING<0.2' -s ${id/jer_/} | bcftools sort -Oz -o 0.tmp/$idd.svtools-sr.vcf.gz && tabix 0.tmp/$idd.svtools-sr.vcf.gz
bcftools view --regions $chrs 0.tmp/$idd.svtools-sr.vcf.gz -Oz -o 0.sv_vcfgs/$idd.svtools-sr.vcf.gz && tabix 0.sv_vcfgs/$idd.svtools-sr.vcf.gz
"
    done

pangenieMerged_vcf=${PROJECT_ROOT}/pangenie_HiFi/3.pan_all/jerHap-pg-2024-12-18_graph_genotyping.merge-biallelic.filter.vcf.gz
cat ${PROJECT_ROOT}/stat_pan/jer.sample |
    while read id; do
        idd=$id
        sbatch -A ${SLURM_ACCOUNT} -J $idd.pangenie \
            -o logs/$idd.pangenie.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
        
bcftools view -s ${id/jer_/} --regions $chrs $pangenieMerged_vcf | awk -f pangenie_remove50.awk  | bcftools sort -Oz -o 0.sv_vcfgs/$idd.pangenie-sr.vcf.gz && tabix 0.sv_vcfgs/$idd.pangenie-sr.vcf.gz
"
    done

pangenieMerged_vcf2=${PROJECT_ROOT}/pangenie_HiFi/3.pan_all/hol-pg2hic-2024-05-22_graph_genotyping.merge-biallelic.filter.vcf.gz
cat ${PROJECT_ROOT}/stat_pan/jer.sample |
    while read id; do
        idd=$id
        sbatch -A ${SLURM_ACCOUNT} -J $idd.pangenie-jerForHol \
            -o logs/$idd.pangenie-jerForHol.out \
            --cpus-per-task=1 \
            --mem-per-cpu=8g \
            --wrap="
chrs="$(echo X | awk '{for(i=1;i<=29;i++)printf i","; print}')"        
bcftools view -s ${id/jer_/} --regions $chrs $pangenieMerged_vcf2 | awk -f pangenie_remove50.awk  | bcftools sort -Oz -o 0.sv_vcfgs/$idd.pangenie-jerForHol.vcf.gz && tabix 0.sv_vcfgs/$idd.pangenie-jerForHol.vcf.gz
"
    done




l 0.sv_vcfgs/*gz | awk '$5 ~ "M"{print $9}'  | while read id ; do basename $id | cut -d. -f 2 | sort ;done | sort | uniq -c

####################################################
# 1.6.sv.com.sh

ref_path=${REF_DIR}
ref_fa=$ref_path/ARS_UCD_v2.0.fa

tools="pav_bp pav_diploid sv-cutesv sv-pbsv sv-sniffles sv-svim sv-svim-asm-bp sv-svim-asm-diploid sv-svision svtools-sr pangenie-sr pangenie-holForJer mcPri"

#rm 1.truvari/
mkdir -p 1.truvari

echo "ID tool cla precision recall f1 gt_concordance" | sed 's/ /\t/g' > truvari.summary.jer.tsv

sa truvari
cat ${PROJECT_ROOT}/stat_pan/jer.sample | 
    while read id; do
        for tool in $tools; do
            idd=$id
            sbatch -A ${SLURM_ACCOUNT} -J $idd.$tool \
                -o logs/$idd.$tool.out \
                --cpus-per-task=1 \
                --mem-per-cpu=8g \
                --wrap="
truvari bench -b 0.sv_vcfgs/$id.pan.vcf.gz -c  0.sv_vcfgs/$id.$tool.vcf.gz \
    --dup-to-ins -o 1.truvari/$id.$tool.truvari -r 2000 -C 2000 \
    --no-ref a --passonly --sizemax 1000000 --sizemin 50
# --passonly
cat 1.truvari/$id.$tool.truvari/summary.json  |  
    ${CONDA_BASE}/bin/jq -r '[.precision, .recall, .f1, .gt_concordance] | @tsv' |
    sed \"s|^|$id\t$tool\tALL\t|\" >> truvari.summary.jer.tsv

rm 1.truvari/$id.$tool.truvari -rf
                "
        done
    done

### different regions
tools="pav_bp pav_diploid sv-cutesv sv-pbsv sv-sniffles sv-svim sv-svim-asm-bp sv-svim-asm-diploid sv-svision svtools-sr pangenie-sr pangenie-holForJer mcPri"


echo "ID tool cla precision recall f1 gt_concordance" | sed 's/ /\t/g' > truvari.regions.jer.tsv

ref_rm=${REF_DIR}/ARS_UCD_v2.0.ref_repeat
rms="None RM DNA LTR Low_complexity LINE srpRNA rRNA Unknown RNA RC scRNA SINE Satellite Simple_repeat tRNA snRNA"

for rm in $rms; do
cat ${PROJECT_ROOT}/stat_pan/jer.sample | 
    while read id; do
        for tool in $tools; do
            sbatch -A ${SLURM_ACCOUNT} -J $id.$tool.$rm \
                -o logs/$id.$tool.$rm.out \
                --cpus-per-task=1 \
                --mem-per-cpu=8g \
                --wrap="
truvari bench \
    -b 0.sv_vcfgs/$id.pan.vcf.gz \
    -c 0.sv_vcfgs/$id.$tool.vcf.gz \
    --dup-to-ins --no-ref a --passonly \
    -r 2000 -C 2000 --sizemax 1000000 --sizemin 50 \
    -o 1.truvari/$id.$tool.truvari.$rm \
    --includebed $ref_rm/rm.$rm.bed

cat 1.truvari/$id.$tool.truvari.$rm/summary.json  |  
    ${CONDA_BASE}/bin/jq -r '[.precision, .recall, .f1, .gt_concordance] | @tsv' |
    sed \"s|^|$id\t$tool\t$rm\t|\" >> truvari.regions.jer.tsv

rm 1.truvari/$id.$tool.truvari.$rm -rf
                "
        done
    done
done

### different size

echo "ID tool cla precision recall f1 gt_concordance" | sed 's/ /\t/g' > truvari.size.jer.tsv

tools="pav_bp pav_diploid sv-cutesv sv-pbsv sv-sniffles sv-svim sv-svim-asm-bp sv-svim-asm-diploid sv-svision svtools-sr pangenie-sr pangenie-holForJer mcPri"

len=(50 200 500 1000 10000 100000 1000000)

cat ${PROJECT_ROOT}/stat_pan/jer.sample | 
    while read id; do
        idd=$id
        seq 1 6 | while read i; do
            tool="pan"
            sbatch -A ${SLURM_ACCOUNT} -J $idd.$tool.$i \
                    -o logs/$idd.$tool.$i.out \
                    --cpus-per-task=1 \
                    --mem-per-cpu=8g \
                    --wrap="
bcftools view -i \"(INFO/SVLEN > -${len[$i]} & INFO/SVLEN < -${len[$i-1]}) | (INFO/SVLEN < ${len[$i]} & INFO/SVLEN > ${len[$i-1]})\" 0.sv_vcfgs/$id.$tool.vcf.gz -Oz -o 0.sv_vcfgs/$id.$tool.size_$i.vcf.gz && tabix 0.sv_vcfgs/$id.$tool.size_$i.vcf.gz
    "
        done
    done

cat ${PROJECT_ROOT}/stat_pan/jer.sample | 
    while read id; do
        idd=$id
        seq 1 6 | while read i; do
            for tool in $tools; do
                idd=$id
                sbatch -A ${SLURM_ACCOUNT} -J $idd.$tool.$i \
                    -o logs/$idd.$tool.$i.out \
                    --cpus-per-task=1 \
                    --mem-per-cpu=8g \
                    --wrap="
bcftools view -i \"(INFO/SVLEN > -${len[$i]} & INFO/SVLEN < -${len[$i-1]}) | (INFO/SVLEN < ${len[$i]} & INFO/SVLEN > ${len[$i-1]})\" 0.sv_vcfgs/$id.$tool.vcf.gz -Oz -o 0.sv_vcfgs/$id.$tool.size_$i.vcf.gz  && tabix 0.sv_vcfgs/$id.$tool.size_$i.vcf.gz

truvari bench -b 0.sv_vcfgs/$id.pan.size_$i.vcf.gz -c0.sv_vcfgs/$id.$tool.size_$i.vcf.gz \
    --dup-to-ins -o0.sv_vcfgs/$id.$tool.size_$i -r 2000 -C 2000 \
    --no-ref a --passonly --sizemax 1000000 --sizemin 50

cat0.sv_vcfgs/$id.$tool.size_$i/summary.json  |  
    ${CONDA_BASE}/bin/jq -r '[.precision, .recall, .f1, .gt_concordance] | @tsv' |
    sed \"s|^|$id\t$tool\t$i\t|\" >> truvari.size.jer.tsv

# rm0.sv_vcfgs/$id.$tool.size_$i -rf
            "
        done
    done
    done

### different type

mkdir 1.truvari

echo "ID tool cla precision recall f1 gt_concordance" | sed 's/ /\t/g' > truvari.type.jer.tsv

tools="pav_bp pav_diploid sv-cutesv sv-pbsv sv-sniffles sv-svim sv-svim-asm-bp sv-svim-asm-diploid sv-svision svtools-sr pangenie-sr pangenie-holForJer mcPri"

cat ${PROJECT_ROOT}/stat_pan/jer.sample | 
    while read id; do
        idd=$id
        tool="pan"
        sbatch -A ${SLURM_ACCOUNT} -J $idd.$tool.type \
                -o logs/$idd.$tool.type.out \
                --cpus-per-task=1 \
                --mem-per-cpu=8g \
                --wrap="
bcftools view -i \"INFO/SVTYPE == 'DEL'\" 0.sv_vcfgs/$id.$tool.vcf.gz -Oz -o 0.sv_vcfgs/$id.$tool.DEL.vcf.gz  && tabix 0.sv_vcfgs/$id.$tool.DEL.vcf.gz
bcftools view -i \"INFO/SVTYPE == 'INS'\" 0.sv_vcfgs/$id.$tool.vcf.gz -Oz -o 0.sv_vcfgs/$id.$tool.INS.vcf.gz  && tabix 0.sv_vcfgs/$id.$tool.INS.vcf.gz
    "
    done

cat ${PROJECT_ROOT}/stat_pan/jer.sample | 
    while read id; do
        idd=$id
            for tool in $tools; do
                idd=$id
                sbatch -A ${SLURM_ACCOUNT} -J $idd.$tool.type \
                    -o logs/$idd.$tool.type.out \
                    --cpus-per-task=1 \
                    --mem-per-cpu=8g \
                    --wrap="
bcftools view -i \"INFO/SVTYPE == 'DEL'\" 0.sv_vcfgs/$id.$tool.vcf.gz -Oz -o 0.sv_vcfgs/$id.$tool.DEL.vcf.gz  && tabix 0.sv_vcfgs/$id.$tool.DEL.vcf.gz
bcftools view -i \"INFO/SVTYPE == 'INS' | INFO/SVTYPE == 'DUP'\" 0.sv_vcfgs/$id.$tool.vcf.gz -Oz -o 0.sv_vcfgs/$id.$tool.INS.vcf.gz  && tabix 0.sv_vcfgs/$id.$tool.INS.vcf.gz

            "
        done
    done
    done

cat ${PROJECT_ROOT}/stat_pan/jer.sample | 
    while read id; do
        idd=$id
        for type in DEL INS; do
            for tool in $tools; do
                idd=$id
                sbatch -A ${SLURM_ACCOUNT} -J $idd.$tool.type \
                    -o logs/$idd.$tool.type.out \
                    --cpus-per-task=1 \
                    --mem-per-cpu=8g \
                    --wrap="
truvari bench -b 1.truvari/$id.pan.$type.vcf.gz -c 1.truvari/$id.$tool.$type.vcf.gz \
    --dup-to-ins -o 1.truvari/$id.$tool.$type -r 2000 -C 2000 \
    --no-ref a --passonly --sizemax 1000000 --sizemin 50

cat 1.truvari/$id.$tool.$type/summary.json  |  
    ${CONDA_BASE}/bin/jq -r '[.precision, .recall, .f1, .gt_concordance] | @tsv' |
    sed \"s|^|$id\t$tool\t$type\t|\" >> truvari.type.jer.tsv
done
            "
        done
    done
    done


scp -r ${CLUSTER_HOST}:${PROJECT_ROOT}/stat_pan/4.truvari_jer/*tsv .

C:\Users\${USER}\OneDrive - University of Maryland\Data\2024-02-07.cattleLR-SR-GWAS\2.analyses\1.6.sv_truvari



#### truvari for pan sv

ref_path=${REF_DIR}
ref_fa=$ref_path/ARS_UCD_v2.0.fa

tools1="pangenie-sr pangenie-jerForHol mcPri"
tools2="pav_bp pav_diploid sv-cutesv sv-pbsv sv-sniffles sv-svim sv-svim-asm-bp sv-svim-asm-diploid sv-svision pangenie-sr pangenie-jerForHol mcPri"

#rm 1.truvari/
mkdir -p 1.truvari

summtsv=truvari.pangenie-jer.tsv
echo "ID tool cla precision recall f1 gt_concordance" | sed 's/ /\t/g' > $summtsv

sa truvari
cat ${PROJECT_ROOT}/stat_pan/jer.sample | 
    while read id; do
        for tool1 in $tools1; do
        for tool2 in $tools2; do
        if [ $tool1 == $tool2 ]; then continue; fi
            sbatch -A ${SLURM_ACCOUNT} -J $id.$tool1.truvari \
                -o logs/$id.$tool1.truvari.out \
                --cpus-per-task=1 \
                --mem-per-cpu=8g \
                --wrap="
if [ -f 1.truvari/$id.$tool2.$tool1.truvari-pan ]; then rm -rf 1.truvari/$id.$tool2.$tool1.truvari-pan; fi
truvari bench -b 0.sv_vcfgs/$id.$tool1.vcf.gz -c  0.sv_vcfgs/$id.$tool2.vcf.gz \
    --dup-to-ins -o 1.truvari/$id.$tool2.$tool1.truvari-pan -r 2000 -C 2000 \
    --no-ref a --passonly --sizemax 1000000 --sizemin 50

cat 1.truvari/$id.$tool2.$tool1.truvari-pan/summary.json  |  
    ${CONDA_BASE}/bin/jq -r '[.precision, .recall, .f1, .gt_concordance] | @tsv' |
    sed \"s|^|$id\t${tool2}.${tool1}\tALL\t|\" >> $summtsv

rm 1.truvari/$id.$tool2.$tool1.truvari -rf
                "
            sleep 0.2
        done
    done
done
