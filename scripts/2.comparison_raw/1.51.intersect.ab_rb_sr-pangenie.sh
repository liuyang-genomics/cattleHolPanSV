# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
: "${REF_DIR:?REF_DIR is unset - see config.sh.example at the repository root}"
# --------------------------


mkdir -p ${PROJECT_ROOT}/stat_pan/3.intersect.ab_rb_sr
cd ${PROJECT_ROOT}/stat_pan/3.intersect.ab_rb_sr
mkdir -p 0.panBed 1.ab_rb_srBed

#### intersect for pan sv
mkdir -p 7.sv.intersect_pangenie
mkdir -p 9.sv.share_pangenie

tools1="pangenie-sr pangenie-holForJer mcPri svtools-sr lumpy-sr cnvnator-sr"
tools2="pav_bp pav_diploid sv-cutesv sv-pbsv sv-sniffles sv-svim sv-svim-asm-bp sv-svim-asm-diploid sv-svision svtools-sr lumpy-sr cnvnator-sr pangenie-sr pangenie-holForJer mcPri"
cat ${PROJECT_ROOT}/stat_pan/hol.sample | 
    while read id; do
        for tool1 in $tools1; do
        for tool2 in $tools2; do
        if [ $tool1 == $tool2 ]; then continue; fi
            idd=$id
            sbatch -A ${SLURM_ACCOUNT} -J $idd.$tool1.$tool2 \
                -o $idd.$tool1.$tool2.out \
                --cpus-per-task=1 \
                --mem-per-cpu=8g \
                --wrap="
bedtools intersect \
    -a 1.ab_rb_srBed/$id.$tool1.bed \
    -b 1.ab_rb_srBed/$id.$tool2.bed \
    -f 0.90 -r -wo > 7.sv.intersect_pangenie/$id.${tool2}_${tool1}.bed
awk -f pan_share_stats.awk 7.sv.intersect_pangenie/$id.${tool2}_$tool1.bed > 9.sv.share_pangenie/$id.$tool2.$tool1.share.stats
                "
            sleep 0.2
        done
    done
done

#### intersect for pan sv
mkdir -p 7.sv.intersect_pangenie
mkdir -p 10.sv.share_pangenie-jer

tools1="pangenie-sr pangenie-jerForHol mcPri svtools-sr lumpy-sr cnvnator-sr"
tools2="pav_bp pav_diploid sv-cutesv sv-pbsv sv-sniffles sv-svim sv-svim-asm-bp sv-svim-asm-diploid sv-svision svtools-sr lumpy-sr cnvnator-sr pangenie-sr pangenie-jerForHol mcPri"
cat ${PROJECT_ROOT}/stat_pan/jer.sample | 
    while read id; do
        for tool1 in $tools1; do
        for tool2 in $tools2; do
        if [ $tool1 == $tool2 ]; then continue; fi
            idd=$id
            sbatch -A ${SLURM_ACCOUNT} -J $idd.$tool1.$tool2 \
                -o $idd.$tool1.$tool2.out \
                --cpus-per-task=1 \
                --mem-per-cpu=8g \
                --wrap="
bedtools intersect \
    -a 1.ab_rb_srBed/$id.$tool1.bed \
    -b 1.ab_rb_srBed/$id.$tool2.bed \
    -f 0.90 -r -wo > 7.sv.intersect_pangenie/$id.${tool2}_${tool1}.bed
awk -f pan_share_stats.awk 7.sv.intersect_pangenie/$id.${tool2}_$tool1.bed > 10.sv.share_pangenie-jer/$id.$tool2.$tool1.share.stats
                "
            sleep 0.2
        done
    done
done


#### truvari for pan sv

ref_path=${REF_DIR}
ref_fa=$ref_path/ARS_UCD_v2.0.fa

tools1="pangenie-sr pangenie-holForJer mcPri"
tools2="pav_bp pav_diploid sv-cutesv sv-pbsv sv-sniffles sv-svim sv-svim-asm-bp sv-svim-asm-diploid sv-svision pangenie-sr pangenie-holForJer mcPri"


#rm 1.truvari/
mkdir -p 1.truvari

summtsv=truvari.pangenie.tsv
echo "ID tool cla precision recall f1 gt_concordance" | sed 's/ /\t/g' > $summtsv

sa truvari
cat ${PROJECT_ROOT}/stat_pan/hol.sample | 
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



### different size
summtsv=truvari.pangenie-regions.tsv

echo "ID tool cla precision recall f1 gt_concordance" | sed 's/ /\t/g' > $summtsv
ref_rm=${REF_DIR}/ARS_UCD_v2.0.ref_repeat
rms="None RM DNA LTR Low_complexity LINE srpRNA rRNA Unknown RNA RC scRNA SINE Satellite Simple_repeat tRNA snRNA"

for rm in $rms; do

cat ${PROJECT_ROOT}/stat_pan/hol.sample | 
    while read id; do
        for tool1 in $tools1; do
        for tool2 in $tools2; do
        if [ $tool1 == $tool2 ]; then continue; fi
            idd=$id
            sbatch -A ${SLURM_ACCOUNT} -J $idd.rm$rm.$tool1 \
                -o logs/$idd.rm$rm.$tool1.out \
                --cpus-per-task=1 \
                --mem-per-cpu=8g \
                --wrap="
if [ -f 1.truvari/$id.$tool2.$tool1.truvari-pan.$rm  ]; then 
    rm -rf 1.truvari/$id.$tool2.$tool1.truvari-pan.$rm ; 
fi
truvari bench -b 0.sv_vcfgs/$id.$tool1.vcf.gz -c  0.sv_vcfgs/$id.$tool2.vcf.gz \
    --dup-to-ins -o 1.truvari/$id.$tool2.$tool1.truvari-pan.$rm -r 2000 -C 2000 \
    --no-ref a --passonly --sizemax 1000000 --sizemin 50 \
    --includebed $ref_rm/rm.$rm.bed

cat 1.truvari/$id.$tool2.$tool1.truvari-pan.$rm/summary.json  |  
    ${CONDA_BASE}/bin/jq -r '[.precision, .recall, .f1, .gt_concordance] | @tsv' |
    sed \"s|^|$id\t${tool2}.${tool1}\t$rm\t|\" >> $summtsv

rm  1.truvari/$id.$tool2.$tool1.truvari.$rm -rf
                "
            sleep 0.2
        done
    done
done
done


# size
summtsv=truvari.pangenie-size.tsv
echo "ID tool cla precision recall f1 gt_concordance" | sed 's/ /\t/g' > $summtsv

len=(50 200 500 1000 10000 100000 1000000)

cat ${PROJECT_ROOT}/stat_pan/hol.sample | 
    while read id; do
       for tool1 in $tools1; do
        for tool2 in $tools2; do
        if [ $tool1 == $tool2 ]; then continue; fi
            seq 1 6 | while read i; do
                sbatch -A ${SLURM_ACCOUNT} -J $id.size$i.tool1.$i \
                    -o logs/$id.size$i.$tool1.$i.out \
                    --cpus-per-task=1 \
                    --mem-per-cpu=8g \
                    --wrap="
if [ -f 1.truvari/  ]; then 
    rm -rf 1.truvari/$id.$tool2.$tool1.size_$i ; 
fi                    
truvari bench -b 3.truvari_size/$id.$tool1.size_$i.vcf.gz -c 3.truvari_size/$id.$tool2.size_$i.vcf.gz \
    --dup-to-ins -o 1.truvari/$id.$tool2.$tool1.size_$i -r 2000 -C 2000 \
    --no-ref a --passonly --sizemax 1000000 --sizemin 50

cat 1.truvari/$id.$tool2.$tool1.size_$i/summary.json  |  
    ${CONDA_BASE}/bin/jq -r '[.precision, .recall, .f1, .gt_concordance] | @tsv' |
    sed \"s|^|$id\t$tool2.$tool1\t$i\t|\" >> $summtsv

rm 1.truvari/$id.$tool2.$tool1.size_$i
            "
        done
    done
    done
    done

### different type

mkdir 4.truvari_type
summtsv=truvari.pangenie-type.tsv
echo "ID tool cla precision recall f1 gt_concordance" | sed 's/ /\t/g' > $summtsv
cat ${PROJECT_ROOT}/stat_pan/hol.sample | 
    while read id; do
        for tool1 in $tools1; do
        for tool2 in $tools2; do
        if [ $tool1 == $tool2 ]; then continue; fi
        idd=$id
        for type in MNP DEL INS; do
                sbatch -A ${SLURM_ACCOUNT} -J $idd.$tool1.type \
                    -o logs/$idd.$tool1.type.out \
                    --cpus-per-task=1 \
                    --mem-per-cpu=8g \
                    --wrap="
if [ -f 1.truvari/  ]; then 
    rm -rf 1.truvari/$id.$tool2.$tool1.$type ; 
fi 
truvari bench -b 4.truvari_type/$id.$tool1.$type.vcf.gz -c 4.truvari_type/$id.$tool2.$type.vcf.gz \
    --dup-to-ins -o  1.truvari/$id.$tool2.$tool1.$type -r 2000 -C 2000 \
    --no-ref a --passonly --sizemax 1000000 --sizemin 50

cat 1.truvari/$id.$tool2.$tool1.$type/summary.json  |  
    ${CONDA_BASE}/bin/jq -r '[.precision, .recall, .f1, .gt_concordance] | @tsv' |
    sed \"s|^|$id\t$tool2.$tool1\t$type\t|\" >> $summtsv
done
            "
            sleep 0.1
        done
    done
    done
    done
