# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
: "${CLUSTER_HOST:?CLUSTER_HOST is unset - see config.sh.example at the repository root}"
# --------------------------

cd ${PROJECT_ROOT}/stat_pan/2.vcf_stats/0.filter

# awk 'FNR==NR{split($0,b," ");printf "ID\t";for(j in b){printf "%s\t", b[j]} ;printf "\n";next} 
#     {printf $8"\t";for(i=10;i<=NF;i++){a[$i]++}; for(j in b){printf "%s\t", a[b[j]]}; delete a;printf "\n"}'  \
#     <(echo '. 0 1 .|. 0|. .|0 1|. .|1 0|0 0|1 1|0 1|1') \
#     <(zcat bovHol-2024-08-13.anno_biallelic.filtered.vcf.gz | grep -v "#") > bovHol-2024-08-13.anno_biallelic.filtered.count


# count by position
# zcat bovHol-2024-08-13.anno_biallelic.filtered.vcf.gz | grep -v "#" | 
#     awk 'BEGIN {print "ID\tCHR\tSTART\tEND\tLENGTH\tWALK\tTYPE\tPUB0\tPUB1\tHOL_0\tHOL_1\tCLASS"} {
#     split($8,bc,/-/);
#     printf $8"\t"$1"\t"$2"\t"$2+bc[5]-1"\t"bc[5]"\t"$3"\t"bc[3]"\t";
#     for(i=10;i<=NF;i++){a[$i]++};
#     hol_1=a["1|."] + a[".|1"] + a["1|0"] + a["0|1"] + a["1|1"];
#     hol_0=a["0|."] + a[".|0"] + a["1|0"] + a["0|1"] + a["0|0"];
#     printf a["0"]"\t"a["1"]"\t"hol_0"\t"hol_1"\t";
#     if(a["1"] > 0 && hol_1 == 0){printf "BovMul"} else if(a["1"] == 0 && hol_1 > 0){printf "Holstein"} else {printf "Shared"};
#     delete a;
#     printf "\n";
#     }' > bovHol-2024-08-13.anno_biallelic.filtered.class

# count by position + GCA_021347905

cat > class.id.awk <<'EOF'
BEGIN{
FS=OFS="\t";
print "ID\tCHR\tSTART\tEND\tLENGTH\tWALK\tTYPE\tPUB0\tPUB1\tHOL_0\tHOL_1\tCLASS"
}
$1 ~ /##/{
    print
    next
}
$1 ~ /#/{
    print "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"
    print "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"
    print "##INFO=<ID=.,Number=1,Type=String,Description=\"aaa\">"
    print;
    next
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
    #print $1,$2,$2+leng,leng,typ,"-",typ"-"NR 
    if(leng < 50){
        typ="s"typ
    }
    $8="SVTYPE="typ";END="$2+leng";SVLEN="leng";"$8


    printf $8"\t"$1"\t"$2"\t"$2+leng-1"\t"leng"\t"$3"\t"typ"\t";
    for(i=10;i<=NF;i++){if(i != 12){a[$i]++}};
    hol_1=a["1|."] + a[".|1"] + a["1|0"] + a["0|1"] + a["1|1"];
    hol_0=a["0|."] + a[".|0"] + a["1|0"] + a["0|1"] + a["0|0"];
    if($12 == "1"){hol_1=hol_1+1} else if($12 == "0"){hol_0=hol_0+1};
    printf a["0"]"\t"a["1"]"\t"hol_0"\t"hol_1"\t";
    if(a["1"] > 0 && hol_1 == 0){printf "BovMul"} else if(a["1"] == 0 && hol_1 > 0){printf "Holstein"} else {printf "Shared"};
    delete a;
    printf "\n";
}
EOF

zcat bovHol-2024-08-13.anno_biallelic.filtered.vcf.gz | grep -v "#"  | awk -f class.id.awk > bovHol-2024-08-13.anno_biallelic.filtered.class


cat > class.stat.awk <<'EOF'
BEGIN{
    print "pan type genome_count genome_length autosome_count autosome_length"
} 

{
    gc[$12,$7]++;gl[$12,$7]+=$5;if($2 <= 29){ac[$12,$7]++;al[$12,$7]+=$5}

}

END {
    for (comb in gc) {
        split(comb,sep,SUBSEP);
        print sep[1], sep[2], gc[sep[1],sep[2]], gl[sep[1],sep[2]], ac[sep[1],sep[2]], al[sep[1],sep[2]]
        gca[sep[1]]+=gc[sep[1],sep[2]];
        gla[sep[1]]+=gl[sep[1],sep[2]];
        aca[sep[1]]+=ac[sep[1],sep[2]];
        ala[sep[1]]+=al[sep[1],sep[2]];
    }
    for (j in gca) {
        print j, "ALL", gca[j], gla[j], aca[j], ala[j]
    }
}
EOF

cat bovHol-2024-08-13.anno_biallelic.filtered.class | awk -F "\t" 'NR > 1'| awk -F "\t" -f class.stat.awk > bovHol-2024-08-13.anno_biallelic.filtered.class.stat

cat bovHol-2024-08-13.anno_biallelic.filtered.class | awk -F "\t" 'NR > 1 && ($9 > 1 || $11 > 1)' | awk -F "\t" -f class.stat.awk > bovHol-2024-08-13.anno_biallelic.filtered.af1class.stat

cat bovHol-2024-08-13.anno_biallelic.filtered.class | awk -F "\t" 'NR > 1 && ($9 == 12 || $11 == 21)' | awk -F "\t" -f class.stat.awk > bovHol-2024-08-13.anno_biallelic.filtered.afaclass.stat

cat bovHol-2024-08-13.anno_biallelic.filtered.class | awk  -F "\t"  '$5 > 50 {print $2,$3,$4,$7,$9,$11,$12}' OFS="\t" > bovHol-2024-08-13.anno_biallelic.filtered.class.sv

awk -F "\t" '$5 > 50 {print $12,$5,$7,$9+$11}' bovHol-2024-08-13.anno_biallelic.filtered.class > bovHol-2024-08-13.anno_biallelic.filtered.class.hist.in


awk 'NR > 1 && $7 != "BovMul"' bovHol-2024-08-13.anno_biallelic.filtered.class.sv | sortBed | mergeBed -i - -c 4,7 -o distinct,distinct > bovHol-2024-08-13.anno_biallelic.filtered.class.holsvr
awk 'NR > 1 && $7 != "Holstein"' bovHol-2024-08-13.anno_biallelic.filtered.class.sv | sortBed | mergeBed -i - -c 4,7 -o distinct,distinct > bovHol-2024-08-13.anno_biallelic.filtered.class.pubsvr

awk 'NR > 1' bovHol-2024-08-13.anno_biallelic.filtered.class.sv | sortBed | mergeBed -i - -c 4,7 -o distinct,distinct > bovHol-2024-08-13.anno_biallelic.filtered.class.allsvr

scp ${CLUSTER_HOST}:${PROJECT_ROOT}/stat_pan/2.vcf_stats/0.filter/bovHol-2024-08-13.anno_biallelic.filtered.af1class.stat .

### 
cat > frq.stat.awk <<'EOF'
BEGIN{
    print "pan type frq genome_count genome_length autosome_count autosome_length"
} 

{
    gc[$12,$7,$9+$11]++;
    gl[$12,$7,$9+$11]+=$5;
    if($2 <= 29){
        ac[$12,$7,$9+$11]++;
        al[$12,$7,$9+$11]+=$5;
        }
}

END {
    for (comb in gc) {
        split(comb,sep,SUBSEP);
        print sep[1], sep[2], sep[3], gc[sep[1],sep[2],sep[3]], gl[sep[1],sep[2],sep[3]], ac[sep[1],sep[2],sep[3]], al[sep[1],sep[2],sep[3]]
    }
}
EOF

cat bovHol-2024-08-13.anno_biallelic.filtered.class | awk -F "\t" 'NR > 1'| awk -F "\t" -f frq.stat.awk > bovHol-2024-08-13.anno_biallelic.filtered.class_frq.stat

cat bovHol-2024-08-13.anno_biallelic.filtered.class | awk -F "\t" '$5 > 50' | awk -F "\t" -f frq.stat.awk > bovHol-2024-08-13.anno_biallelic.filtered.class_frq.sv.stat

