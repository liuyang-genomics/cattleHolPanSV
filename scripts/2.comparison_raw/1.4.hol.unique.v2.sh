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


"C:\Users\${USER}\OneDrive - University of Maryland\Data\2024-02-07.cattleLR-SR-GWAS\2.analyses\1.2x.anno.stats"
scp ${CLUSTER_HOST}:${PROJECT_ROOT}/stat_pan/2.vcf_stats/0.filter/bovHol-2024-08-13.anno_biallelic.filtered.class_frq.stat .
scp ${CLUSTER_HOST}:${PROJECT_ROOT}/stat_pan/2.vcf_stats/0.filter/bovHol-2024-08-13.anno_biallelic.filtered.class_frq.sv.stat .


vg chunk -r 1026333:1026338 -x hol-pg2hic-2024-05-22.gbz -O gfa -t 8 > hol-pg2hic-2024-05-22.node.1026333-1026338.gfa

install graphviz

vg chunk -r 1026333:1026338 -x bovHol-2024-08-13.gbz -O gfa -t 8 > bovHol-2024-08-13.node.1026333-1026338.gfa

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa \
    | vg view -d - \
    | dot -Tsvg -o x.svg
chromium-browser x.svg

bovHol
4681 4716
162121 162166
162000 162200

hol
1027060 1027074

vg chunk -r 4681:4716 -x ${PROJECT_ROOT}/minigraph-cactus/bovHol-2024-08-13/bovHol-2024-08-13.gbz -O gfa -t 8 > bovHol-2024-08-13.node.4681-4716.gfa
vg chunk -r 1027060:1027074 -x ${PROJECT_ROOT}/minigraph-cactus/hol-pg2hic-2024-05-22/hol-pg2hic-2024-05-22.gbz -O gfa -t 8 > hol-pg2hic-2024-05-22.node.1027060-1027074.gfa


head bovHol-2024-08-13.anno_biallelic.filtered.class.allsvr
1       160002  160071  COMPLEX BovMul,Shared
1       160769  161363  COMPLEX,INS     BovMul,Holstein,Shared
1       162077  162337  INS     BovMul,Holstein,Shared ########
1       163724  164075  DEL,INS Holstein,Shared
1       164396  164464  COMPLEX BovMul
1       164486  164871  DEL     BovMul
1       164985  165236  DEL     BovMul
1       166650  166742  COMPLEX Shared
1       167328  167419  DEL     Holstein
1       171515  171673  DEL     BovMul


cat bovHol-2024-08-13.anno_biallelic.filtered.class | awk '$2 == 1 && $3 >= 162077 && $4 <= 162337'
AT=>4674>4679,>4674>4678>4679;ID=1-162077-INS->4674>4678>4679-1 1       162077  162077  1       >4674>4679      INS     7       1       19      0       BovMul
AT=>4674>4679,>4674>4675>4677>4678>4679;ID=1-162077-INS->4674>4675>4677>4678>4679-115   1       162077  162191  115     >4674>4679      INS     2       6       8       16      Shared
AT=>4674>4679,>4674>4675>4676>4677>4678>4679;ID=1-162077-INS->4674>4675>4676>4677>4678>4679-138 1       162077  162214  138     >4674>4679      INS     8               16      7       Holstein
AC=1;AF=0.0227273;AN=44;AT=>4679>4680>4681,>4679>4681;NS=27;LV=0;ID=1-162101-DEL->4679>4681-1   1       162101  162101  1       >4679>4681      DEL     7       1       19      0       BovMul
AC=1;AF=0.0227273;AN=44;AT=>4681>4682>4683,>4681>4683;NS=27;LV=0;ID=1-162105-DEL->4681>4683-1   1       162105  162105  1       >4681>4683      DEL     7       1       19      0       BovMul
AC=1;AF=0.0227273;AN=44;AT=>4683>4685,>4683<4684>4685;NS=27;LV=0;ID=1-162121-INS->4683<4684>4685-1      1       162121  162121  1       >4683>4685      INS     7       1       19      0       BovMul
AT=>4685>4711>4713,>4685>4693<4695>4696<4698>4699>4701<4703>4704>4706>4708>4710<4712>4713;ID=1-162157-INS->4685>4693<4695>4696<4698>4699>4701<4703>4704>4706>4708>4710<4712>4713-58     1       162157  162214      58      >4685>4713      INS     7       1       19      0       BovMul
AT=>4685>4711,>4685>4693>4694>4696>4697>4699>4700>4701>4702>4704>4705>4706>4707>4708>4709>4710>4711;ID=1-162154-INS->4685>4693>4694>4696>4697>4699>4700>4701>4702>4704>4705>4706>4707>4708>4709>4710>4711-68
        1       162154  162221  68      >4685>4713      INS     2       6       8       16      Shared
AT=>4685>4711,>4685>4686>4688>4689>4691>4692>4693>4694>4696>4697>4699>4700>4701>4702>4704>4705>4706>4707>4708>4709>4710>4711;ID=1-162154-INS->4685>4686>4688>4689>4691>4692>4693>4694>4696>4697>4699>4700>4701>4702>4704>4705>4706>4707>4708>4709>4710>4711-136     1       162154  162289  136     >4685>4713      INS     8               16      7       Holstein
AT=>4685>4711,>4685>4686>4687>4689>4690>4692>4693>4694>4696>4697>4699>4700>4701>4702>4704>4705>4706>4707>4708>4709>4710>4711;ID=1-162154-INS->4685>4686>4687>4689>4690>4692>4693>4694>4696>4697>4699>4700>4701>4702>4704>4705>4706>4707>4708>4709>4710>4711-184     1       162154  162337  184     >4685>4713      INS     7       1       19      0       BovMul
AC=1;AF=0.0227273;AN=44;AT=>4714>4715>4716,>4714>4716;NS=27;LV=0;ID=1-162166-DEL->4714>4716-3   1       162166  162168  3       >4714>4716      DEL     7       1       19      0       BovMul
AC=1;AF=0.0227273;AN=44;AT=>4716>4717>4718,>4716>4718;NS=27;LV=0;ID=1-162175-DEL->4716>4718-1   1       162175  162175  1       >4716>4718      DEL     7       1       19      0       BovMul
AC=1;AF=0.0227273;AN=44;AT=>4718>4719>4721,>4718<4720>4721;NS=27;LV=0;ID=1-162178-SNV->4718<4720>4721-1 1       162178  162178  1       >4718>4721      SNV     7       1       19      0       BovMul
AT=>4721>4722>4723>4726,>4721<4724>4725>4726;ID=1-162191-COMPLEX->4721<4724>4725>4726-5 1       162191  162195  5       >4721>4726      COMPLEX 7       1       19      0       BovMul
AT=>4722>4723>4726,>4722>4725>4726;ID=1-162195-SNV->4722>4725>4726-1    1       162195  162195  1       >4721>4726      SNV     1       7       1       19      Shared
AC=1;AF=0.0227273;AN=44;AT=>4726>4727>4729,>4726<4728>4729;NS=27;LV=0;ID=1-162210-SNV->4726<4728>4729-1 1       162210  162210  1       >4726>4729      SNV     7       1       19      0       BovMul
AC=1;AF=0.0227273;AN=44;AT=>4729>4730>4731,>4729>4731;NS=27;LV=0;ID=1-162223-DEL->4729>4731-1   1       162223  162223  1       >4729>4731      DEL     7       1       19      0       BovMul
AC=42;AF=0.954545;AN=44;AT=>4731>4733>4734,>4731>4732>4734;NS=27;LV=0;ID=1-162231-SNV->4731>4732>4734-1 1       162231  162231  1       >4731>4734      SNV     1       7       1       19      Shared
AC=1;AF=0.0227273;AN=44;AT=>4734>4735>4736,>4734>4736;NS=27;LV=0;ID=1-162236-DEL->4734>4736-1   1       162236  162236  1       >4734>4736      DEL     7       1       19      0       BovMul
AC=1;AF=0.0227273;AN=44;AT=>4736>4738,>4736<4737>4738;NS=27;LV=0;ID=1-162249-INS->4736<4737>4738-1      1       162249  162249  1       >4736>4738      INS     7       1       19      0       BovMul
AC=1;AF=0.0227273;AN=44;AT=>4738>4739>4741,>4738<4740>4741;NS=27;LV=0;ID=1-162265-COMPLEX->4738<4740>4741-5     1       162265  162269  5       >4738>4741      COMPLEX 7       1       19      0       BovMul
AC=1;AF=0.0227273;AN=44;AT=>4741>4742>4744,>4741<4743>4744;NS=27;LV=0;ID=1-162273-SNV->4741<4743>4744-1 1       162273  162273  1       >4741>4744      SNV     7       1       19      0       BovMul


vg chunk -r 4674:4741 -x ${PROJECT_ROOT}/minigraph-cactus/bovHol-2024-08-13/bovHol-2024-08-13.gbz -O gfa -t 8 > bovHol-2024-08-13.node.4674-4741.gfa
4674-4741

vg chunk -r 4007561:4007769 -x ${PROJECT_ROOT}/minigraph-cactus/bovinePan-2024-08-13/bovinePan-2024-08-13.gbz -O gfa -t 8 > bovHol-2024-08-13.node.4007561-4007769.gfa
4007561:4007769

#### 

cat bovHol-2024-08-13.anno_biallelic.filtered.class.allsvr | awk '$3-$2 + 1 > 100kb && $5 ~ ","' | head
1       160002  160071  COMPLEX BovMul,Shared
1       160769  161363  COMPLEX,INS     BovMul,Holstein,Shared
1       162077  162337  INS     BovMul,Holstein,Shared # 
1       163724  164075  DEL,INS Holstein,Shared
1       173457  178857  COMPLEX,DEL,INS BovMul,Holstein,Shared
1       246328  248968  INS     BovMul,Shared
1       383105  438752  COMPLEX,DEL,INS BovMul,Holstein,Shared
1       450233  457365  INS     Holstein,Shared
1       474270  474665  INS     Holstein,Shared
1       480376  481323  INS     Holstein,Shared

cat bovHol-2024-08-13.anno_biallelic.filtered.class | awk '$2 == 1 && $3 >= 162077 && $4 <= 162337'

cat bovinePan-2024-08-13.anno_biallelic.filtered.svtype | awk '$2 == 1 && $3 >= 162077 && $3 + $6 -1 <= 162337'

cat hol-pg2hic-2024-05-22.anno_biallelic.filtered.svtype | awk '$2 == 1 && $3 >= 162077 && $3 + $6 -1 <= 162337'
>1027060>1027063 1 162093 >1027060>1027061>1027063 SNV 1
>1027063>1027071 1 162105 >1027063<1027069>1027071 COMPLEX 186
>1027063>1027071 1 162107 >1027064>1027066>1027067 SNV 1
>1027071>1027074 1 162195 >1027071<1027072>1027074 SNV 1
>1027074>1027077 1 162231 >1027074<1027075>1027077 SNV 1






cat ${PROJECT_ROOT}/minigraph-cactus/bovinePan-2024-08-13/bovinePan-2024-08-13.anno_biallelic.svtype | awk '$2 == 1 && $3 >= 162077 && $3 + $6 -1 <= 162337'
>4007561>4007741 1 162077 >4007709<4007710>4007711 INS 1
>4007561>4007741 1 162101 >4007711>4007713 DEL 1
>4007561>4007741 1 162105 >4007713>4007715 DEL 1
>4007561>4007741 1 162121 >4007715<4007716>4007717 INS 1
>4007561>4007741 1 162157 >4007717>4007718>4007719>4007721<4007723>4007724<4007726>4007727>4007729<4007731>4007732>4007734>4007736>4007738<4007740>4007741 INS 58
>4007742>4007744 1 162166 >4007742>4007744 DEL 3
>4007744>4007746 1 162175 >4007744>4007746 DEL 1
>4007746>4007749 1 162178 >4007746<4007748>4007749 SNV 1
>4007749>4007754 1 162191 >4007749<4007751>4007752>4007754 COMPLEX 5
>4007749>4007754 1 162195 >4007750>4007752>4007754 SNV 1
>4007754>4007757 1 162210 >4007754<4007756>4007757 SNV 1
>4007757>4007759 1 162223 >4007757>4007759 DEL 1
>4007759>4007762 1 162231 >4007759>4007760>4007762 SNV 1
>4007762>4007764 1 162236 >4007762>4007764 DEL 1
>4007764>4007766 1 162249 >4007764<4007765>4007766 INS 1
>4007766>4007769 1 162265 >4007766<4007768>4007769 COMPLEX 5
>4007769>4007772 1 162273 >4007769<4007771>4007772 SNV 1

cat ${PROJECT_ROOT}/minigraph-cactus/hol-pg2hic-2024-05-22/hol-pg2hic-2024-05-22.anno_biallelic.svtype | awk '$2 == 1 && $3 >= 162077 && $3 + $6 -1 <= 162337'
>1027060>1027063 1 162093 >1027060>1027061>1027063 SNV 1
>1027063>1027071 1 162105 >1027063<1027069>1027071 COMPLEX 186
>1027063>1027071 1 162107 >1027064>1027066>1027067 SNV 1
>1027071>1027074 1 162195 >1027071<1027072>1027074 SNV 1
>1027074>1027077 1 162231 >1027074<1027075>1027077 SNV 1


cat ${PROJECT_ROOT}/minigraph-cactus/bovHol-2024-08-13/bovHol-2024-08-13.anno_biallelic.svtype | awk '$2 == 1 && $3 >= 162077 && $3 + $6 -1 <= 162337'
>4674>4679 1 162077 >4674>4678>4679 INS 1
>4674>4679 1 162077 >4674>4675>4677>4678>4679 INS 115
>4674>4679 1 162077 >4674>4675>4676>4677>4678>4679 INS 138
>4679>4681 1 162101 >4679>4681 DEL 1
>4681>4683 1 162105 >4681>4683 DEL 1
>4683>4685 1 162121 >4683<4684>4685 INS 1
>4685>4713 1 162157 >4685>4693<4695>4696<4698>4699>4701<4703>4704>4706>4708>4710<4712>4713 INS 58
>4685>4713 1 162154 >4685>4693>4694>4696>4697>4699>4700>4701>4702>4704>4705>4706>4707>4708>4709>4710>4711 INS 68
>4685>4713 1 162154 >4685>4686>4688>4689>4691>4692>4693>4694>4696>4697>4699>4700>4701>4702>4704>4705>4706>4707>4708>4709>4710>4711 INS 136
>4685>4713 1 162154 >4685>4686>4687>4689>4690>4692>4693>4694>4696>4697>4699>4700>4701>4702>4704>4705>4706>4707>4708>4709>4710>4711 INS 184
>4714>4716 1 162166 >4714>4716 DEL 3
>4716>4718 1 162175 >4716>4718 DEL 1
>4718>4721 1 162178 >4718<4720>4721 SNV 1
>4721>4726 1 162191 >4721<4724>4725>4726 COMPLEX 5
>4721>4726 1 162195 >4722>4725>4726 SNV 1
>4726>4729 1 162210 >4726<4728>4729 SNV 1
>4729>4731 1 162223 >4729>4731 DEL 1
>4731>4734 1 162231 >4731>4732>4734 SNV 1
>4734>4736 1 162236 >4734>4736 DEL 1
>4736>4738 1 162249 >4736<4737>4738 INS 1
>4738>4741 1 162265 >4738<4740>4741 COMPLEX 5
>4741>4744 1 162273 >4741<4743>4744 SNV 1




### test length
