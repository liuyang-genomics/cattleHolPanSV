




${PANEL_DIR}/her_20201802/2.assembly/her_20201802.quast-hifiasm-bp/contigs_reports

ls ${PANEL_DIR}/*/2.assembly/*.busco-hifiasm-bp/short_summary.specific.mammalia_odb10.*.busco-hifiasm-bp.json | 
    while read id ; do 
        idd=$(echo $id | cut -d/ -f6 )
        cat $id | grep "C:" | sed "s/$/\t${idd}/g"
    done


ls ${PANEL_DIR}/*/2.assembly/*.busco-hifiasm-*/short_summary.specific.mammalia_odb10.*.busco-hifiasm-*.json | 
    while read id ; do 
        idd=$(echo $id | cut -d/ -f8 )
        cat $id | grep "Contigs N50" | grep -v "of" | sed "s/$/\t${idd}/g"
    done


ls ${PANEL_DIR}/*/2.assembly/*.busco-hifiasm-*/short_summary.specific.mammalia_odb10.*.busco-hifiasm-*.json | 
    while read id ; do 
        idd=$(echo $id | cut -d/ -f8 )
        cat $id | grep "Total length" | grep -v "of" | sed "s/$/\t${idd}/g"
    done


ls *hifi_hi2c.hic*.fa.gfastats | 
    while read id ; do 
        idd=$(echo $id | sed 's/.p_ctg.fa.gfastats//' )
        cat $id | grep "Contig N50" | grep -v "of" | sed "s/$/\t${idd}/g"
    done

ls *hifi_hi2c.hic*.fa.gfastats | 
    while read id ; do 
        idd=$(echo $id | sed 's/.p_ctg.fa.gfastats//' )
        cat $id | grep "Total scaffold length" | grep -v "of" | sed "s/$/\t${idd}/g"
    done


ls ${PROJECT_ROOT}/hifi-hic-assembly/sample_*.hifi_hic.hic*/run_mammalia_odb10/short_summary.json | 
    while read id ; do 
        idd=$(echo $id | cut -d/ -f7 )
        cat $id | grep "Contigs N50" | grep -v "of" | sed "s/$/\t${idd}/g"
    done

ls ${PROJECT_ROOT}/hifi-hic-assembly/sample_*.hifi_hic.hic*/run_mammalia_odb10/short_summary.json | 
    while read id ; do 
        idd=$(echo $id | cut -d/ -f7 )
        cat $id | grep "Total length" | grep -v "of" | sed "s/$/\t${idd}/g"
    done





ls ${PANEL_DIR}/*/2.assembly/*.inspector-hifiasm-bp/summary_statistics | 
    while read id ; do 
        idd=$(echo $id | cut -d/ -f6 )
        cat $id | grep "Mapping rate in large contigs /%" | sed "s/$/\t${idd}/g"
    done

ls ${PANEL_DIR}/*/2.assembly/*.inspector-hifiasm-bp/summary_statistics | 
    while read id ; do 
        idd=$(echo $id | cut -d/ -f6 )
        cat $id | grep "Depth in large conigs" | sed "s/$/\t${idd}/g"
    done


ls ${PANEL_DIR}/*/2.assembly/*.inspector-hifiasm-bp/summary_statistics | 
    while read id ; do 
        idd=$(echo $id | cut -d/ -f6 )
        cat $id | grep "QV" | sed "s/$/\t${idd}/g"
    done

ls ${PANEL_DIR}/*/2.assembly/*.inspector-hifiasm-bp/summary_statistics | 
    while read id ; do 
        idd=$(echo $id | cut -d/ -f6 )
        cat $id | grep "N50" | grep -v "of" | sed "s/$/\t${idd}/g"
    done



${PANEL_DIR}/her_20201802/2.assembly/her_20201802.quast-hifiasm-bp/contigs_reports

# total length / total length without N's / maximal covered length: 

ls ${PANEL_DIR}/s*/2.assembly/s*.quast-hifiasm*/genome_stats/genome_info.txt |
    while read id; do
        idd=$(echo $id | cut -d/ -f8)
    cat $id  | grep "X (total length: " | cut -d " " -f 4,10,15 | sed "s|^|$idd\t|"
    done

