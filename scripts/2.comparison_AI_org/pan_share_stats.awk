#!/bin/awk -f
# Count and compare intersected SVs from pan vs tool BED
# Stratify by SV type, length bin, and whether region is exactly matched

BEGIN {
    lenCate = "200,1000,10000,100000,1000000"
}

function count_num(bed, chr, type, lenCate, extp) {
    array["count", bed, "ALL", "ALL", lenCate, extp]++
    array["count", bed, chr, "ALL", lenCate, extp]++
    array["count", bed, "ALL", type, lenCate, extp]++
    array["count", bed, chr, type, lenCate, extp]++
}

function count_len(bed, chr, type, lenCate, extp, len) {
    array["length", bed, "ALL", "ALL", lenCate, extp] += len
    array["length", bed, chr, "ALL", lenCate, extp] += len
    array["length", bed, "ALL", type, lenCate, extp] += len
    array["length", bed, chr, type, lenCate, extp] += len
}

function count_point(bed, chr, type, lenCate, extp, len) {
    count_num(bed, chr, type, lenCate, "no")
    count_len(bed, chr, type, lenCate, "no", len)
    if (extp) {
        count_num(bed, chr, type, lenCate, "yes")
        count_len(bed, chr, type, lenCate, "yes", len)
    }
}

function each_count_point(bed, chr, type, lenCate, extp, len) {
    count_point(bed, chr, type, "ALL", extp, len)
    split(lenCate, bins, ",")
    for (i = 1; i <= length(bins); i++) {
        if (len <= bins[i]) {
            count_point(bed, chr, type, bins[i], extp, len)
        }
    }
}

{
    if ($5 != $12) $5 = $5 "-" $12
    if ($15 < 0) $15 = 1
    uniq_pan[$7]++
    uniq_tool[$14]++
    extp = ($2 == $9 && $3 == $10)
    chr = $1
    type = $5
    len = $15

    if (uniq_pan[$7] == 1) {
        bed = "pan"
        each_count_point(bed, chr, type, lenCate, extp, len)
    }
    if (uniq_tool[$14] == 1) {
        bed = "tool"
        each_count_point(bed, chr, type, lenCate, extp, len)
    }
}

END {
    for (key in array) {
        split(key, fields, SUBSEP)
        for (i = 1; i <= length(fields); i++) {
            printf "%s	", fields[i]
        }
        print array[key]
    }
}
