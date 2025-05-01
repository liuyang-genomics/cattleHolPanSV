#!/bin/awk -f
# Count the number and cumulative length of SVs per chromosome, type, and length category

BEGIN {
    lenCate = "200,1000,10000,100000,1000000"
}

function count_num(bed, chr, type, lenCate, extp) {
    array["count", bed, "ALL", "ALL", lenCate, extp]++
    array["count", bed, chr, "ALL", lenCate, extp]++
    array["count", bed, "ALL", type, lenCate, extp]++
    array["count", bed, chr, type, lenCate, extp]++
}

function count_len(bed, chr, type, lenCate, extp, length) {
    array["length", bed, "ALL", "ALL", lenCate, extp] += length
    array["length", bed, chr, "ALL", lenCate, extp] += length
    array["length", bed, "ALL", type, lenCate, extp] += length
    array["length", bed, chr, type, lenCate, extp] += length
}

function count_point(bed, chr, type, lenCate, extp, length) {
    count_num(bed, chr, type, lenCate, "no")
    count_len(bed, chr, type, lenCate, "no", length)
}

function each_count_point(bed, chr, type, lenCate, extp, length) {
    count_point(bed, chr, type, "ALL", extp, length)
    split(lenCate, bins, ",")
    for (i = 1; i <= length(bins); i++) {
        if (length <= bins[i]) {
            count_point(bed, chr, type, bins[i], extp, length)
        }
    }
}

{
    bed = ""
    extp = ""
    chr = $1
    type = $5
    shareLen = $4
    each_count_point(bed, chr, type, lenCate, extp, shareLen)
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
