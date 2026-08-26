#!/bin/awk -f
# Annotate VCF with SVTYPE, END, SVLEN, and remove SVs <= 50bp

BEGIN {
    FS = OFS = "\t"
}

{
    if ($1 ~ /^##/) {
        print
    } else if ($1 ~ /^#/) {
        print "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
        print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"
        print "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"
        print
    } else {
        ref_len = length($4)
        alt_len = length($5)
        len = 0

        if (alt_len > ref_len && ref_len == 1) {
            type = "INS"
            len = alt_len - 1
        } else if (alt_len < ref_len && alt_len == 1) {
            type = "DEL"
            len = ref_len - 1
        } else if (alt_len == ref_len && ref_len == 1) {
            type = "SNP"
            len = 1
        } else {
            type = "MNP"
            len = alt_len - 1
        }

        if (len > 50) {
            $8 = "SVTYPE=" type ";END=" ($2 + len) ";SVLEN=" len ";" $8
            print
        }
    }
}
