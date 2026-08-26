#!/bin/awk -f
# summarize_svtype.awk
# Count number and cumulative length of SVs by type, total and autosomes only

BEGIN {
    FS = OFS = "\t"
}

{
    type = $5
    len = $6
    chr = $2

    gc[type]++        # genome-wide count
    gl[type] += len   # genome-wide length

    if (chr <= 29) {
        ac[type]++    # autosome count
        al[type] += len
    }
}

END {
    for (t in gc) {
        print idd, t, gc[t], gl[t], ac[t], al[t]
    }
}
