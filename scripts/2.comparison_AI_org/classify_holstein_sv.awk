BEGIN {
    print "ID	CHR	START	END	LENGTH	WALK	TYPE	PUB0	PUB1	HOL_0	HOL_1	CLASS"
}
!/^#/ {
    split($8, bc, "-")
    chr = $1
    start = $2
    len = bc[5]
    end = start + len - 1
    type = bc[3]
    walk = $3
    id = $8

    for (i = 10; i <= NF; i++) {
        a[$i]++
    }

    pub0 = a["0"]
    pub1 = a["1"]
    hol0 = a["0|."] + a[".|0"] + a["0|0"] + a["0|1"] + a["1|0"]
    hol1 = a["1|."] + a[".|1"] + a["1|1"] + a["1|0"] + a["0|1"]

    class = (pub1 > 0 && hol1 == 0) ? "PubOnly" : ((pub1 > 0 && hol1 > 0) ? "Shared" : "HolOnly")

    print id, chr, start, end, len, walk, type, pub0, pub1, hol0, hol1, class
    delete a
}
