# Re-write a sumstats file to conform with the strand orientation in a
# genotype file

BEGIN {
    OFS = "\t"
    nflipped = 0
    nunknown = 0

    if (length(outfile) == 0) {
        print "No output file specified"
        exit
    }
}

{

    # Process map file, which defines the correct A1 and A2 alleles
    if (FNR == NR) {
        id = $1
        a1[id] = $3
        a2[id] = $4
    }

    # Process score file
    else {

        # Header
        if (FNR == 1) print > outfile

        id = $1
        a1_tmp = $4
        a2_tmp = $5

        # Match
        if (a1_tmp == a1[id] && a2_tmp == a2[id]) print > outfile

        # Flip
        else if (a1_tmp == a2[id] && a2_tmp == a1[id]) {
            $4 = a2_tmp
            $5 = a1_tmp
            print > outfile
            nflipped++
        }

        # Complete mismatch
        else {
            nunknown++
        }
    }
}

END {
    close(outfile)
    print "Output written to \"" outfile "\""
    print "Flipped SNPs:", nflipped, "Unknown SNPs (skipped):", nunknown
}
