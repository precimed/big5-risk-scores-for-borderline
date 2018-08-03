# Merge 23andMe info and sumstat files

# First input file should be all_snp_info-4.0.txt, which contains the following fields
# all.data.id	gt.data.id	im.data.id	assay.name	scaffold	position	alleles	ploidy	strand	cytoband	gene.context	is.v1	is.v2	is.v3	is.v4	h550	omni

# Second file should be *-4.0.dat file, which contains
# all.data.id	src	pvalue	effect	stderr	pass	dose.b.0	dose.b.1	AA.0	AB.0	BB.0	AA.1	AB.1	BB.1

BEGIN {
    OFS = "\t"
    npass = 0

    if (length(outfile) == 0) {
        print "No output file specified"
        exit
    }

    if (length(remove_MHC) == 0 || remove_MHC > 1 || remove_MHC < 0) {
        print "Invalid option for including (1) or excluding (0) MHC region: " remove_MHC
        exit
    }
    if (remove_MHC == 1) print "Excluding MHC region"
    else print "Including MHC region"

    print "SNP", "CHR", "BP", "A1",  "A2", "PVAL", "EFFECT", "STDERR" > outfile
}


{
    # Check that we are operating on the correct files
    if (FNR == NR) {
        if (! FILENAME ~ /all_snp_info/) {
            print "First file is " FILENAME " - check inputs"
            exit
        }
    }
    else {
        if (! FILENAME ~ /.dat$/) {
            print "Second file is " FILENAME " - check inputs"
            exit
        }
    }
    

    # Skip header
    if (FNR == 1) {
        print "Processing " FILENAME "..."
        next
    }
    
    # Process SNP info file
    if (FNR == NR) {
    
        # 'all.data.id' is the array key for all fields.
        id = $1
        snpid[id] = $4
        chrname = $5; sub(/^chr/, "", chrname); chr[id] = chrname
        bp[id] = $6

        alleles = $7
        split(alleles, a1a2, "/")
        a1[id] = a1a2[1]
        a2[id] = a1a2[2]

    }

    # Process GWAS result file
    else {

        id = $1
        pval = $3
        effect = $4
        err = $5
        pass = $6
        
        # Check 23andMe "pass" and rs number
        if (pass != "Y") {nfail["no-pass"]++; next}
        if (! (snpid[id] ~ /^rs/)) {nfail["no-rs"]++; next}
        
        # No p-value or no effect size
        if (! (pval ~ /^[+,-]?[0-9]+/)) {print; nfail["no-pval"]++; next}
        if (! (effect ~ /^[+,-]?[0-9]+/)) {nfail["no-effect"]++; next}

        # Remove ambiguous SNPs
        if ((a1[id] == "A" && a2[id] == "T") || (a1[id] == "T" && a2[id] == "A")) {nfail["ambig"]++; next}
        if ((a1[id] == "C" && a2[id] == "G") || (a1[id] == "G" && a2[id] == "C")) {nfail["ambig"]++; next}

        # Regression fail
        if (effect > 10 || effect < -10) {nfail["regress"]++; next}

        if (remove_MHC) {
            if (chr[id] == 6 && bp[id] > 26000000 && bp[id] < 34000000) {nfail["in-mhc"]++; next}
        }

        # Write output
        print snpid[id], chr[id], bp[id], a1[id], a2[id], pval, effect, err > outfile
        npass++

    }
}


END {
    close(outfile)
    print "Output written to \"" outfile "\""
    print "SNPs passed:", npass #, "SNPs failed:", nfail
    print "SNPs failed:"
    for (f in nfail) {
        print "  " f ":", nfail[f]
    }

    # Exit with nonzero code in case of obvious trouble
    if (npass < 1) {
        print "Error: No SNPs passed selection"
        exit 1
    }
}
