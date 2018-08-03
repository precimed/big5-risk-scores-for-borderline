#! bin/bash

set -o errexit

basedir=$(pwd)

# The big-5 traits
traits=("agreeableness" "conscientiousness" "extraversion" "neuroticism" "openness")
#traits=("agreeableness")

# Exclude SNPs in MHC region 
remove_MHC_from_sumstat="1"

# p-value thresholds for scoring
pval_thresholds="${basedir}/pvalue-ranges.txt"

# Genotype file name pattern
genotypes="${basedir}/genotypes/bpd_bom3_eur_sr-qc.hg19.ch.fl.bg"

# 1000 Genome reference files
reference_1kg="/Users/steffma/work/1000Genome"

# For macOS where gawk is not standard
if [[ ! -x "$(command -v gawk)" ]]; then
    export PATH=$PATH:/Users/steffma/software/gawk-4.2.1/bin
fi



# Get .map files from genotype input (to check for strand consistency)
if [[ ! -f ${genotypes}.recode.chr-1.map ]]; then
    plink -bfile ${genotypes} --recode beagle --out ${genotypes}.recode --allow-no-sex
    rm ${genotypes}.recode.chr-*.dat
fi



# Loop over traits and do the calculation
for trait in ${traits[@]}; do

    mkdir -p ${trait}; cd ${trait}

    scorefile="personality_${trait}_sumstats_merged.tsv"
    

    # -------------------------------------------------------------------------
    # Preparation
    # -------------------------------------------------------------------------

    # Merge and format 23andMe data files
    if [[ ! -f $scorefile ]]; then
        echo "Preparing score file for ${trait}"
        gawk -f ${basedir}/prepare-23andMe-sumstats.awk \
            -v outfile=${scorefile}.tmp -v remove_MHC=${remove_MHC_from_sumstat} \
            ${basedir}/Lo_2016_data/all_snp_info-4.0.txt \
            ${basedir}/Lo_2016_data/personality_${trait}-4.0.dat
        
        sort -k 1,1 ${scorefile}.tmp > ${scorefile}
        rm ${scorefile}.tmp
    fi


    # Prepare reference
    if [[ ! -f ${trait}.reference.bim ]]; then
        
        echo "Preparing reference file for ${trait}"

        # List of reference files
        merge_list="1kg_reference_merge_list.txt"
        find ${reference_1kg}/*.bed | cut -f 1 -d "." | sort > ${merge_list}

        # List of SNPs to be extracted from reference
        snplist="${trait}_snplist.txt"
        tail -n +2 ${scorefile} | cut -f 1 > ${snplist}

        # Merge binary files
        echo ""
        plink \
            --merge-list ${merge_list} \
            --extract ${snplist} \
            --make-bed \
            --allow-no-sex \
            --out ${trait}.reference
        echo ""

    fi


    # -------------------------------------------------------------------------
    # Clumping
    # -------------------------------------------------------------------------

    final_scorefile="$(basename ${scorefile} .tsv)_final.tsv"
    if [[ ! -f ${trait}.clumped ]]; then
        echo "Performing clumping for ${trait}"

        echo ""
        plink \
            --bfile ${trait}.reference \
            --clump ${scorefile} \
            --clump-field PVAL \
            --clump-p1 1 \
            --clump-p2 1 \
            --clump-r2 0.25 \
            --clump-kb 500 \
            --allow-no-sex \
            --out ${trait}
        echo ""


        coresnps="${trait}.clumped.coresnps.txt"
        awk '{if (NR > 1 && NF) print $3}' ${trait}.clumped | sort > ${coresnps}

        # Intersect with score file
        tail -n +2 $scorefile | join - ${coresnps} > ${final_scorefile}.tmp
        (echo $(head -n 1 ${scorefile}) && cat ${final_scorefile}.tmp) > ${final_scorefile}_strand_uncorr.tmp
        rm ${final_scorefile}.tmp
    
    fi



    # Correct strand flips
    echo "Correcting strand flips"
    cat ${genotypes}.recode.chr-*.map | awk -f ${basedir}/enforce-strand-consistency.awk \
        -v outfile=${final_scorefile} \
        - ${final_scorefile}_strand_uncorr.tmp 

    rm ${final_scorefile}_strand_uncorr.tmp


    # -------------------------------------------------------------------------
    # Scoring
    # -------------------------------------------------------------------------
    echo "Computing scores for ${trait}"
    
    echo ""
    plink \
        --bfile ${genotypes} \
        --allow-no-sex \
        --score ${final_scorefile} 1 5 7 header \
        --q-score-range ${pval_thresholds} ${final_scorefile} 1 6 header \
        --out scores-${trait}
    # Note to PLINK options:
    # --score needs the column numbers to read from: 
    #   first the SNP id, then allele 1, then effect size
    # Since 23andMe results gives the effect of the *B* allele, set allele
    # column number to 5, which corresponds to A2
    echo ""
    
    echo "Completed ${trait}"
    cd ${basedir}

done

echo "Done."
