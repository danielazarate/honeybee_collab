#!/bin/bash -l

# this script calculates pi for a population and list of bams

# to run: ./pi_within_ancestry.sh AR01 R01.bams

POP="${1}"
BAM_LIST="${2}"
DIR_OUT="results/pi"
BAM_LIST="$DIR_OUT"/"$POP".bams

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

REF="../data/honeybee_genome/Amel_4.5_scaffolds.fa"

echo Calculating pi for POP: "$POP"

# make directory to store output (if doesn't yet exist)
mkdir -p "$DIR_OUT"

echo "finding site allele frequencies"
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) Use all sites to estimate site allele frequency
# based on this FAQ, we don't do folding at this stage yet: https://github.com/ANGSD/angsd/issues/259
angsd -out "$DIR_OUT/$POP" \
-anc "$REF" \
-fold 0 \
-underFlowProtect 1 \
-bam "$BAM_LIST" \
-remove_bads 1 -minMapQ 30 -minQ 20 \
-GL 1 \
-doSaf 1 \
-P 2 # cores

echo "done with SAF! calculating SFS"
# try folding here?
realSFS "$DIR_OUT/$POP.saf.idx" -P 2 -fold 1 > "$DIR_OUT/$POP.sfs"

echo "done with SFS! calculating within-pop diversity 'thetas'"
# try folding here?
angsd -out "$DIR_OUT/$POP" \
-anc "$REF" \
-doThetas 1 \
-doSaf 1 \
-fold 1 \
-pest "$DIR_OUT/$POP.sfs" \
-underFlowProtect 1 \
-bam "$BAM_LIST" \
-remove_bads 1 -minMapQ 30 -minQ 20 \
-GL 1 \
-P 2

echo "summarizing thetas for region and windows" # could also calculate pi directly from sfs
thetaStat do_stat "$DIR_OUT"/"$POP".thetas.idx -outnames "$DIR_OUT"/"$POP".thetasAll
# ignore windows for now
#thetaStat do_stat "$DIR_OUT"/"$POP".thetas.idx -win 5000 -step 1000 -outnames "$DIR_OUT"/"$POP".thetasWindows


# options
# basic quality filtering for reads
# anc polarizes SFS by reference genome
# underFlowProtect is necessary for large #s of bams
