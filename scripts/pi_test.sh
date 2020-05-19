!/bin/bash
# this script calculates pi for a population and list of bams

# to run: ./pi_within_ancestry.sh AR01 R01.bams

POP="${1}"
BAM_LIST="${2}"
DIR_OUT="results/pi"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

REF="/media/data/dazarate/European_genome_scaffolds/Amel_4.5_scaffolds.fa"

echo Calculating pi for POP: "$POP"

# make directory to store output (if doesn't yet exist)
mkdir -p "$DIR_OUT"

echo "finding site allele frequencies"
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) Use all sites to estimate site allele frequency
# based on this FAQ, we don't do folding at this stage yet: https://github.com/ANGSD/angsd/issues/259
singularity run ~dazarate/scripts/angsd.v930.simg -out "$DIR_OUT/$POP" \
-r Group1.1 \
-anc "$REF" \
-remove_bads 1 \
-minMapQ 30 \
-minQ 20 \
-fold 0 \
-bam "$BAM_LIST" \
-GL 1 \
-doSaf 1 \
-underFlowProtect 1 \
-P 2 # cores

echo "done with SAF! calculating SFS"
# try folding here? ### supply reference genome as "anc" to polarize by reference allele 
singularity exec angsd.v930.simg /opt/angsd/angsd/misc/realSFS "$DIR_OUT/$POP.saf.idx" -P 2 -fold 1 > "$DIR_OUT/$POP.folded.sfs"

echo "done with SFS! calculating within-pop diversity 'thetas'"
# try folding here?
singularity exec angsd.v930.simg /opt/angsd/angsd/misc/realSFS  saf2theta "$DIR_OUT/$POP.saf.idx" -fold 1 -outname "$DIR_OUT/$POP"  -sfs "$DIR_OUT/$POP.folded.sfs" 

echo "summarizing thetas for region and windows" # could also calculate pi directly from sfs
singularity exec angsd.v930.simg /opt/angsd/angsd/misc/thetaStat do_stat "$DIR_OUT/$POP.thetas.idx" -outnames "$DIR_OUT/$POP.thetasAll"
# ignore windows for now
#thetaStat do_stat "$DIR_OUT"/"$POP".thetas.idx -win 5000 -step 1000 -outnames "$DIR_OUT"/"$POP".thetasWindows 


# options
# basic quality filtering for reads
# anc polarizes SFS by reference genome
# underFlowProtect is necessary for large #s of bams
