#!/bin/bash
# the master nucleotide diversity calculation script
# run with parallel 
# parallel \
# --noswap --joblog logs/my_log.log --jobs 4 \
# ./my_script.sh {2}.list {1} :::: all_scaffolds.list :::: all_pops.list

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# define your variables
POP=$1 # first argument is a list of BAMS
SCAFFOLD=$2
minInd=13
OUT=${POP}.${SCAFFOLD}.out

# reference genome
REF=/media/data/dazarate/European_genome_scaffolds/Amel_4.5_scaffolds.fa


# Step 1: Finding a 'global estimate' of the SFS
# -doSaf 1 : Calculate the Site allele frequency likelihood based on individual genotype likelihoods assuming HWE
#  for now, we can use the REF as ANC...ancestral state needs to be supplied for the full SFS
## but you can use the -fold 1 to estimate the folded SFS and then use the reference as ancestral.

angsd.v928 \
-bam $POP \
-r $SCAFFOLD \
-doSaf 1 \
-anc $REF \
-GL 2 \
-P 10 \
-out $OUT \
-minQ 20 \
-minMapQ 30 \
-minInd $minInd

# Step 2: Using realSFS, fold the SAF to create the SFS
echo "Obtaining maximum likelihood estimate of the SFS using realSFS"
INPUT=${BAM}.${SCAFFOLD}.out.saf.idx
singularity exec /media/data/dazarate/singularity.images/angsd.v928.sif /opt/angsd/angsd/misc/realSFS $INPUT -fold 1 -P 10 > $OUT.sfs


#Step 3: Calculate the thetas for each site
angsd.v928 -bam $POP -out $OUT -r $SCAFFOLD: -doThetas 1 -doSaf 1 -pest $OUT.sfs  -anc $REF -GL 2
