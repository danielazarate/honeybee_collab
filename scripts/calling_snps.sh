### SNP CALLING PIPELINE USING SAMTOOLS TO CALCULATE NUCLEOTIDE DIVERSITY
### Adapted from Dan Bolser " A Simple SNP Calling Pipeline" dbolser@ebi.ac.uk

## Version of Samtools: v1.10
## Version of VCFtools: v0.1.16

# trim the FASTQ reads for quality and length
## URS4
perl /projects/ps-dazarate/popoolation_1.2.2/basic-pipeline/trim-fastq.pl \
--input1 /projects/ps-dazarate/Wallberg.Ref.HoneyBees/scutellata.URS4.raw.reads.fastq \
--fastq-type sanger \
--output /oasis/tscc/scratch/dazarate/scutellata.URS4.trim \
--quality-threshold 25 --min-length 40

# Index reference genome in BWA and map the individual reads to reference 
module load bwa
bwa index -a is /projects/ps-dazarate/European_genome_scaffolds/Amel_4.5_scaffolds.fa

# Write as a BAM file (binary)
module load samtools 
samtools view -bS /oasis/tscc/scratch/dazarate/scutellata.B4.map > /oasis/tscc/scratch/dazarate/scutellata.B4.bam

# Remove low quality reads
module load samtools
samtools view -q 20 -bS /oasis/tscc/scratch/dazarate/scutellata.B4.bam > /oasis/tscc/scratch/dazarate/scutellata.B4.bam2

# Sort the BAM files individually 
samtools sort /oasis/tscc/scratch/dazarate/scutellata.B4.bam2 > /oasis/tscc/scratch/dazarate/scutellata.B4.bam3

# Index the genome assembly (again!)
samtools faidx my.fasta

# Run 'mpileup' to generate VCF format
samtools mpileup -g -f my.fasta my-sorted1.bam my-sorted-2.bam my-sorted-n.bam > myraw.bcf

# Call SNPs
bcftools view -bvcg my-raw.bcf > my-var.bcf

# Remove Duplicates with bcftools
bcftools norm --remove-duplicates input.bcf > fil.output.bcf

# Filter SNPs
bcftools view my.var.bcf | vcfutils.pl varFilter - > my.var-final.vcf

# Remove Duplicates with bcftools
bcftools norm --remove-duplicates input.bcf > fil.output.bcf

## Include only biallelic sites, sites with a MAF > or equal to 0.05, and sites with mean depth coverage with a max of X and a min of Y 
vcftools \
--vcf file1.vcf \
--min-alleles 2 \
--max-alleles 2 \
--maf 0.05 \
--min-meanDP <float> \
--max-meanDP <float>



# Use VCFtools to calculate nucleotide diversity
zcat input_file.vcf.gz | vcftools --vcf - --site-pi --positions SNP_list.txt --out nucleotide_diversity
