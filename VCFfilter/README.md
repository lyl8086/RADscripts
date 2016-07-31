VCFfilter a script for filtering VCF files suit for RAD data.

Filter options include coverage for each population, depth for each individual, Ho, MAF and Fis etc.  

Prerequisites:
=
#VCFTOOLS https://github.com/vcftools/vcftools

#bgzip https://github.com/samtools/htslib


Usage:
=
Vcf_filter.sh <gzvcf file> <coverage> <minDepth> <Fis> <maf> <out>
