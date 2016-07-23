#!/bin/bash
## VCF filter V0.3
##Usage: Vcf_filterV3.sh <gzvcf file> <coverage> <minDepth> <Fis> <out>
f_vcf=$1
cov=$2
dep=$3
Fis=$4
MAP=$5
prefix=$6
if [ $# -lt 3 ] || [ $1 = -h ];
then 
	echo "Usage: Vcf_filterV3.sh <gzvcf file> <coverage> <minDepth> <Fis> <mac> <out>"
	exit 1
fi
if ( ! ls *.pop.txt) >/dev/null 2>&1;
then
	echo "Pop files does not exist! Please rename as <name>.pop.txt"
	exit 1
fi

echo -e "
File name: $f_vcf
Coverage : $cov
Missing  : 0.9
minGQ    : 20
MinDepth : $dep
MaxHo    : 0.5
MinFis   : -$Fis
MAF      : $MAP
Out name : $prefix.vcf.gz" >$prefix.Parameters.log

pop=(`ls *.pop.txt`)                     #Population files including individuals' names
len=${#pop[@]}
##Statistics Function
function vcf_stat()
{
	for ((i=0; i<$len; i++));
	do vcftools --gzvcf $1 --keep ${pop[i]} --hardy  -c |cut -f1,2,3 |sed 's/\//\|/g' >${pop[i]}.het.tmp
	vcftools --gzvcf $1 --keep ${pop[i]} --site-pi  -c |cut -f3 >${pop[i]}.pi.tmp
	cut -f3 ${pop[i]}.het.tmp|awk -F '|' '{if (NR > 1) {n=$1+$2+$3;Ho=$2/n;p=(2*$1+$2)/(2*n); \
	He=(n/(n-1))*((2*p*(1-p))-(Ho/(2*n))); \
	if (He==0) {fis="NA"} else {fis=1-(Ho/He)}; \
	print Ho"\t"He"\t"fis} else {print "Ho""\t""Nei_He""\t""Fis"}}' \
	|paste  ${pop[i]}.het.tmp - ${pop[i]}.pi.tmp>${pop[i]%".pop.txt"}.het
	done
	rm *.tmp
}
## End

##Calculate maf
function calc_maf()
{
	
	for f in *.het;	#Providing Statistics files firstly.
	do cut -f3 $f|awk -F '|' '{if (NR > 1) {a1=2*$1+$2; a2=2*$3+$2; n=$1+$2+$3; if (a1<a2) {maf=a1/(2*n);flag=1} else {maf=a2/(2*n);flag=2}; \
	if (maf<'$1') {print flag} else {print "NA"}} else {print "'${f%".het"}'"}}' >${f%".het"}.maf.tmp
	done
	
	cut -f1,2 $f |paste - *.maf.tmp >final.MAF.tmp
	awk '{if (NR==1) {print $0"\tallele1\tallele2"} else {c=3;A=0;B=0;while (c<=NF) {if ($c==1) {A=A+1} else if ($c==2) \
	{B=B+1};c++} print $0"\t"A"\t"B}}' final.MAF.tmp >maf.cnt
	awk '$'$(($len+3))'== '$len' || $'$(($len+4))'== '$len' {print $0}' maf.cnt |cut -f1,2 >black.maf
	rm *.tmp
	
}
##end

## Fist filter
echo "First filtering..."
vcftools --gzvcf ${f_vcf} --max-missing 0.9 --minDP $dep --minQ 30 --minGQ 15 --min-alleles 2 --max-alleles 2 \
--recode --recode-INFO-all -c --out 1st.filter |bgzip -c >tmp.vcf.gz
##

## Coverage filter
echo "Filtering coverage..."
for ((i=0; i<$len; i++));
do vcftools --gzvcf tmp.vcf.gz --keep ${pop[i]} --max-missing $cov --removed-sites -c >>black.lst
done
vcftools --gzvcf tmp.vcf.gz --exclude-positions black.lst --recode --recode-INFO-all -c --out 2nd.filter|bgzip -c >vcf.cov.gz
##

## Remove Ho gt 0.5
echo "Filtering Ho..."
vcf_stat vcf.cov.gz
for f in *.het;do awk '$4-0.5>0 {print $0}' $f | cut -f1,2 >>black.2nd ##exclude HO>0.5 modified 
done
if [ -n "$Fis" ];
then
	echo "Filtering Fis..."
	awk '$6+'$Fis'<0 {print $0}' *.het|cut -f1,2 >>black.2nd ##exclude abs(Fis) modified
fi

if [ -n "$MAP" ];
then
	echo "Filtering MAF..."
	calc_maf $MAP
	cat black.maf >>black.2nd
fi

vcftools --gzvcf vcf.cov.gz --exclude-positions black.2nd --recode --recode-INFO-all -c --out 3rd.filter|bgzip -c >${prefix:=rmHom}.vcf.gz
rm *.het vcf.cov.gz
###

## Output statistics files
echo "Outputting Statistics..."
vcf_stat ${prefix:=rmHom}.vcf.gz
echo "Final retained `zgrep -v '^$\|^#' ${prefix:=rmHom}.vcf.gz|wc -l` sites ^_^" 
##
rm *tmp*
mkdir Log Het
mv *.log *black* ./Log
mv *.het ./Het
echo "Done!"