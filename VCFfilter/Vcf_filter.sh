#!/bin/bash
## VCF filter V0.3
## Author: Yulong Li.

#f_vcf=$1
#cov=$2
#minDP=$3
#maxDP=$4
#GQ=$5
#prefix=$6
#Fis=$7
#GMAF=$8
#LMAF=$9
####### Get options.
Usage='Vcf_filterGQ.sh [-i vcf file] [-c coverage] [-m minDepth] [-M maxDP] [-q minGQ] [-o out] [-f Fis] [-G global maf] [-L local maf]'

while getopts "i:c:m:M:Qq:o:f:G:L:gh" arg
do
	case $arg in
		i)
			f_vcf=$OPTARG  	#input file.
			;;
		c)
			cov=$OPTARG     #coverage.
			;;
		m)
			minDP=$OPTARG   #minimum Depth for a genotype.
			;;
		M)
			maxDP=$OPTARG   #minimum Depth for a genotype.
			;;
		Q)
			minQ="--minQ 30"           #minimum Quality for a site.
			;;
		q)
			if [ -n $OPTARG ];
			then
				GQ="--minGQ $OPTARG"   #minimum Quality for a genotype.
			fi
			;;
		o)
			prefix=$OPTARG  #prefix of out file.
			;;
		f)
			Fis=$OPTARG     #Fis.
			;;
		G)
			GMAF=$OPTARG    #global maf.
			;;
		L)
			LMAF=$OPTARG    #local maf.
			;;
		g)
			gflag=1         #global filtering for Fis and Ho.
			;;
		h)
			echo "$Usage"
			exit 1
			;;
		?)
			echo "unknown argument"
			exit 1
			;;
	esac
done

####### Test part.
if [ $# -lt 6 ] || [ $1 = -h ];
then 
	echo "$Usage"
	exit 1
fi
if ( ! ls *.pop.txt) >/dev/null 2>&1;
then
	echo "Pop files does not exist! Please rename as <name>.pop.txt"
	exit 1
fi

#######Test if files already exist.
if [ -a ${prefix:=rmHom}.vcf.gz ];then echo -e "${prefix:=rmHom}.vcf.gz file exists! Are you sure to overwrite it??";exit 1;fi
if [ -a black.lst ];then echo "black.lst file exists! Please remove it!!!";exit 1;fi  
if [ -a black.2nd ];then echo "black.2nd file exists! Please remove it!!!";exit 1;fi 
#######

pop=(`ls *.pop.txt`) #Population files including individuals' names
len=${#pop[@]}

####### Statistics Function
function vcf_stat()
{
	for ((i=0; i<$len; i++));
	do vcftools --gzvcf $1 --keep ${pop[i]} --hardy  -c |cut -f1,2,3 |sed 's/\//\|/g' >${pop[i]}.het.tmp
	vcftools --gzvcf $1 --keep ${pop[i]} --site-pi  -c |cut -f3 >${pop[i]}.pi.tmp
	cut -f3 ${pop[i]}.het.tmp|awk -F '|' '{if (NR > 1) {n=$1+$2+$3;Ho=$2/n;p=(2*$1+$2)/(2*n); \
	He=2*p*(1-p); if (He==0) {fis="NA"} else {fis=1-(Ho/He)}; \
	print Ho"\t"He"\t"fis} else {print "Ho""\t""He""\t""Fis"}}' \
	|paste  ${pop[i]}.het.tmp - ${pop[i]}.pi.tmp>${pop[i]%".pop.txt"}.het
	done
	rm *.tmp
}
####### End

####### Calculate maf
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
####### end

####### Test if vcf file is gzipped.
if [ ${f_vcf:0-2} == 'gz' ];
then
	in_file="--gzvcf ${f_vcf}"
elif [ ${f_vcf:0-3} == 'vcf' ];
then 
	in_file="--vcf ${f_vcf}"
else
	echo -e "Input file ${f_vcf} is not vcf format?"
	exit 1
fi
#######

####### Output parameters.
echo -e "
File name  : $f_vcf
Coverage   : $cov
Missing    : $cov
MinGQ      : $GQ
MinDepth   : $minDP
MaxDP      : $maxDP
MaxHo      : 0.5
MinFis     : -$Fis
Global MAF: $GMAF
Local MAF  : $LMAF
Out name : $prefix.vcf.gz" >${prefix:=rmHom}.Parameters.log

####### Fist filter
echo "First filtering..."
vcftools $in_file --max-missing $cov --minDP $minDP --maxDP $maxDP $minQ $GQ --min-alleles 2 --max-alleles 2 \
--remove-indels --remove-filtered-all --remove-filtered-geno-all --recode -c --out 1st.filter |bgzip -c >tmp.vcf.gz
#######

####### Coverage filter
echo "Filtering coverage..."
for ((i=0; i<$len; i++));
do vcftools --gzvcf tmp.vcf.gz --keep ${pop[i]} --max-missing $cov --removed-sites -c >>black.lst
done
vcftools --gzvcf tmp.vcf.gz --exclude-positions black.lst --recode -c --out cov.filter|bgzip -c >vcf.cov.gz
#######

####### Remove Ho gt 0.5
echo "Filtering Ho..."
if [ gflag==1 ];
then
	#Calculating statistics on all individuals.
	vcftools --gzvcf vcf.cov.gz --hardy  -c |cut -f1,2,3 |sed 's/\//\|/g' >all.het.tmp
	cut -f3 all.het.tmp|awk -F '|' '{if (NR > 1) {n=$1+$2+$3;Ho=$2/n;p=(2*$1+$2)/(2*n); \
	He=2*p*(1-p); \
	if (He==0) {fis="NA"} else {fis=1-(Ho/He)}; \
	print Ho"\t"He"\t"fis} else {print "Ho""\t""He""\t""Fis"}}' \
	|paste  all.het.tmp - >all.het
else
	#Calculating statistics on each population.
	vcf_stat vcf.cov.gz
fi

for f in *.het;do awk '$4-0.5>0 {print $0}' $f | cut -f1,2 >>black.2nd ##exclude HO>0.5;
done
#######

if [ -n "$Fis" ];
then
	echo "Filtering Fis..."
	awk 'sqrt($6*$6)+0>'$Fis'+0 {print $0}' *.het|cut -f1,2 >>black.2nd ##exclude Fis < - $fis or Fis > $fis;
fi

vcftools --gzvcf vcf.cov.gz --exclude-positions black.2nd --recode -c --out HET.filter|bgzip -c >het.vcf.gz
rm *.het vcf.cov.gz
#######

####### Filter MAF
if [ -n "$LMAF" ];
then
	echo "Filtering global and local MAF..."
	vcftools --gzvcf het.vcf.gz --maf $GMAF --removed-sites -c >maf.removed #Get sites of maf  
	vcftools --gzvcf het.vcf.gz --positions maf.removed --recode -c |bgzip -c >maf.vcf.gz #Get the vcf file of maf
	vcf_stat maf.vcf.gz
	calc_maf $LMAF
	vcftools --gzvcf het.vcf.gz --exclude-positions black.maf --recode -c --out MAF.filter|bgzip -c >${prefix:=rmHom}.vcf.gz
	rm het.vcf.gz maf.vcf.gz *.het
elif [ -n "$GMAF" ];
then 
	echo "Filtering global MAF..."
	vcftools --gzvcf het.vcf.gz --maf $GMAF --recode -c --out MAF.filter |bgzip -c >${prefix:=rmHom}.vcf.gz
else	
	mv het.vcf.gz ${prefix:=rmHom}.vcf.gz
fi
#######

####### Output statistics files
echo "Outputting Statistics..."
vcf_stat ${prefix:=rmHom}.vcf.gz
echo "Final retained `zgrep -v '^$\|^#' ${prefix:=rmHom}.vcf.gz|wc -l` sites ^_^" 
#######
rm *tmp*
mkdir Log Het
mv *.log *black* *cnt *removed ./Log
mv *.het ./Het
echo "Done!"