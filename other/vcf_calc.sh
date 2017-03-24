####### Read PopMap
f_vcf=$1
PopMap=$2
if [ $# -lt 2 ];
then
	echo "$0 [vcf] [popmap]"
	exit 1
fi
####### Test if vcf file is gzipped.
if [ ${f_vcf:0-2} == 'gz' ];
then
	in_file="--gzvcf ${f_vcf}"
elif [ ${f_vcf:0-3} == 'vcf' ];
then 
	in_file="--vcf ${f_vcf}"
elif [ ${f_vcf:0-2} == '7z' ];
then
	in_file="--vcf -"
	lzma='y'
else
	echo -e "Input file ${f_vcf} is not vcf format?"
	exit 1
fi
####################################
declare -A pop

while read line
do
	key=`echo $line | awk '{print $2}'`
	val=`echo $line | awk '{print $1}'`
	#echo $line
	pop[$key]+="--indv $val " #For vcftools.
done < $PopMap

len=${#pop[@]}
#######

####### Statistics Function
function vcf_stat()
{
	for i in ${!pop[@]}; #Change to Shell Associative Array.
	do vcftools $1 ${pop[$i]} --hardy  -c |cut -f1,2,3 |sed 's/\//\|/g' >$i.het.tmp
	#vcftools $1 ${pop[$i]} --site-pi  -c |cut -f3 >$i.pi.tmp
	cut -f3 $i.het.tmp|awk -F '|' '{if (NR > 1) {n=$1+$2+$3;Ho=$2/n;p=(2*$1+$2)/(2*n); \
	He=2*p*(1-p); if (He==0) {fis="NA"} else {fis=1-(Ho/He)}; \
	if (p<=0.5) {maf=p} else {maf=1-p};
	if (maf>0) {print Ho"\t"He"\t"fis"\t"maf}} else {print "Ho""\t""He""\t""Fis""\tMAF"}}' >$i.maf
	#|paste  $i.het.tmp - $i.pi.tmp>$i.het
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

vcf_stat "$in_file"
