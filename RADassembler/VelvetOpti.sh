#!/bin/bash
##VelvetOptimiser for multi files
##Author: Yulong Li

f_in=$1
f_out=$2
T=$3	##number of threads according to number of cpus
if [ $# -lt 3 ] || [ $1 = -h ];
then 
	echo "Usage: VelvetOpti.sh <in_path> <out_path> <threads>"
	exit 1
fi

if [ -d ${f_out} ];
then 
	echo -e "Dir ${f_out} exist, please remove it!!!"
	exit 1
else
	mkdir ${f_out}
fi

##split fasta files into several groups 
cd ${f_in}
ls -f *.fa >fasta.txt
cd -
mv -i ${f_in}/fasta.txt ./

split -d -n l/$T fasta.txt F_vel

## start the iteration
for line in F_vel*;
do (a=(`cat $line`); 
	for ((i=0; i<${#a[@]}; i++)); 
		do (VelvetOptimiser.pl -s 19 -e 55 -x 4 -k 'max' -c 'max' \
			-p ${f_out}/${a[i]}.out -f "-fasta -shortPaired ${f_in}/${a[i]}" &>/dev/null;

##Test if the assembly is successful.
#count=`find ${f_out} -name "${a[i]}.out_*" -type d|wc -l`; #Deprecated
#echo -e "${a[i]}\t$count"

			if [ $?==0 ];
			then
				sed "s/^>/>${a[i]}./g" ${f_out}/${a[i]}.out*/*.fa >>${f_out}/$line.collected.fa
				rm -rf ${f_out}/${a[i]}.out*
			else
				rm -rf ${f_out}/${a[i]}.out*
			fi
			)
		done
	) &
done
##end

wait
##collect the final contigs
cd ${f_out}
cat *.collected.fa >velvetopti.final.fa
rm -rf *.collected.fa
cd - && rm -rf F_vel* fasta.txt
echo "Done!!!"
