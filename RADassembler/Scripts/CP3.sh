#!/bin/bash
##CAP3 for multi files
##Author: Yulong Li
In_path=$1
Out_path=$2
T=$3	##number of threads according to number of cpus
if [ $# -lt 3 ] || [ $1 = -h ];
then 
	echo "Usage: CP3.sh <in_path> <out_path> <threads>"
	exit 1
fi
## Test the out path, if not exit, create it.
if [ -d ${Out_path} ];
then 
	echo -e "Dir ${Out_path} exist, please remove it!!!"
	exit 1
else
	mkdir ${Out_path}
fi

##split fasta files into several groups 
cd ${In_path}
ls -f *.fa >fasta.txt
cd - && mv -i ${In_path}/fasta.txt ./
split -d -n l/$T fasta.txt F_cap

## start the iteration
for line in F_cap*;
do (a=(`cat $line`); 
	for ((i=0; i<${#a[@]}; i++)); 
		do (cap3 ${In_path}/${a[i]} &>/dev/null;	##do cap3 assembly on each locus
			sed "s/^>/>${a[i]}./g" ${In_path}/${a[i]}.cap.contigs >>${Out_path}/$line.collected.fa;	##collect contigs
			rm -rf ${In_path}/${a[i]}.cap.*; 
			)
		done
	) &
done
##	end
wait
##collect the final contigs
cat ${Out_path}/*.collected.fa >${Out_path}/cap3.final.fa
rm -rf ${Out_path}/*.collected.fa F_cap* fasta.txt
echo "Done!"
