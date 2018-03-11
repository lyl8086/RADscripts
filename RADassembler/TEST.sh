#!/bin/bash
rm -rf tmp_out
mkdir tmp_out
threads=$1
if [ ! "$threads" ]; then threads=`grep -c 'processor' /proc/cpuinfo`; fi
v=`ustacks --version 2>&1 |head -1|awk '{print $2}'`
s=0
f=0
# test choose M
echo "Testing functions for RADassembler..."
echo "Stacks version: $v"
echo
echo "  1. Testing chooseM function for RADassembler..."
t0=$(date +%s)
t1=$(date +%s)

RADassembler chooseM samples/read1/fishsim-1.fq.gz tmp_out/chooseM 10 $threads 5 \
>tmp_out/logs 2>&1

if [ $? == 0 ]; then let s++; t2=$(date +%s)
echo "     chooseM successful, time $(( $t2 - $t1 ))s, see tmp_out/chooseM."
else let f++; echo "     failed in chooseM function!"; fi

# test chooseN
echo "  2. Testing single run function for RADassembler..."
t1=$(date +%s)

RADassembler -i samples/read1 -s samples/read2 -P samples/PopMap -M 5 -m 5 -t $threads \
-o tmp_out/chooseN -R ustacks >>tmp_out/logs 2>&1

if [ $? == 0 ]; then let s++; t2=$(date +%s)
echo "     Single run of Stacks successful, time $(( $t2 - $t1 ))s, see tmp_out/chooseN/stacks."
else let f++; echo "     failed to run Stacks on each individual!"; fi

echo "  3. Testing chooseN function for RADassembler..."
t1=$(date +%s)

RADassembler chooseN tmp_out/chooseN/stacks tmp_out/chooseN samples/PopMap 10 $threads \
>>tmp_out/logs 2>&1

if [ $? == 0 ]; then let s++; t2=$(date +%s) 
echo "     chooseN successful, time $(( $t2 - $t1 ))s, see tmp_out/chooseN."
else let f++; echo "     failed in chooseN function!"; fi

# test full assembly
echo "  4. Testing full de novo assembly function for RADassembler..."
t1=$(date +%s)

RADassembler -i samples/read1 -s samples/read2 -P samples/PopMap -M 5 -m 5 -n 5 -t $threads \
-o tmp_out/Assembly_out -D 10:100 >>tmp_out/logs 2>&1

if [ $? == 0 ]; then let s++; t2=$(date +%s)
echo "     Full de novo assembly successful, time $(( $t2 - $t1 ))s, see tmp_out/Assembly_out."
else let f++; t2=$(date +%s); echo "     failed in full de novo assembly!"; fi
echo "==============================================================="
echo "Total 4 runs, $s success, $f fail, time $(( $t2 - $t0 ))s."
echo "All done!"
