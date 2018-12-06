#!/bin/bash
rm -rf tmp_out
mkdir tmp_out
threads=$1
if [ ! "$threads" ]; then threads=`grep -c 'processor' /proc/cpuinfo`; fi
v=`ustacks --version 2>&1 |head -1|awk '{print $2}'`
s=0
f=0
t=0
# test choose M
let t++
echo Testing functions for `RADassembler -v`... | tee tmp_out/logs 
echo "Stacks version: $v" | tee -a tmp_out/logs 
echo
echo "  $t. Testing chooseM function for RADassembler..." | tee -a tmp_out/logs 
t0=$(date +%s)
t1=$(date +%s)

RADassembler chooseM samples/read1/fishsim-1.fq.gz tmp_out/chooseM 10 $threads 3 \
>>tmp_out/logs 2>&1

if [ $? == 0 ]; then let s++; t2=$(date +%s)
echo "     chooseM successful, time $(( $t2 - $t1 ))s, see tmp_out/chooseM/."
else let f++; echo "     failed in chooseM function!"; fi

# test chooseN
let t++
echo "  $t. Testing single run function for RADassembler..."  | tee -a tmp_out/logs 
t1=$(date +%s)

RADassembler -i samples/read1 -s samples/read2 -P samples/PopMap -m 6 -t $threads \
-o tmp_out/chooseN -R ustacks >>tmp_out/logs 2>&1

if [ $? == 0 ]; then let s++; t2=$(date +%s)
echo "     Single run of Stacks successful, time $(( $t2 - $t1 ))s, see tmp_out/chooseN/stacks/."
else let f++; echo "     failed to run Stacks on each individual!"; fi

let t++
echo "  $t. Testing chooseN function for RADassembler..." | tee -a tmp_out/logs
t1=$(date +%s)

RADassembler chooseN tmp_out/chooseN/stacks tmp_out/chooseN samples/PopMap 10 $threads \
>>tmp_out/logs 2>&1

if [ $? == 0 ]; then let s++; t2=$(date +%s) 
echo "     chooseN successful, time $(( $t2 - $t1 ))s, see tmp_out/chooseN/."
else let f++; echo "     failed in chooseN function!"; fi

# test r80
let t++
echo "  $t. Testing r80 function for RADassembler..." | tee -a tmp_out/logs
t1=$(date +%s)

RADassembler r80 samples/read1/ tmp_out/r80_out/ samples/PopMap 10 10 $threads \
>>tmp_out/logs 2>&1

if [ $? == 0 ]; then let s++; t2=$(date +%s) 
echo "     r80 successful, time $(( $t2 - $t1 ))s, see tmp_out/r80_out/."
else let f++; echo "     failed in r80 function!"; fi

# test full assembly
let t++
echo "  $t. Testing full de novo assembly function for RADassembler..." | tee -a tmp_out/logs
t1=$(date +%s)

RADassembler -i samples/read1 -s samples/read2 -P samples/PopMap -m 6 -n 6 -t $threads \
-o tmp_out/Assembly_out -D 10:100 >>tmp_out/logs 2>&1

if [ $? == 0 ]; then let s++; t2=$(date +%s)
echo "     Full de novo assembly successful, time $(( $t2 - $t1 ))s, see tmp_out/Assembly_out/."
else let f++; t2=$(date +%s); echo "     failed in full de novo assembly!"; fi
echo "==============================================================="
echo "Total $t runs, $s success, $f fail, time $(( $t2 - $t0 ))s."
echo "All done!"
