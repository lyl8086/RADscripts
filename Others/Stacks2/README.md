A modified version of [Stacks2](http://catchenlab.life.illinois.edu/stacks/stacks_v2.php)
=
modifications (experimental features): 

1. support more fastq header formats [tsv2bam].

 * type I: end with 1 or 2 seperated by a delimator (|, / or _).
 
   e.g. `@4_1101_13393_18/1`

 * type II: seperated by space, containing index informations. 

   e.g. `@HWI-ST1276:71:C1162ACXX:1:1101:1208:2458 1:N:0:CGATGT`

2. support the minimum depth (`-m`) option [populations].

 * apply the filter on each SNP for each loci of each sample.

   if all the GT of a sample for a loci are missing, delete it.
 
 * likelihood filter (`-c`) added [populations].

3. other minor bugs.

Note: `-r` and `-p` options are applied on loci, not snp...

How to install
---
```
tar -xvf stacks-2.0Beta7-modified.tar.gz
./configure
make
```

Copy the new compiled program to the PATH, then run!
