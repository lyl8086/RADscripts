# RADassembler

<b>A Pipeline For Assembly of RAD-seq (RPE) from Multiple Individuals</b>

Note: Only for Paired-end RAD-seq reads with random sheared ends (the original RAD protocol).

## PREREQUISITES

* [CAP3](http://seq.cs.iastate.edu/cap3.html)

* [STACKS](http://catchenlab.life.illinois.edu/stacks/)

* Perl >= 5.16 and modules including
  
  [Parallel::ForkManager](https://metacpan.org/release/Parallel-ForkManager)
  
  [Array::Shuffle](https://metacpan.org/release/Array-Shuffle)

## How to Install

```
# wget https://github.com/lyl8086/RADscripts/releases/download/V1.11/RADassembler.tar.gz
# tar -xvf RADassembler.tar.gz
# cd RADassembler
# bash INSTALL.sh
# if you want to run a test, run bash TEST.sh
```

## How to Run
important:

`-i` individuals' files (reads with enzyme cut sites) should be name as <b>`name.fq.gz`</b>, [examples](samples/read1)

`-s` individuals' files (with randomly sheared reads) should be name as <b>`name_2.fq.gz` or `name_1.fq.gz`</b>, [examples](samples/read2)

PopMap [example](samples/PopMap)

```
# RADassembler

Usage:
    -i: path to reads with enzyme cut sites of each individual, i.e. read1 files? 
        name as (individual name).fq[fa][.gz].
    -o: path to output.
    -s: path to paired-end reads of each individual, i.e. read2 files?
        name as (individual name_[12]).fq[fa][.gz].
    -P: PopMaP file.
    -f: type of input files. "fasta", "fastq", "gzfastq", "gzfasta".
    -M: minimum stack depth.
    -D: minimum read depth of a locus to export for assembly, also accept [lower:upper].
    -m: mismatch for ustacks.
    -n: mismatch for cstacks.
    -g: perform gapped assembly in stacks, 1 for ustacks, 2 for cstacks, 3 for all, 0 turn off.
    -c: individual coverage for a locus.
    -a: assembler, either "cap3" or "velvet", default cap3.
    -k: hash length for velvet, default 27.
    -x: do not verify multiple matches when exporting fasta files.
    -A: turn off assembly.
    -R: run a single component, accept "ustacks", "cstacks" or "assembly".
    -t: number of threads.
    
    chooseM: Similarity threshold selection within individual [ustacks].
    chooseN: Similarity threshold selection across individuals [cstacks].

```

## How to select similarity thresholds within and across individuals

* within individual(ustacks)

```
Usage : RADassembler chooseM [infile] [outpath] [max mismatch] [threads] [minDP] [gapped] [replot] [yrange]

infile : a single individual reads file with enzyme cut site for clustering.
outpath : output path for this run.
max mismatch : maximum mismatch to run ustacks, RADassembler will run a set of mismatches by a step of one.
threads : number of threads used for ustacks.
minDP : minimum stack depth for ustacks.
gapped : gapped assembly for ustacks (1: on, 0: off).
replot : just re-plot graphs using different parameters.
yrange : ranges for y-axis used in plot, minimum:maximum.
```

* across individuals(cstacks)

```
Usage: RADassembler chooseN [inpath] [outpath] [popmap] [max mismatch] [threads] [gapped] [replot] [yrange]

inpath: input path for cstacks, containing ustacks files.
outpath: output path for this run.
popmap: a [PopMap](samples/PopMaP) file, e.g. `individual\tpopulation`.
max mismatch: maximum mismatch for cstacks, RADassembler will run a set of mismatches, i.e. 1,2...maximum.
threads: number of threads used in cstacks.
gapped: turn on gapped assembly? 1:on, 0:off.
replot: just re-plot graphs.
yrange: ranges for y-axis used in plot, minimum:maximum.
```

<b>see</b> [samples](samples) for an example run.

## Tutorial

[Tutorial](Tutorial.md)

## About the output

Folder [Assembly](samples/Assembly_out/Assembly)

* <b>collected_final.fa</b>: the final assembled contigs, with enzyme cut sites at 3' end.
* <b>assembly_1st</b>: folder contains fasta files of the first assembly
* <b>assembly_2nd</b>: folder contains fasta files of the second assembly
* <b>log</b>: folder contains run paramerts

Folder [stacks](samples/Assembly_out/stacks)

* results of stacks runs

Folder [reads_export](samples/Assembly_out/reads_export)

* exported fasta files for assembly.

## Contact

Yulong Li <liyulong12@mails.ucas.ac.cn>

Please consider to cite our paper:

> Li YL, Xue DX, Zhang BD, Liu JX. 2018 An optimized approach for local de novo assembly of overlapping 
paired-end RAD reads from multiple individuals. Royal Society Open Science. 5:171589. [[link]](http://dx.doi.org/10.1098/rsos.171589)
