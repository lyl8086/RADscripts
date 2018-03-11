RADassembler
===

<b>A Pipeline For Assembly of RAD-seq (RPE) from Multiple Individuals</b>

Note: Only for Paired-end RAD-seq reads with random sheared ends (the original RAD protocol).

PREREQUISITES
---
* [CAP3](http://seq.cs.iastate.edu/cap3.html)

* [STACKS](http://catchenlab.life.illinois.edu/stacks/)

* Perl module [Parallel::ForkManager](https://metacpan.org/release/Parallel-ForkManager)

How to Install
---
```
# wget https://github.com/lyl8086/RADscripts/releases/download/V1.01/RADassembler.tar.gz
# tar -xvf RADassembler.tar.gz
# cd RADassembler
# bash INSTALL.sh
# if you want to run a test, run bash TEST.sh
```
How to Run
---
important:

`-i` individuals' files (reads with enzyme cut sites) should be name as <b>`name.fq.gz`</b>, [examples](samples/read1)

`-s` individuals' files (with randomly sheared reads) should be name as <b>`name_2.fq.gz` or `name_1.fq.gz`</b>, [examples](samples/read2)

PopMap [example](samples/PopMap)


```
# RADassembler

Usage:
    -i: input path. Clean reads containing enzyme site, read1?
        e.g. (individual name).fq[fa][.gz].
    -o: out path.
    -s: paired-end path. Paired-end reads, read2?
        e.g. (individual name_[12]).fq[fa][.gz].
    -f: input file type. "fasta", "fastq", "gzfastq", "gzfasta".
    -P: PopMaP file.
    -M: minimum stacks depth.
    -D: minimum read depth of a locus to export for assembly,
        can be [lower:upper].
    -m: mismatch for ustacks.
    -n: mismatch for cstacks.
    -g: perform gapped assembly in stacks, 1 for ustacks, 2 for cstacks, 3 for all.
    -c: individual coverage for a locus.
    -A: turn off assembly.
    -R: run a single component, accept "ustacks", "cstacks", "assembly".
    -t: number of threads.
    chooseM: Similarity threshold selection within individual [ustacks].
    chooseN: Similarity threshold selection across individuals [cstacks].
```
How to select similarity thresholds within and across individuals
---
* within individual [ustacks]
```
Usage : RADassembler chooseM [infile] [outpath] [max mismatch] [threads] [minDP] [gapped] [replot] [yrange]

infile : a single individual reads file with enzyme cut site for clustering.
outpath : output path for this run.
max mismatch : maximum mismatch to run ustacks, RADassembler will run a set of mismatches by a step of one.
threads : number of threads used for ustacks.
minDP : minimum stack depth for ustacks.
gapped : gapped assembly for ustacks (1: on, 0: off).
replot : just re-plot graphs using different parameters.
yrange : ranges for axis used in plot, minimum:maximum.
```

* across individuals[cstacks]
```
Usage: RADassembler chooseN [inpath] [outpath] [popmap] [max mismatch] [threads] [gapped] [replot] [yrange]

inpath: input path for cstacks, containing ustacks files.
outpath: output path for this run.
popmap: a [PopMap](samples/PopMaP) file, e.g. `individual\tpopulation`.
max mismatch: maximum mismatch for cstacks, RADassembler will run a set of mismatches, i.e. 1,2...maximum.
threads: number of threads used in cstacks.
gapped: turn on gapped assembly? 1:on, 0:off.
replot: just re-plot graphs.
yrange: ranges for axis used in plot, minimum:maximum.
```
<b>see</b> [samples](samples) for an example run.

Tutorial
---
Working on...

About the output
---
Folder [Assembly](samples/Assembly_out/Assembly)
* <b>collected_final.fa</b>: the final assembled contigs
* <b>assembly_1st</b>: folder contains fasta files of the first assembly
* <b>assembly_2nd</b>: folder contains fasta files of the second assembly
* <b>log</b>: folder contains run paramerts

Folder [stacks](samples/Assembly_out/stacks)
* results of stacks runs

Folder [reads_export](samples/Assembly_out/reads_export)
* exported fasta files for assembly.

