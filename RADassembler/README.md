RADAssembler
===

<b>A Pipeline For Assembly of RAD-seq (RPE) from Multiple Indidviduals</b>

Note: Only for Paired-end RAD-seq reads with random sheared ends (the original RAD protocol).

PREREQUISITES
---
* [CAP3](http://seq.cs.iastate.edu/cap3.html)

* [STACKS](http://catchenlab.life.illinois.edu/stacks/)

* Perl module [Parallel::ForkManager](https://metacpan.org/release/Parallel-ForkManager)

How to Install
---
```
# wget https://github.com/lyl8086/RADscripts/raw/master/RADassembler/RADassembler.tar.gz
# tar -xvf RADassembler.tar.gz && cd RADassembler && bash INSTALL.sh
```
How to Run
---
important:

`-i` individuals' files (reads with enzyme cut sites) should be name as <b>`name.fq.gz`</b>

`-s` individuals' files (with randomly sheared reads) should be name as <b>`name_2.fq.gz` or `name_1.fq.gz`</b>

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
    -c: individual coverage for a locus.
    -A: turn off assembly.
    -R: run a single component.
    -t: number of threads.
    chooseM: Similarity threshold selection within individual [ustacks].
    chooseN: Similarity threshold selection across individuals [cstacks].
```
How to select similarity thresholds within and across individuals
---
* within individual [ustacks]
```
Usage: RADassembler chooseM [infile] [outpath] [max mismatch] [threads] [minDP] [replot] [yrange]
```
* across individuals[cstacks]
```
Usage: RADassembler chooseN [inpath] [outpath] [popmap] [max mismatch] [threads] [replot] [yrange]
```
<b>see</b> [samples](samples) for an example run.

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

