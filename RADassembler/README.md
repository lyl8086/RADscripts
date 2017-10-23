<h1>RADAssembler</h1>

<b>A Pipeline For Assembly of RAD-seq from Multiple Indidviduals</b>

---


PREREQUISITES
---
* [CAP3](http://seq.cs.iastate.edu/cap3.html)

* [STACKS](http://catchenlab.life.illinois.edu/stacks/)

* Perl module [Parallel::ForkManager](https://metacpan.org/release/Parallel-ForkManager)

How to Install
---
```
# sh INSTALL.sh
```
How to Run
---

```
# RADassembler

Usage:
  -i input path. Clean reads containing enzyme site, read1. e.g. (indiv name).fq[fa][.gz]
  -o out path. 
  -s paired-end path. Paired-end reads, read2. e.g. (indiv name_[12]).fq[fa][.gz]
  -f input file type. "fasta", "fastq", "gzfastq", "gzfasta"
  -P PopMaP file
  -M minimum stacks dapth
  -m mismatch for ustacks
  -n mismatch for cstacks
  -D minimum read depth of a locus to export for assembly
  -t number of threads
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

