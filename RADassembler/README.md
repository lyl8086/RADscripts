<h1>RADAssembler</h1>

<b>A Pipeline For Assembly of RAD-seq from Multiple Indidviduals</b>

---


PREREQUISITES
---
* [CAP3](http://seq.cs.iastate.edu/cap3.html)

* [STACKS](http://catchenlab.life.illinois.edu/stacks/)

How to install
---
```
#sh INSTALL.sh
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
About the output
---
<ul>
<li><b>collected_final.fa</b>: the final assembled contigs</li>
<li><b>assembly_1st</b>: folder contains fasta files of the first assembly</li>
<li><b>assembly_2nd</b>: folder contains fasta files of the second assembly</li>
<li><b>log</b>: folder contains run paramerts</li>
</ul>


