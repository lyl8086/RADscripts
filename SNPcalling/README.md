<strong>Use bounded-error SNP calling model from [Stacks](http://catchenlab.life.illinois.edu/stacks/) to call SNPs from paired-end RAD data.</strong>

Please cite our paper if you find it useful to your work.

>Yulong Li <em>et al</em>. (in prep) An optimized highly efficient approach for local de novo assembly of multiple individuals’ reads based on overlapping paired-end RAD sequencing

PREREQUISITES
---
* [samtools](https://sourceforge.net/projects/samtools/files/samtools/0.1.19/)

* [bgzip](https://github.com/samtools/htslib)

* perl


Usage
---
```shell
  samtools mpileup -b [bam list] -f [reference] -B | Call_snp.pl -c [bam list] | bgzip -c >out.vcf.gz

  bam list: List of input BAM files, one file per line
```

Reference
---
J. Catchen, P. Hohenlohe, S. Bassham, A. Amores, and W. Cresko. Stacks: an analysis tool set for population genomics. Molecular Ecology. 2013. [[reprint]](http://dx.doi.org/10.1111/mec.12354)

