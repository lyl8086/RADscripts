#PREREQUISITES
---
* [CAP3](http://seq.cs.iastate.edu/cap3.html)

* [VelvetOptimiser](https://github.com/tseemann/VelvetOptimiser)

* Perl module "Parallel::ForkManager"

* Linux split cut sed etc.

> Using "sort_read_pairs.pl" to export paired or single fasta file for each locus, [see](http://catchenlab.life.illinois.edu/stacks/pe_tut.php)

#Usage
---
CP3_Opti.pl
---	
Optimal script for local assembly of RAD paired reads based on CAP3, using two steps. 

Options:
```
	CP3_Opti.pl [options]
	
	-i input path of fasta files
	-o output path for assembly files
	-l catalog file from stacks [batch_XXX.catalog.tags.tsv.gz]
	-t number of threads [1]
	-c collect the assembly results for the second reads [0]
	-p parameters parsed to CAP3
	-h help
```

Please cite our paper if you find it useful to your work.

>Yulong Li <em>et al</em>. (in prep) An optimized highly efficient approach for local de novo assembly of multiple individuals’ reads based on overlapping paired-end RAD sequencing.

