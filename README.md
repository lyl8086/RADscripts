PREREQUISITES
=
CAP3 available from http://seq.cs.iastate.edu/cap3.html

VelvetOptimiser https://github.com/tseemann/VelvetOptimiser

Perl module "Parallel::ForkManager"

Linux split cut sed etc.

Using "sort_read_pairs.pl" to export paired or single fasta file for each locus, see http://catchenlab.life.illinois.edu/stacks/pe_tut.php

Usage
=
 

Notation：please stop previous processes of CAP3 before running this script, same for the latter shell scripts. You can use such cmd:

	 % ps -A |grep "cap3" |grep -v "grep"|xargs kill -9


CP3_Opti.sh
=
Optimised script for local assembly of paired reads for RAD, using two steps.
	
Options:
	
	% CP3_Opti.sh [in_path] [out_path] [consensus] [threads]
	
	[in_path] directory contains multi fasta files of the second reads (Illumina read1) grouped by the first reads (Illumina read2).
	
	[out_path] Out path for the assembly files. Also contain assembled contigs for each locus.

	[consensus] refers to the file of consensus for the first reads, please reverse compliment it !!!
	You can use my Perl script "exract_consensus.pl" to extact consensus sequence from Stacks' results.
	
	[threads] number of threads.
	
VelvetOpti.sh
=	
A script for local assembly of RAD paired reads using VelvetOptimiser for multi-files.
	
Options:
	
	% VelvetOpti.sh [in_path] [out_path] [threads]
	
	[in_path] directory contains multi fasta files of the paired reads grouped by the first reads.
	
	[out_path] Out path for the assembly files.
	
	[threads] number of threads.

CP3_Opti.pl
=	
Optimised script for local assembly of paired reads for RAD by CAP3, using two steps. Perl version, should be faster than the shell script. Please install the Perl module "Parallel::ForkManager".

Options:
	
	% CP3_Opti.pl [path to file] [out] [loci] [threads] [collect]
	
	[path to file] directory contains multi fasta files of the second reads (Illunina read1) grouped by the first reads (Illumina read2).
	
	[out] out path for the assembly files.
	
	[loci] consensus file, see above for details.
	
	[threads] number of threads.
	
	[collect] turn on (any non-null character, such as "y") this option if you want to collect the assembled contigs for the second reads.

