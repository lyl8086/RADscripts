#!/usr/bin/env perl
##	Author:Yu-Long Li
use strict; 
use warnings; 
use Parallel::ForkManager;
use Getopt::Long;


my (@files, $file, $out_fh, $in_fh, $in_tags, $sub_lst);
my $in_path;
my $out_path;
my $loci;
my $num_threads = 1;
my $cmd = '';	##Parameters parsed to CAP3 for the fist assembly.
my $collect;
my $help = 0;
my $usage = 
'	
	Options [default values]:
	
	-i input path of fasta files
	-o output path for assembly files
	-l catalog file from stacks [batch_XXX.catalog.tags.tsv.gz]
	-t number of threads [1]
	-c collect the assembly results for the second reads [0]
	-p parameters parsed to CAP3
	-f sub file list for assembly 
	-h help
';

GetOptions ("par=s" 		=> \$cmd,
			"in=s" 			=> \$in_path,
			"out=s" 		=> \$out_path,
			"loci=s" 		=> \$loci,
			"threads=i" 	=> \$num_threads,
			"collect:1" 	=> \$collect,
			"file:s" 	    => \$sub_lst,
			"help:1"		=> \$help
			)
			or die ("Error in command line arguments\n");
	
die "$usage\n" if !($in_path && $out_path && $loci) or ($help);

## If the outpath does not exist, create it;

if (! -e "$out_path/assembly_1st") {
	
	`mkdir -p $out_path/assembly_1st`;
	`mkdir -p $out_path/log`
	
}

##delete any tmp files;

`find $in_path/ -name "*cap*" |xargs rm -rf`;


print "Starting the first run...\n";
unless($cmd) {
	#Parameters for short reads.
	$cmd = '-r 0 -i 30 -j 31 -o 18 -s 300 -p 85';
}

`echo -e "cap3 $cmd" > $out_path/log/assembly.1par`;


##	extract consensus from stacks catalogs;

if (substr($loci, -2) eq 'gz') {
	open($in_tags, "gunzip -c $loci|") or die "$!";
	} 
elsif (substr($loci, -2) eq '7z') {
	open($in_tags, "7zcat.sh $loci|") or die "$!";
	}
else {
	open($in_tags, "$loci") or die "$!";
}
my %seqs;
while (<$in_tags>) {
	
	next if /^#/;
	my @parts = split(/\t/, $_);
	my $id    = $parts[2];
	my $seq   = $parts[9];
	die "catalog file is wrong!!\n" if ($id+1 ne $.); 
	$seq = uc(reverse $seq);
	$seq =~ tr/ATCGN/TAGCN/;
	$seqs{$id} = $seq;
		
}

close $in_tags;

## if provided sub file list...
if ($sub_lst) {
	open(my $in_fh, "$sub_lst") or die "$!";
	while(<$in_fh>) {
		next if /^$/;
		chomp;
		push @files, $_;
	}
	close $in_fh;
} else {
	get_flist($in_path,$out_path,"reads_1st","fa");
}

##	multi-threading part for the 1st assembly;
my $thread_m = new Parallel::ForkManager($num_threads);
my $total = @files;
my $i;
foreach $file (@files){
	
	$i++;
	$thread_m->start and next;
	
	run_assembly($file, $in_path, "$out_path/assembly_1st", 1, $cmd);
	print STDERR "Assembling locus $i of $total \r";
	
	$thread_m->finish;
	
	
}

$thread_m->wait_all_children;
`find $in_path/ -name "*cap*" |xargs rm -rf`;

## 2nd assembly;
print "\nStarting the second run...\n";
#######

`echo -e "cap3 $cmd" > $out_path/log/assembly.2par`;
#######

if (! -e "$out_path/assembly_2nd") {
	
	`mkdir -p $out_path/assembly_2nd`;
	
}


get_flist("$out_path/assembly_1st/", $out_path, "reads_2nd", "fa");

foreach $file (@files){
	
	
	
	#print "start one threads" . "\n";

	$thread_m->start and next;
	
	run_assembly($file, "$out_path/assembly_1st", "$out_path/assembly_2nd", 2, $cmd);
	
	$thread_m->finish;
	
	
}

$thread_m->wait_all_children;
`find $out_path/assembly_1st/ -name "*cap*" |xargs rm -rf`;

##collect the final results, and delet useless files;
print "Starting to collect the final contigs...\n";

get_flist("$out_path/assembly_2nd", $out_path, "2nd_assembled","fa");	# get the list of second assembly files;

parse_fasta(\@files, "$out_path/assembly_2nd", "$out_path", "collected_final.fa");



##	if collect contigs for reads1;
if ($collect) {
	
	my $in_fa = $out_path . "/" . "assembly_1st";
	my $out_fa = $out_path;
	my $name_fa = "collected_reads1.fa";
	
	
	get_flist("$out_path/assembly_1st",$out_path,"collect_reads1","fa");
	parse_fasta(\@files, $in_fa, $out_fa, $name_fa);
	
	`sed -i /Consensus/d  $out_fa/$name_fa`;
	
}


sub parse_fasta {

	my ($files, $in_fa, $out_fa, $name_fa) =@_;
	my ($in_fh, $out_fh, $file);
	
	open ($out_fh, ">$out_fa/$name_fa") or die "$!";
	
	foreach $file (@$files) {
		open ($in_fh, "<$in_fa/$file") or die "$!";
		while (<$in_fh>) {
		
		chomp if /Consensus/;
		$_ =~ s/^>/>$file\./;
		print $out_fh $_;
		
		}
		
		close $in_fh;
	}
	
	#print $out_fh "\n";	#the last line;
	close $out_fh;
	
	
}

sub run_assembly {

	
	my ($file, $in_path, $out_path, $run, $cmd) = @_;
	
	`cap3 $in_path/$file $cmd 1>/dev/null 2>>$out_path/error.log`; #Changed, maybe some bugs here.
	unless (-z "$in_path/$file.cap.contigs") {
		
		`mv $in_path/$file.cap.contigs $out_path/$file`;
		
	}
	`rm -rf $in_path/$file.cap.*`;
	
	
	if ($run ==1) {
		my $locus = substr($file, 0, -3);
		my $seq = $seqs{$locus};
		open (my $out_fh, ">>$out_path/$file");
		print $out_fh 
		">" . $locus . "_Consensus" . "\n" . $seq . "\n";
		
		close $out_fh;
	}
	
	
	
	
	
}


sub get_flist {
	my ($in_path, $out_path, $name, $type) = @_;
	
	`mkdir -p $out_path/log/`;
	
	
	undef (@files);
	opendir (D, $in_path) or die "$!";
	open(my $out_list, ">$out_path/log/$name\.txt") or die "$!";

	##	get the file list;
	
	while (($file = readdir(D))) {
		next if $file !~ /.+\.$type$/;
		print $out_list $file . "\n";
		push (@files, $file);
	}
	
	return @files;
	close $out_list;
	close D;
}