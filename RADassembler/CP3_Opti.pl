#!/usr/bin/env perl
##	Author:Yulong Li <liyulong12@mails.ucas.ac.cn>.
use strict; 
use warnings; 
use Parallel::ForkManager;
use Getopt::Long;


my (@files, $file, $out_fh, $in_fh);
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
	-h help
';

GetOptions ("par=s" 		=> \$cmd,
			"in=s" 			=> \$in_path,
			"out=s" 		=> \$out_path,
			"loci=s" 		=> \$loci,
			"threads=i" 	=> \$num_threads,
			"collect:1" 	=> \$collect,
			"help:1"		=> \$help
			)
			or die ("Error in command line arguments\n");
	
die "$usage\n" if !($in_path && $out_path && $loci) or ($help);

## If the outpath does not exist, create it;

if (! -e "$out_path/assembly_1st") {
	
	`mkdir -p $out_path/assembly_1st`;
	
}

##delete any tmp files;

`find $in_path/ -name "*cap*" |xargs rm -rf`;

get_flist($in_path,$out_path,"reads_1st","fa");


print "Starting the first run...\n";
`echo cap3 $cmd > $out_path/log/assembly.1par`;


##	extract consensus from stacks catalogs;
open (my $in_tags, "gunzip -c $loci |") or die "$!";
my %seqs;
while (<$in_tags>) {
	
	next if $_ =~ /^#/;
	my $seq = (split(/\t/, $_))[9];
	my $id = (split(/\t/, $_))[2];
	die "catalog file is wrong!!\n" if ($id+1 ne $.); 
	$seq = uc(reverse $seq);
	$seq =~ tr/ATCGN/TAGCN/;
	$seqs{$id} = $seq;
		
}

close $in_tags;

##	multi-threading part for the 1st assembly;
my $thread_m = new Parallel::ForkManager($num_threads);

foreach $file (@files){
	
	$thread_m->start and next;
	run_assembly($file, $in_path, "$out_path/assembly_1st", 1, $cmd);
	$thread_m->finish;
	
	
}

$thread_m->wait_all_children;
`find $in_path/ -name "*cap*" |xargs rm -rf`;

## 2nd assembly;
print "Starting the second run...\n";
#######
$cmd = '-r 0 -k 0'; #Parameters for the second assembly.
`echo cap3 $cmd > $out_path/log/assembly.2par`;
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
	
	foreach $file (@files) {
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
	
	system("cap3 $in_path/$file $cmd 1>/dev/null 2>>$out_path/error.log && cp $in_path/$file.cap.contigs $out_path/$file"); #Changed, maybe some bugs here.
		
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