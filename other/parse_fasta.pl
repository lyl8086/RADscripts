#!/usr/bin/env perl
#filter length, extract seqs given names of fasta...
#Author: Yulong Li
use strict; use warnings;
die "Usage: parse_fasta.pl <file> <min length>\n" unless @ARGV >= 2;
#Optional usage: parse_fasta.pl <file> <min length> [y for retain EcoRI].
open(my $in, "<$ARGV[0]") or die "error openning $ARGV[0]";
my $min_len = $ARGV[1];
my ($buf, $id, $ctg_id, %seqs, @ctgs);
$buf = "";
while (<$in>) {
	chomp;
	if (/^>/) {
		if (length($buf) > 0) {
			$seqs{$id}  = $buf;
			$buf = "";
		}
		$id = $_;
	}
	else {
	    $buf .= $_;
	}
}
#The last one
$seqs{$id}  = $buf;
#Output results
if (defined $ARGV[2] && $ARGV[2] ne 'y') {
	open($ctg_id, "<$ARGV[2]");
	while (<$ctg_id>) {
		chomp;
		my $ctg = ">".$_;
		push(@ctgs, $ctg);
	}
	close $ctg_id;
	foreach $id (@ctgs) {
	if ( exists $seqs{$id}) {
		print $id, "\n", $seqs{$id}, "\n";	
		}
	}
}

	
else {
	foreach $id (keys %seqs) {
	next if (length($seqs{$id}) < $min_len);
	next if ($seqs{$id} !~ /GAATT$/ && defined $ARGV[2] && $ARGV[2] eq 'y');
	print $id, "\n", $seqs{$id}, "\n";
	}
}

close $in