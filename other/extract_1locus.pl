#!/usr/bin/env perl
use strict;
use warnings;

my ($fh, $fa, %seqs, $buf, $id, %loci, $ref);
die "Usage: perl extract_1locus.pl <fasta file>\n" unless @ARGV >= 1;
$fa = $ARGV[0];
$buf = "";
open($fh, "<$fa") or die "$!";
while(<$fh>) {
	next if /^$/;
	chomp;
	if (/^>/) {
		
		if (length($buf) > 0) {
			$seqs{$id}  = $buf;
			$buf = "";
		}
		$id = $_;
		my $locus = (split(/\./, $_, 2))[0];
		$locus =~ s/>//; #Locus name.
		#print $locus, "\n";
		$loci{$locus}++; #Use hash to record number;
		$ref->{$locus} = $id; #Just retain One Contig name for each loci;
	}
	else {
	    $buf .= $_;
	}
}
#The last one
$seqs{$id}  = $buf;

foreach (sort {$a<=>$b} keys %loci) {

	my $ctg_id = $ref->{$_};
	next if ($seqs{$ctg_id} !~ /GAATT$/ && defined $ARGV[1] && $ARGV[1] eq 'y');	
	print $ctg_id, "\n", $seqs{$ctg_id}, "\n" if $loci{$_} == 1;
		

}

		
