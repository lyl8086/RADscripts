#!/bin/env perl
#Extract seq flanking specified length around SNP.
##Author: Yulong Li

use strict;
use warnings;

my $in_genome = $ARGV[0];
my $in_snps   = $ARGV[1];
my $span      = $ARGV[2] ? $ARGV[2] : 5000;
my ($snps, $genome, %seqs, $id);
my $buf = "";
die "Usage: perl cut_genome_pos.pl <fasta file> <SNP positions> <length>\n" unless @ARGV >= 2;

open ($genome, "$in_genome") or die "$!";
while (<$genome>) {
	
	next if /^$/;
	chomp;
	if (/^>/) {
		if (length($buf) > 0) {
			$seqs{$id}  = $buf;
			$buf = "";
		}
		$id = (split(/ /, $_))[0];
	}
	else {
	    $buf .= $_;
	}
}
#The last one
$seqs{$id}  = $buf;

open ($snps, "$in_snps") or die "$!";

while (<$snps>) {

	next if /^$|^#/;
	chomp;
	my @parts   = split;
	my $contig  = ">" . $parts[0];
	next if not exists $seqs{$contig};
	my $pos     = $parts[1];
	my ($pos_low, $pos_high, $pos_now, $len, $cut, @seq, $seq, $sub_seq);
	$len        = length($seqs{$contig});
	#@seq        = split(//,$seqs{$contig});
	$seq        = $seqs{$contig};
	
	
=cut
	if (($pos - 2500) > 0) {
		$pos_low = ($pos - 2500);
		if (($pos + 2500) < $len) {
			$pos_high = ($pos + 2500);
		} else {
			$pos_high = $len;
			#We need adjust lower position.
			$cut = $pos + 2500 - $len;
			$pos_low = ($pos - 2500 - $cut) > 0 ? ($pos - 2500 - $cut) : 0;
		}
	} else {
		
		$pos_low = 0;
		#We need adjust higher position.
		$cut = 2500 - $pos;
		$pos_high = ($pos + 2500 + $cut) < $len ? ($pos + 2500 + $cut) : $len;
		
	}
=cut
	$pos_low  = ($pos - $span) > 0 ? ($pos - $span) : 0;
	$pos_high = ($pos + $span) < $len ? ($pos + $span) : ($len - 1);
	$pos_now  = ($pos - $span) > 0 ? $span : $pos;
	
=cut
	##When the high bound is short, adjust the lower bound.
	if (($pos - 2500) > 0 && ($pos_low + 5000) > $len) {
		$pos_low = ($pos_high - 5000) > 0 ? ($pos_high - 5000) : 0;
	}
	
	my $pattern = 'N'x200
	while (index($seq, $pattern, $pos_low)) {
		$pos_low += 1;
	}

=cut
	my $seq_len           = $pos_high - $pos_low; 
	$sub_seq->{'seq'}     = substr($seq, $pos_low, $seq_len);
	$sub_seq->{'id'}      = $contig;
	$sub_seq->{'pos'}     = $pos;
	$sub_seq->{'pos_now'} = $pos_now;
	$sub_seq->{'len'}     = length($sub_seq->{'seq'});
	## Output results.
	print_fasta($sub_seq);
	

	
}

sub print_fasta {
    my ($seq) = @_;

    my $tmp = $seq->{'seq'};
	
    print $seq->{'id'}, ';Pos_new:', $seq->{'pos_now'}, ';Pos_ori:', $seq->{'pos'}, ';Length:', $seq->{'len'}, "\n";

    while (length($tmp) > 125) {
        print substr($tmp, 0, 125), "\n";
        $tmp = substr($tmp, 125);
    }

    print $tmp, "\n" if (length($tmp) > 0);
}


	
			
	
	
	
