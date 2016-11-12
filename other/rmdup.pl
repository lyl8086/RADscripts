#!/usr/bin/env perl

##############               Readme                ###############
#	                                                         #
#                                                                #
#     Remove PCR duplications according the identical sequences. #
#     Usage: perl <fq file1> <fq file2> <if write results>       #
#                                                                #
#     Note: Leave 3rd option to blank if you just want to obtain #
#     the duplication rate.                                      #
#	                                                         #
#     Author: Yulong Li <liyulong12@mails.ucas.ac.cn>            #
##################################################################

use strict;
use warnings;

my $fa1 = $ARGV[0];
my $fa2 = $ARGV[1];
my $rw  = $ARGV[2] ? 'y' : '';
my (%seqs, %ids, %quals, $in_fh, $reads, $dup, $dup_rate);
die "Usage: $0 <fq file1> <fq file2> <write results>\n" unless @ARGV >= 2;

#open (my $seq_1, "gunzip -c $fa1 |") or die "$!";
#open (my $seq_2, "gunzip -c $fa2 |") or die "$!";
$in_fh = guessfmt($fa1);
open (my $seq_1, "$in_fh") or die "$!";
$in_fh = guessfmt($fa2);
open (my $seq_2, "$in_fh") or die "$!";

while (<$seq_1>) {

	my $id  = $_;
	
	my $seq  = <$seq_1>;	#seq
			   <$seq_1>;	#identifier
	my $qual = <$seq_1>;	#quality
					
			   <$seq_2>;	#id
	  $seq  .= <$seq_2>;	#Combined seq
			   <$seq_2>;	#identifier
	  $qual .= <$seq_2>;	#Combined quality
	  
	  $seqs{$seq}++;
	  $ids{$seq}   = $id; #Retain only one id for the same paired reads.
	  $quals{$seq} = $qual;#same as above.
	  $reads++;
}
close($seq_1);
close($seq_2);

$fa1 =~ s/\.gz//;
$fa2 =~ s/\.gz//;
open (my $out_1, "|bgzip -c >$fa1.rmclone.gz") or die "$!" if ($rw eq 'y');
open (my $out_2, "|bgzip -c >$fa2.rmclone.gz") or die "$!" if ($rw eq 'y');

foreach (keys %seqs) {

	$dup += ($seqs{$_} - 1) if $seqs{$_} > 1;
	if ($rw eq 'y') {
	
		my ($id_same, $id_diff) = split(/\ /, $ids{$_}, 2);
		my $id_1 = $id_same . ' ' . $id_diff;
		$id_diff =~ tr /12/21/; #ID for paired read.
		my $id_2 = $id_same . ' ' . $id_diff;
		my ($fq_1, $fq_2) = split(/\n/, $_, 2);
		my ($qual_1, $qual_2) = split(/\n/, $quals{$_}, 2);
		#Reads1:
		print $out_1 $id_1, $fq_1, "\n", "+\n", $qual_1, "\n";
		#Reads2:
		print $out_2 $id_2, $fq_2, "+\n", $qual_2;
	}
	
}

$dup_rate = sprintf "%.2f", $dup*100/$reads;

print "Duplication rate is $dup_rate%, Total paired reads are $reads.\n";

close($out_1) if ($rw eq 'y');
close($out_2) if ($rw eq 'y');

sub guessfmt {
	
	my $infile = shift;
	my $fmt;
	if ($infile eq '-c') {
		$fmt = 'STDIN'; #Read STDIN.
	} elsif (substr($infile, -2) eq 'gz') {
		$fmt = "gunzip -c $infile|" #Gzipped file.
	} else {
		$fmt = "$infile"
	}
	return $fmt;
}


	
