#!/usr/bin/env perl

##############		           Readme 			        ##############
#	                                                                 #
#                                                                    #
#	  Convert vcf to DIYABC snp file                                 #
#     assuming all the loci to be autosome with equal sex ratio      #
#	  Usage: perl Vcf2diyabc.pl vcf > outfile                        #
#                                                                    #
#	                                                                 #
#     Author: Yulong Li <liyulong12@mails.ucas.ac.cn>                #
######################################################################
use strict;
use warnings;
die "Usage: perl Vcf2diyabc.pl [vcf file] > outfile\n" unless @ARGV == 1;
my (@header, @parts, $cnt, $in_fh);
my $vcf = $ARGV[0];

if (substr($vcf, -2) eq 'gz') {
	open($in_fh, "gunzip -c $vcf|") or die "$!";
	} else {
	open($in_fh, "$vcf") or die "$!";
}

while (<$in_fh>) {
	
		next if /^##|^$/;
		if (/^#/) {
			@header = split;
			next;
		}
		$cnt++;
		push @parts, [ split ];
}

close $in_fh;

print 'RAD SNPs file for DIYABC <NM=1NF> <MAF=hudson>', "\n", join("\t", 'IND', 'SEX', 'POP', split(//, 'A'x$cnt)),"\n";

for (my $i=9; $i<@header; $i++) {
	
	##print data...
	print join("\t", $header[$i], '9', substr($header[$i],0, 3));
	for (my $j=0; $j<$cnt; $j++) {
		
		my $allele;
		my $GT = (split(/:/, $parts[$j][$i]))[0];
		if ( $GT eq '0/0') {
			$allele = 2
		} elsif ( $GT eq '0/1') {
			$allele = 1
		} elsif ( $GT eq '1/1') {
			$allele = 0
		} else {
			$allele = 9
		}
		
		print "\t", $allele;
	}
	print "\n";
}

print STDERR "Done!! Total number of SNPs is $cnt\n";
 

			