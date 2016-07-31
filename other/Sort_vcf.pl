#!/usr/bin/env perl


##############		           Readme 			   ##############
#	                                                                        #
#                                                                               #
#	  sort the order of individuals in vcf according to the popmap file     #
#	  Usage: perl Sort_vcf.pl vcf popmap > outfile                          #
#         popmap file: indv1	pop1                                            #
#		       indv2	pop1                                            #
#		       indv3	pop2                                            #
#	                                                                        #
#         Author: Yulong Li <liyulong12@mails.ucas.ac.cn>                       #
#################################################################################
use strict;
use warnings;
die "Usage: perl Sort_vcf.pl [vcf file] [popmap file] > outfile\n" unless @ARGV == 2;
my (@header, @vcf, %samples, %pops, @id, $in_fh);
my $cnt = 0;
my $infile = $ARGV[0];

if (substr($infile, -2) eq 'gz') {
	open($in_fh, "gunzip -c $infile|") or die "$!";
	} else {
	open($in_fh, "$infile") or die "$!";
}


while (<$in_fh>) {
	next if /^##|^$/;
	chomp;
	if (/^#/) {
		push @header, split;
	}
	$cnt++;
	push @vcf, [ split ];
}

close $in_fh;

#Create index for samples;
for (my $i=9; $i<=$#header; $i++) {
	$samples{$header[$i]} = $i;
}

#Read population files;
open(my $pop, "$ARGV[1]") or die "$!";

while (<$pop>) {
	next if /^##|^$/;
	chomp;
	my @part = split;
	
	push @{$pops{$part[1]}}, $samples{$part[0]};
}
close $pop;

#Sort population index;
foreach my $key (sort keys %pops) {
	push @id, @{$pops{$key}};
}

#Print results;	
for (my $i=0; $i<$cnt; $i++) {
	
		print join("\t", @{$vcf[$i]}[0..8,@id]), "\n";
	
}
