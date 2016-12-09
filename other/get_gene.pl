#!/usr/bin/env perl
#Extract gene from blast results of japanese eel.
use strict;
use warnings;
my %pos;
my $usage = "get_pos.pl <quey results> <gene names>\n";
open(my $in_fh, "$ARGV[0]") or die "$!";

while (<$in_fh>) {
	next if /^$|^#/;
	chomp;
	my @parts   = split;
	my $q_name  = (split(/;/,$parts[0]))[0];
	my $q_pos   = substr((split(/;/,$parts[0]))[1],8);
	my $pos_ori = substr((split(/;/,$parts[0]))[2],8); #Original position.
	my $seq_len = substr((split(/;/,$parts[0]))[3],7);
	my $q_start = $parts[6];
	my $q_end   = $parts[7];
	my $len     = $q_pos - $q_start;
	my $s       = $parts[1]; #Name of scaffold.
	my $overlap = ($q_end - $q_start + 1)/$seq_len; #Mapped rate of sequence.
	#next if $overlap < 0.8; #skip short mapped query.
	my $s_start = $parts[8];
	my $s_end   = $parts[9];
	my $s_pos   = $s_end > $s_start ? ($s_start + $len) : ($s_start - $len); #approximate position for query seq in subject.
	#Test print.
	#print join("\t", $q_pos, $q_start, $len, $s_start, $s_end, $s_pos), "\n";
	$pos{$s}->{'new'}    = $s_pos; # Subject => Postion.
	$pos{$s}->{'old'}    = $q_name . '_' . $pos_ori;	

}
close $in_fh;
open($in_fh, "$ARGV[1]") or die "$!";

##Print header.
print join("\t", 'scaffold_name', 'start', 'end', 'ori_pos', 'gene_name'), "\n"; 
	
while (<$in_fh>) {

	next if /^$|^#/;
	chomp;
	my @parts = split;
	my $name  = (split(/_/, $parts[0]))[0];
	my $start = $parts[3]; #for each gene in the scaffold, extract from gff file.
	my $end   = $parts[4];
	next if not exists $pos{$name};
	my $pos_n = $pos{$name}->{'new'};
	my $dist  = abs($start - $pos_n) < abs($end - $pos_n) ? abs($start - $pos_n) : abs($end - $pos_n);
	#my @gene  = split(/\|/, $parts[2]);
	my $gene  = $parts[2];
	
	if ($dist <= 5000 or ($pos_n >= $start && $pos_n <= $end)) {
		
		print join("\t", $name, $start, $end, $pos{$name}->{'old'}, $gene), "\n";
		
	}
	
	
}
