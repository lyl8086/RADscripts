#!/usr/bin/env perl
## Fst estimator for VCF using Nei, Weir and Cockerham, Arlequin, Stacks and Reich formulas.
## Author: Yulong Li <liyulong12@mails.ucas.ac.cn>

use strict;
use warnings;
use List::MoreUtils qw{uniq};
use Getopt::Long qw(:config no_ignore_case bundling);

my ($in_fh, $out_fh, @header, $cnt, $pops, @order, %samples, $sum_fst);
my ($infile, $outfile, $pop, $indivs, $popmap, $sum_a, $sum_t, $cnt_fst);
$infile = $ARGV[0];
$popmap = $ARGV[1];
die "Usage: Fst_estimator.pl <vcf> <popmap>\n" unless @ARGV >= 2;

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
		for (my $i=9; $i<=$#header; $i++) {
			$samples{$header[$i]} = $i;
		}
		
		##################
		
		# Read PopMap file;
		
		# Popmap file:
		# indiv1 pop1
		# indiv2 pop2
		# indiv3 pop1 ...
		
		##################
		
		open(my $pop, "$popmap") or die "No PopMap file!";

		while (<$pop>) {
			next if /^##|^$/;
			chomp;
			my @part = split;
			push @order, $part[1];
			push @{$indivs->{$part[1]}}, $samples{$part[0]}; # Pop name => @indiv rank. 
			
		}
		close $pop;
		die "Header is wrong!" if @order != (@header -9);
		@order = uniq @order;
		#print @order, "\n";
		next;
		
	}
	
	#################################################################
	
	# VCF format:
	# CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Indiv_genotypes ...
	# Genotypes:GT:PL:DP:SP:GQ ...
	
	#################################################################
	
	my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genos) = split;
	my ($fmt, @gts_total, $AA, $Aa, $aa);
	undef @gts_total;
	$cnt++; # SNP number.
	
	my @formats = split(/:/, $format);                                # order of format:
	for (my $i=0; $i<=$#formats; ++$i) { $fmt->{$formats[$i]} = $i;}; # Geno => Order.
	
	# Iteration each population.
	foreach my $pop_name (@order) {
	
		$AA->{$pop_name} = 0;
		$Aa->{$pop_name} = 0;
		$aa->{$pop_name} = 0;
		foreach my $rank(@{$indivs->{$pop_name}}) {
			# each individual of every population.
			my @geno = split(/:/, $genos[$rank-9]);
			my $GT   = $geno[$fmt->{'GT'}];
			$AA->{$pop_name}++ if $GT eq '0/0';
			$Aa->{$pop_name}++ if $GT eq '0/1';
			$aa->{$pop_name}++ if $GT eq '1/1';
		}
		
	}
	
	# Pairwise populations.
	for (my $i=0; $i<$#order; ++$i) {
	
		my $pop1 = $order[$i];
		
		# n=>number of alelles; p=>number of A alleles; het=>number of heteozygotes.
		
		my $n1  = $AA->{$pop1} + $Aa->{$pop1} + $aa->{$pop1};
		   $n1 *= 2;
		my $p1   = 2 * $AA->{$pop1} + $Aa->{$pop1};
		my $het1 = $Aa->{$pop1};
		
		for (my $j=$i+1; $j<=$#order; ++$j) {
			
			
			my $pop2 = $order[$j];
			my $n2   = $AA->{$pop2} + $Aa->{$pop2} + $aa->{$pop2};
			   $n2  *= 2;
			my $p2   = 2 * $AA->{$pop2} + $Aa->{$pop2};
			my $het2 = $Aa->{$pop2};
			my ($Fst, @Fsts, $a, $sum);
			($Fst->{'wc'},  $a->{'wc'},  $sum->{'wc'})  = weir_and_cockerham_fst($n1, $n2, $p1, $p2, $het1, $het2);
			($Fst->{'ar'},  $a->{'ar'},  $sum->{'ar'})  = arlequin_fst($n1, $n2, $p1, $p2); # Arlequin.
			($Fst->{'frq'}, $a->{'frq'}, $sum->{'frq'}) = weir_and_cockerham_freq($n1, $n2, $p1, $p2); #
			$Fst->{'st'}    = stacks_fst($n1, $n2, $p1, $p2); # Stacks.
			$Fst->{'nei'}   = nei_fst($n1, $n2, $p1, $p2);
			($Fst->{'rei'}, $a->{'rei'}, $sum->{'rei'}) = Reich_fst($n1, $n2, $p1, $p2); #Reich's Fst.
			my @para = ('wc', 'ar', 'frq', 'st', 'nei', 'rei');
			foreach my $i (@para) {
				$cnt_fst->{$pop1.$pop2}->{$i}++ if $Fst->{$i} ne 'NAN'; # Count of Fst.
				$sum_fst->{$pop1.$pop2}->{$i} += $Fst->{$i} if $Fst->{$i} ne 'NAN'; # Sum of Fst.
				
				$sum_a->{$pop1.$pop2}->{$i} += $a->{$i} if defined $a->{$i};
				$sum_t->{$pop1.$pop2}->{$i} += $sum->{$i} if defined $sum->{$i};
				push @Fsts, sprintf("%9.6f", $Fst->{$i});
				
			}
			
			my $out_name = $pop1 . "-" . $pop2 . '.fst';
			if ($cnt == 1) {
				open($out_fh, ">$out_name") or die "$!";
				printf $out_fh "%-9s\t"x8,('CHROM', 'POS', 'wcFst', 'Arlequin', 'Freq', 'Stacks', 'Nei', 'Reich');
				print $out_fh "\n";
				close $out_fh;
			}
			open($out_fh, ">>$out_name") or die "$!";
			print $out_fh join("\t", $chrom, $pos, @Fsts), "\n";
		}
	}
	
	
}
close $out_fh;
# Output matrix.
print "Weir and Cockerham Weighted Fst:\n";
print_matrix('weight', 'wc');
# 
print "\nArlequin Weight Fst:\n";
print_matrix('weight', 'ar');
#
print "\nFrequency Weight Fst:\n";
print_matrix('weight', 'frq');
#
print "\nStacks Mean Fst:\n";
print_matrix('mean', 'st');
#
print "\nNei Mean Fst:\n";
print_matrix('mean', 'nei');
#
print "\nReich wieght Fst:\n";
print_matrix('weight', 'rei');	
######

sub print_matrix {
	
	my $flag_1 = shift;
	my $flag_2 = shift;
    printf "%10s", " ";
	map {printf "%10s", $_} @order; # First line.
	print "\n";
	for (my $i=0; $i<=$#order; ++$i) {
	
		my $pop1 = $order[$i];
		for (my $j=0; $j<=$i; ++$j) {
			my $pop2 = $order[$j];
			printf "%10s", $pop1 if $j == 0; 
			if ($j == $i) {printf "%10.5f\n",'0.00000';next;}
			if ($flag_1 eq 'weight') {
				my $weit  = $sum_a->{$pop2.$pop1}->{$flag_2};
				   $weit /= $sum_t->{$pop2.$pop1}->{$flag_2};
				printf "%10.5f", $weit;
			} elsif ($flag_1 eq 'mean') {
				my $mean  = $sum_fst->{$pop2.$pop1}->{$flag_2};
				   $mean /= $cnt_fst->{$pop2.$pop1}->{$flag_2};
				printf "%10.5f", $mean;
			}
		}
	}	
}
	
sub weir_and_cockerham_fst {
	
	my ($n1, $n2, $p1, $p2, $het1, $het2) = @_;
	my $n_bar = ($n1 + $n2)/4;
	my $nc    = 2*$n_bar - (($n1/2)**2 + ($n2/2)**2)/(2*$n_bar);
	my $p_bar = ($p1 + $p2)/($n1 + $n2);
	my $h_bar = 2*($het1 + $het2)/($n1 + $n2);
	my $ssa   = $n1/2 * ($p1/$n1 - $p_bar)**2;
	   $ssa  += $n2/2 * ($p2/$n2 - $p_bar)**2;
	   $ssa  /= $n_bar;
	my $a     = $p_bar*(1 - $p_bar) - 0.5*$ssa - 0.25*$h_bar;
	   $a    *= 1/($n_bar - 1);
	   $a     = ($ssa - $a) * ($n_bar/$nc);
	my $b     = (2*$n_bar-1)/$n_bar;
	   $b     = $p_bar*(1 - $p_bar) - 0.5*$ssa - $b*0.25*$h_bar;
	   $b    *= $n_bar/($n_bar - 1);
	my $c     = 0.5*$h_bar;
	my $sum   = $a + $b + $c;
	#
	if ($sum == 0) {
		my $fst = 'NAN'; 
		return ($fst, $a, $sum);
		next;
	}
	#
	my $fst   = $a/$sum;
	
	return ($fst, $a, $sum);
}

sub weir_and_cockerham_freq {

	my ($n1, $n2, $p1, $p2) = @_;
	my $n_bar = ($n1 + $n2)/4;
	my $nc    = 2*$n_bar - (($n1/2)**2 + ($n2/2)**2)/(2*$n_bar);
	my $p_bar = ($p1 + $p2)/($n1 + $n2);
	my $ssa   = $n1/2 * ($p1/$n1 - $p_bar)**2;
	   $ssa  += $n2/2 * ($p2/$n2 - $p_bar)**2;
	   $ssa  /= $n_bar;
	my $a     = ($p_bar*(1-$p_bar) - 0.5*$ssa)/(2*$n_bar - 1);
	   $a     = $ssa - $a;
	my $sum   = (2*$nc - 1)*$p_bar*(1 - $p_bar)/(2*$n_bar - 1);
	   $sum  += (1 + 2*($n_bar - $nc)/(2*$n_bar - 1))*0.5*$ssa;
	if ($sum == 0) {return('NAN', $a, $sum);next;}
	my $fst   = $a / $sum;
	return ($fst, $a, $sum);
}
	

sub arlequin_fst {

	my ($n1, $n2, $p1, $p2) = @_;
	my $b       = ($n1 + $n2);
	my $n       = $b; # It is 2N indeed. 
       $b      -= ($n1**2 + $n2**2)/$n;
	my $ssd_T   = ($p1 + $p2) * ($n1 - $p1 + $n2 - $p2)/$n;
	my $ssd_WP  = $p1*($n1 - $p1)/$n1;
	   $ssd_WP += $p2*($n2 - $p2)/$n2;
	my $wp      = ($n-2)*$ssd_T - ($n-1)*$ssd_WP;
	my $total   = ($n-2)*$ssd_T - ($n-1-$b)*$ssd_WP;
	if ($total == 0) {return('NAN', $wp, $total); next;}
	my $fst     = $wp/$total;
	
	return ($fst, $wp, $total);

}

sub stacks_fst {

	my ($n1, $n2, $p1, $p2) = @_;
	my $q1     = $n1 - $p1;
	my $q2     = $n2 - $p2;
	my $p      = $p1 + $p2;
	my $q      = $n1 + $n2 - $p;
	my $pi_1   = calc_pi($p1, $q1, $n1);
	my $pi_2   = calc_pi($p2, $q2, $n2);
	my $pi_all = calc_pi($p, $q, ($n1+$n2));
	my $coef_1 = binom_coeff($n1, 2);
	my $coef_2 = binom_coeff($n2, 2);
	my $coef   = $coef_1 + $coef_2;
	if ($pi_all == 0 or $coef == 0) {
		return 'NAN';
		next;
	}
	
	my $fst    = $pi_1*$coef_1 + $pi_2*$coef_2;
	   $fst   /= $pi_all*$coef;
	   $fst    = 1 - $fst;
	return $fst;
	
}

sub calc_pi {
	
	# p=>number of A alleles; n=>number of total alleles.
	my ($p, $q, $n) = @_;
	my $pi  = binom_coeff($p, 2) + binom_coeff($q, 2);
	   $pi /= binom_coeff($n, 2);
	   $pi  = 1 - $pi;
	return $pi;
}

sub binom_coeff {
	
	my ($n, $k) = @_;
	if ($n < $k) {return 0; next;}
	
	# From Stacks:
    # Compute the binomial coefficient using the method of:
    # Y. Manolopoulos, "Binomial coefficient computation: recursion or iteration?",
    # ACM SIGCSE Bulletin, 34(4):65-67, 2002.
    
    my $r = 1;
    my $s = $k < ($n - $k) ? ($n - $k + 1) : ($k + 1);

    for (my $i = $n; $i >= $s; $i--){
		$r = $r * $i / ($n - $i + 1);
	}
    return $r;
}

sub nei_fst {

	my ($n1, $n2, $p1, $p2) = @_;
	my $p    = ($p1 + $p2)/($n1 +$n2);
	   $p1  /= $n1;
	   $p2  /= $n2;
	if ($p == 0 || $p == 1) {return 'NAN'; next;}
	my $fst  = $n1 * $p1 * (1 - $p1);
	   $fst += $n2 * $p2 * (1 - $p2);
	   $fst /= ($n1 + $n2) * $p * (1 - $p);
	   $fst  = 1 - $fst;
	return $fst;
	
}

sub Reich_fst {
	
	#
	#David Reich et al. (2009) Reconstructing Indian population history. 
	#Nature 461:489-494, especially Supplement 2.
	#
	#n=> num of alleles for a pop; p=> num of A allele.
	my ($n1, $n2, $p1, $p2) = @_;
	my ($N, $D, $h1, $h2);
	$h1 = $p1*($n1 - $p1)/($n1*($n1-1));
	$h2 = $p2*($n2 - $p2)/($n2*($n2-1));
	$N  = ($p1/$n1 - $p2/$n2)**2;
	$N -= ($h1/$n1 + $h2/$n2);
	$D  = $N + $h1 + $h2;
	
	if ($D == 0) {return ('NAN', 0, 0); next;}
		
	my $fst = $N/$D;
	
	return($fst, $N, $D);
	
}
