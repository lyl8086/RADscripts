#!/bin/env perl
##use models from stacks to call snp from pileup files.
##Author: Yulong Li <liyulong12@mails.ucas.ac.cn>

use warnings;
use strict;
my ($in_fh, $in_bam, @indivs, $n1, $n2, $n3, $n4, $allele1, $allele2, $total);
my $bound_low  = 0.001;
my $bound_high = 0.1;
my $het_limit  = -3.84; #Set P = 0.05 
my $hom_limit  = 3.84;
my $stack      = 2; #Minimum reads for a locus to be processed.
my $infile     = $ARGV[0];
my $bam        = $ARGV[1];
my $usage      = 'Usage: Call_snp.pl <pileup files> <Pop files>';

die "$usage\n" unless @ARGV >= 2;

if ($infile eq '-c') {
	$in_fh = 'STDIN'; #Read STDIN.
} elsif (substr($infile, -2) eq 'gz') {
	open($in_fh, "gunzip -c $infile|") or die "$!";
} else {
	open($in_fh, "$infile") or die "$!";
}

open($in_bam, "$bam") or die "$!";

while (<$in_bam>) {
	next if /^$/;
	push @indivs, $_;
}
close($in_bam);

print join("\t", '#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT');
foreach (@indivs) {
	chomp;
	print "\t", $_;
}
print "\n";

while (<$in_fh>) {

	next if /^$/;
	my @parts  = split(/\t/, $_);
	my $chrom  = $parts[0];
	my $pos    = $parts[1];
	my $ref    = uc($parts[2]);
	my $alt    = '';
	my $filter = 0;
	my $miss   = 0;
	my @genos;
	undef @genos;
	for (my $i=3; $i<=$#parts; $i+=3) {
		
		my $dp    = $parts[$i]; #Depth of covered reads.
		my $base  = uc($parts[$i+1]);
		my $baseq = $parts[$i+2];
		my ($GT, $GL, $DP, $error_hom, $error_het, $ln_het, $ln_hom, $l_ratio, $geno);
		
		if ($dp < $stack or not $base) {
			##Do not have enough reads.
			$GT = './.';
			$GL = -255;
			$DP = 0;
			$miss += 1; #Recode missing;
			$geno = join (':', "$GT","$DP","$GL");
			push @genos, $geno;
			next;
			
		}
		
		basecount($base, $ref);
		
		
		if ($total == 0) {
			##Locus is Unknown.
			$GT = './.';
			$GL = -255;
			$DP = 0;
			$miss += 1; #Recode missing;
		} else {
			##SNP calling models, adapt from source code of stacks. Original Author: Julian Catchen <jcatchen@illinois.edu>.
			##Method of Paul Hohenlohe <hohenlohe@uidaho.edu>.
			
			##Sequence errors.
			$error_hom  = (4/3) * (($total - $n1) / $total);
			$error_het  = 2* (($n3 + $n4) / $total);
		
			##Bounded errors.
			if ($error_hom < $bound_low) {
				$error_hom = $bound_low;
			} elsif ($error_hom > $bound_high) {
				$error_hom = $bound_high;
			}
		
			if ($error_het < $bound_low) {
				$error_het = $bound_low;
			} elsif ($error_het > $bound_high) {
				$error_het = $bound_high;
			}
		
			##Calculate the log likelihood for the homozygous and heterozygous genotypes.
			
			#ln L(1/1) = ln(n! / n_1!n_2!n_3!n_4!) + 
			#    			n_1 * ln(n_1 / n) + 
            #   		   (n - n_1) * ln(n - n_1 / 3n)
			$ln_hom = $n1 * log(1 - ((3/4) * $error_hom));
			$ln_hom += $error_hom > 0 ? (($n2 + $n3 + $n4) * log($error_hom / 4)) : 0;
		

			#ln L(1/2) = ln(n! / n_1!n_2!n_3!n_4!) + 
			#              (n_1 + n_2) * ln(n_1 + n_2 / 2n) + 
			#              (n_3 + n_4) * ln(n_3 + n_4 / 2n)
			$ln_het = ($n1 + $n2) * log(0.5 - ($error_het / 4));
			$ln_het += $error_het > 0 ? (($n3 + $n4) * log($error_het / 4)) : 0;
		
		
			##Calculate the likelihood ratio.
			$l_ratio = 2 * ($ln_hom - $ln_het);
		
			if ($l_ratio <= $het_limit) {
				##Locus is Heterozygote.
				$GT = '0/1';
				$GL = sprintf "%.2f", $ln_het;
				$DP = $n1 + $n2; #Consider only the first two bases as the Heterozygote Depth.
				
				if ($allele1 eq $ref) {
				
					if ($alt) { 
						#Test if the alt alleles was consistent, we consider only biallelic loci.
						if ($alt ne $allele2) {
							$filter = 1;
							last;
						}
					}
					
					$alt = $allele2;
					
				} else {
				
					if ($alt) { 
						#Test if the alt alleles was consistent, we consider only biallelic loci.
						if ($alt ne $allele1) {
							$filter = 1;
							last;
						}
					}
					
					$alt = $allele1;
				}
			} elsif ($l_ratio >= $hom_limit) {
				##Locus is Homozygote.
				if ($allele1 eq $ref) {
					$GT = '0/0';
				} else {
					$GT = '1/1';
					
					if ($alt) { 
						#Test if the alt alleles was consistent, we consider only biallelic loci.
						if ($alt ne $allele1) {
							$filter = 1;
							last;
						}
						
					}
				}
				$GL = sprintf "%.2f", $ln_hom;
				$DP = $n1; #Consider only the first base as the Homozygote Depth.
			
			} else {
				##Locus is Unknown.
				$GT = './.';
				$GL = -255;
				$DP = 0;
				$miss += 1; #Recode missing;
			}
		}
		$geno = join (':', "$GT","$DP","$GL");
		push @genos, $geno;		
	}
	
	unless ($alt) {
		$alt = $ref; # Mark as monomorphic loci.
	}
	next if $filter == 1;
	next if $ref eq $alt; #Skip non variant site.
	next if $miss == @indivs; #Skip site with all missing.
	
	print join("\t", $chrom, $pos, '.', $ref, $alt, '.', '.', '.', 'GT:DP:GL', @genos), "\n";
	
}
close($in_fh);

sub basecount {
	my $base = shift;
	my $ref  = shift;
	
	if ($base =~ /[\-\+\>\<]|^$|\*/) {
		##Consider simple match and mismatch, skip indels and others, return unknown.
		$n1 = $n2 = $n3 = $n4 = 0;
		$allele1 = $allele2 = 'N';
	} else {
		$base =~ s/\.|\,/$ref/g; #Replace ', .' with reference base;
		my $bases->{'A'} = ($base =~ tr /Aa/Aa/);
		   $bases->{'C'} = ($base =~ tr /Cc/Cc/);
		   $bases->{'G'} = ($base =~ tr /Gg/Gg/);
		   $bases->{'T'} = ($base =~ tr /Tt/Tt/);
		my @numbers = sort {$b <=> $a} values %{$bases};#
		my @alleles = sort {$bases->{$b} <=> $bases->{$a}} keys %{$bases};#
		   $n1 = $numbers[0];
		   $n2 = $numbers[1];
		   $n3 = $numbers[2];
		   $n4 = $numbers[3];
		   $allele1 = $alleles[0];
		   $allele2 = $alleles[1];
	}
	$total = $n1 + $n2 + $n3 + $n4;
	return $n1, $n2, $n3, $n4, $allele1, $allele2, $total;
}

		
	
	
