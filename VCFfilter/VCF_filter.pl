#!/usr/bin/env perl
#Author: Yulong Li <liyulong12@mails.ucas.ac.cn>
use strict;
use warnings;
use List::MoreUtils qw{uniq};
use Getopt::Long qw(:config no_ignore_case bundling);
my (@header, $cnt, $pops, @order, $in_fh, $out_fh, %samples);
my ($cmd, $infile, $outfile, $popmap, $minDP, $maxDP, $het, $fis);
my ($minGQ, $minQ, $l_maf, $g_maf, $global, $num_threads, $help);
my $cov = 0;
GetOptions (
			"in=s" 		=> \$infile,
			"out=s" 	=> \$outfile,
			"Popmap=s"      => \$popmap,
			"minDP=s" 	=> \$minDP,
			"MaxDP=s" 	=> \$maxDP,
			"Het=f"         => \$het,
			"Fis=f"         => \$fis,
			"GQ=i"          => \$minGQ,
			"Q=i"           => \$minQ,
			"localMAF=f"    => \$l_maf,
			"globalMAF=f"   => \$g_maf,
			"filter:1"      => \$global,
			"threads=i" 	=> \$num_threads,
			"coverage=f" 	=> \$cov,
			"help:1"	=> \$help
			)
			or die ("Error in command line arguments\n");

my $usage = 
'	
	Options [default values]:
	
	--i input vcf file.
	--o output file.
	--P PopMap file.
	--c coverage for each population.
	--m min depth of each individual.
	--M max depth of each individual.
	--H max Ho of each population.
	--F abs Fis value.
	--G min Genotype quality.
	--Q min quality for each SNP.
	--g Global MAF.
	--l Local MAF.
	--f apply the filter (Ho and Fis) on each site instead of each population.
	--t number of threads.
	--h help.
';
	
die "$usage\n" if $help or !($infile && $outfile && $popmap);

if (substr($infile, -2) eq 'gz') {
	open($in_fh, "gunzip -c $infile|") or die "$!";
	} else {
	open($in_fh, "$infile") or die "$!";
}

if (substr($outfile, -2) eq 'gz') {
	open($out_fh, "|bgzip -c >$outfile") or die "$!";
	} else {
	open($out_fh, ">$outfile") or die "$!";
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
		
		#Read PopMap file;
		
		#Popmap file:
		#indiv1	pop1
		#indiv2 pop2
		#indiv3 pop1 ...
		
		##################
		
		open(my $pop, "$popmap") or die "No PopMap file!";

		while (<$pop>) {
			next if /^##|^$/;
			chomp;
			my @part = split;
			push @order, $part[1];
			push @{$pops->{$part[1]}}, $samples{$part[0]}; #Pop name => @indiv rank. 
			
		}
		close $pop;
		die "Header is wrong!" if @order != (@header -9);
		@order = uniq @order;
		
		###### Print header ######
		print $out_fh join("\t", @header[0..8]);
		foreach $pop (@order) {print $out_fh "\t", join("\t", @header[@{$pops->{$pop}}])};
		print $out_fh "\n";
		next;
		
	}
	
	#################################################################
	
	#VCF format:
	#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Indiv_genotypes ...
	#Genotypes:GT:PL:DP:SP:GQ ...
	
	#################################################################
	
	my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genos) = split;
	my ($fmt, @filtered, @gts_total);
	undef @filtered;
	undef @gts_total;
	my $cnt_Fis          = 0;
	my $cnt_Ho           = 0;
	my $cnt_gmaf         = 0;
	my $cnt_lmaf->{$ref} = 0;
	   $cnt_lmaf->{$alt} = 0;
	next if $alt =~ /,|\./;                                           #Skip non-biallelic loci.
	next if (defined $minQ && $qual < $minQ);                         #minmum quality.
	my $recode = 0;                                                   #Populations number of non-enough coverage. 
	my @formats = split(/:/, $format);                                #order of format:
	for (my $i=0; $i<=$#formats; ++$i) { $fmt->{$formats[$i]} = $i;}; #Geno => Order.
	
	#Iteration each population.
	foreach my $name (@order) {

		my $total     = @{$pops->{$name}};
		my $miss      = 0;
		my $cov_ratio = 1;
		my @gts;
		undef @gts;
		
		
		foreach my $rank(@{$pops->{$name}}) {
			#each individual of every population.
			my @geno = split(/:/, $genos[$rank-9]);
			my $GT   = $geno[$fmt->{'GT'}];
			my $DP   = $geno[$fmt->{'DP'}];
			my $GQ   = $geno[$fmt->{'GQ'}] if $fmt->{'GQ'};
			
			#Other filter...
			###### Depth ######
			if ((defined $minDP && $DP < $minDP) or (defined $maxDP && $DP > $maxDP)) {
				
				$geno[$fmt->{'GT'}] = './.';
				$miss++;
				push @filtered, join(":", @geno);
				next;
			}
			###### GQ ######
			if (defined $minGQ && $GQ < $minGQ) {
				$geno[$fmt->{'GT'}] = './.';
				$miss++;
				push @filtered, join(":", @geno);
				next;
			}
			
			$miss++ if $GT eq './.';
			push @filtered, join(":", @geno);
			push (@gts, $GT) if $GT ne './.';                            #For calculate statistics.
			push (@gts_total, $GT) if $GT ne './.' && ($global or $g_maf); #Global genotypes.
		}
		###### Coverage ######
		$cov_ratio = (1-$miss/$total);
		$recode++ if $cov_ratio < $cov;
		last if $cov_ratio < $cov;
		
		###### Local ######
		if ((!$global && ($het or $fis)) or $l_maf) {
		
			my ($Ho, $He, $Fis, $maf, $flag) = calc_stat(\@gts, $ref, $alt); #Each population.
		
			###### Het ######
			if (defined $het && $Ho > $het) {
				$cnt_Ho++;
				last;
			}
		
			###### Fis ######
			if (defined $fis && abs($Fis) > $fis) {
				$cnt_Fis++;
				last;
			}
		
			###### MAF ######
			$cnt_lmaf->{$flag}++ if defined $l_maf && $maf < $l_maf;
		}
		
	}
	
	next if $recode  > 0;
	next if $cnt_Ho  > 0;
	next if $cnt_Fis > 0;
	
	###### Global ######
	if ($global or $g_maf) {
		
		my ($Ho, $He, $Fis, $maf, $flag) = calc_stat(\@gts_total, $ref, $alt);
		
		###### Het ######
		if ($global && $het && $Ho > $het) {
			$cnt_Ho++;
			last;
		}
		
		###### Fis ######
		if ($global && $fis && abs($Fis) > $fis) {
			$cnt_Fis++;
			last;
		}
		
		###### MAF ######
		$cnt_gmaf++ if defined $g_maf && $maf < $g_maf;
		if (not defined $l_maf) {
			$cnt_lmaf->{$ref} = @order;
			$cnt_lmaf->{$alt} = @order;
		}
		
	}
	
	next if $cnt_Ho   > 0;
	next if $cnt_Fis  > 0;
	next if $cnt_gmaf > 0 && ($cnt_lmaf->{$ref} == @order or $cnt_lmaf->{$alt} == @order);
	
	###### Print each snp sites ######
	print $out_fh join("\t", $chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format), "\t";
	print $out_fh join("\t", @filtered), "\n";
	
}

close $out_fh;

sub calc_stat {

	my ($gts, $ref, $alt) = @_;
	my $allele->{$ref} = 0;
	   $allele->{$alt} = 0;
	   $allele->{'het'}= 0;
	
	foreach my $gt (@$gts) {

		$allele->{$ref} += 2 if $gt eq '0/0';
		$allele->{$alt} += 2 if $gt eq '1/1';
		if ($gt eq '0/1') {
			$allele->{$ref}++; 
			$allele->{$alt}++;
			$allele->{'het'}++;
		}
	}
	my $n    = ($allele->{$ref} + $allele->{$alt}) / 2;
	my $p    = $allele->{$ref} / (2 * $n);
	my $Ho   = $allele->{'het'} / $n;
	my $He   = 2 * $p * (1 - $p);
	my $Fis  = $He > 0 ? (1 - $Ho/$He) : 'NAN';
	my $flag = $p < 0.5 ? $ref : $alt;
	my $maf  = $p < 0.5 ? $p : (1 - $p);
	
	return($Ho, $He, $Fis, $maf, $flag);
	
}
