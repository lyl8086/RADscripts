#!/usr/bin/env perl
#@for RADassembler.
#@Yulong Li <liyulong12@mails.ucas.ac.cn>

use strict; 
use warnings;
use Getopt::Long;
use Storable;
my ($buf, $id, $ctg_id, %seqs, $ctgs, %loci, $ref, $one, $help, $flag);
my ($in_file, $out_file, $format);
my $min_len = 0;
my $width   = 0;
my $retain  = '';
my $usage = '
    -i: input fasta file.
    -o: outfile name.
    -m: minimus contig length.
    -s: sub contig list (one per line).
    -R: Only reatin enzyme cut sites contigs, overhang pattern.
    -l: One contig per locus.
    -w: width for print fasta.
    -f: format of input file, "dat" or "fasta", default is dat from RADassembler.
    -h: Help message.
';

GetOptions (
			"in=s" 			=> \$in_file,
			"out=s" 		=> \$out_file,
            "min=i"         => \$min_len,
            "sublst=s"      => \$ctg_id,
            "Retain=s"      => \$retain,
            "locus:1"       => \$one,
            "width:i"       => \$width,
            "format=s"      => \$format,
			"help:1"		=> \$help
			)
			or die ("Error in command line arguments\n");
	
die "$0 $usage\n" if !($in_file && $out_file) or ($help);
$format  = ($format && ($format eq 'fasta' || $format eq 'fa')) ? 'fasta' : 'dat';
# parse sub list.  
if ($ctg_id) {
    
    $flag = 1; #flag for sublist.
	open(my $in, "$ctg_id");
	while (<$in>) {
		chomp;
        $_ =~ s/\r//g;
		my $ctg = $_;
		$ctgs->{$ctg} = $ctg;
	}
	close $in;
} else {
    $flag = 0;
}

# parse data.

if ($format eq 'fasta') {
    (my $seqs, $ref) = parse_fasta($in_file);
    %seqs = %$seqs;
} else {
    (my $seqs, $ref) = retrieve_dat($in_file);
    %seqs = %$seqs;
}
    
#Output results
open(my $out, ">$out_file") or die "$!";
open(my $stat, ">$out_file.stat") or die "$!";
my ($num_ctg, $num_loci, $sum, @lengths);
my %base = (
    "A"=>0,
    "T"=>0,
    "G"=>0,
    "C"=>0,
    "N"=>0
);
print STDERR "Now starting outputting...\n";
if ($one) {

    foreach my $locus (sort {$a cmp $b} keys %$ref) {
        $num_loci++;
        my $out_seq = '';
        my $out_id  = '';
        my @id  = @{$ref->{$locus}};
        if (@id > 1) {
            my $cnt  = 0;
            my (%tmp1, %tmp2);
            foreach (@id) {
                my $seq = $seqs{$_};
                my $len = length($seq);
                next if ($len < $min_len);
                next if ($retain && $seq !~ /$retain$/);
                if ($seq =~ /$retain$/) {
                    $cnt++;
                    $tmp2{$_} = $seq;
                }
                $tmp1{$_} = $seq;
            }
            next if scalar(keys %tmp1) == 0;
            if ($cnt == 1) {
                # one contig with cut sites.
                foreach (keys %tmp2) {
                    my $seq = $tmp2{$_};
                    print_fasta($out, $_, $seq);
                    $out_seq = $seq;
                    $out_id  = $_;
                }
            } elsif ($cnt > 1) {
                # multi contigs with cut sites.
                my ($id, $seq) = pick_max(\%tmp2);
                print_fasta($out, $id, $seq);
                $out_seq = $seq;
                $out_id  = $id;
            } else {
                my ($id, $seq) = pick_max(\%tmp1);
                print_fasta($out, $id, $seq);
                $out_seq = $seq;
                $out_id  = $id;
            }  
        } else {
            my $seq = $seqs{$id[0]};
            die "$id[0]\n" if !$seq;
            my $len = length($seq);
            next if ($len < $min_len);
            next if ($retain && $seq !~ /$retain$/);
            print_fasta($out, $id[0], $seq);
            $out_seq = $seq;
            $out_id  = $id[0];
        }
        # stat.
        my $len = length($out_seq);
        print $stat "$out_id\t$len\n";
        $sum += $len;
        push(@lengths, $len);
        my @curSeq = split(//, uc($out_seq));
        map {$base{$_}++;} @curSeq;
    }
    $num_ctg = $num_loci;
} else {
    my %loci;
    foreach $id (sort {$a cmp $b} keys %seqs) {
        my $tmp = $seqs{$id};
        my $len = length($tmp);
        next if ($len < $min_len);
        next if ($retain && $tmp !~ /$retain$/);
        print_fasta($out,$id, $tmp);
        print $stat "$id\t$len\n";
        $sum += $len; $num_ctg++;
        push(@lengths, $len);
        my @curSeq = split(//, uc($tmp));
        map {$base{$_}++;} @curSeq;
        my $locus = (split(/\./, $id))[0];
        $loci{$locus}++;
	}
    $num_loci = scalar(keys %loci);  
}

die "No sequence left!\n" unless $sum;
my $mean = int($sum / $num_ctg);
my $add    = 0;
my $num    = 0;
@lengths   = sort {$b <=> $a} @lengths;
while ($add < ($sum/2)) { $add += $lengths[$num++]; }
my $n50    = $lengths[$num-1];
my $median = $lengths[int($#lengths/2)];
my $gc     = sprintf "%.4f", (($base{'G'}+$base{'C'})/$sum * 100);
my $n      = sprintf "%.4f", ($base{'N'}/$sum * 100);
my $min    = $lengths[-1];
my $max    = $lengths[0];
print STDERR "Final retained $num_ctg contigs, $num_loci loci, $sum bases, GC: $gc%, N%: $n%,
      mean: $mean, median: $median, n50: $n50, max: $max, min: $min\n";

sub pick_max {
    # return a seq with maximum length.
    my $tmp  = shift;
    my %seqs = %$tmp;
    my $prev = 0;
    my $seq  = '';
    my $id   = '';
    foreach (keys %seqs) {
        my $tmp = $seqs{$_};
        my $len = length($tmp);
        next if $len <= $prev;
        $prev = $len;
        $id  = $_;
        $seq = $tmp;
    }
    return ($id, $seq);
}

sub print_fasta {

    my ($out, $id, $seq) = @_;

    my $tmp = $seq;
	
    print $out $id, "\n";
    if ($width) {
        # very slow, deprecated.
        while (length($tmp) > $width) {
            print $out substr($tmp, 0, $width), "\n";
            $tmp = substr($tmp, $width);
        }
    } else {
        print $out $tmp, "\n";
    }
}	

sub retrieve_dat {
    
    my $in_dat   = shift;
    my (%seqs, $ref);
    my $out_seqs = retrieve($in_dat) or die "$!";
    foreach my $locus (sort keys %{$out_seqs}) {
        next if $flag && !$ctgs->{$locus}; # white list for locus.
        my $seq = $out_seqs->{$locus};
        my $id  = '';
        my @lines = split(/\R/,$seq);
        foreach my $line (@lines) {
            next if $line =~ /^$/;
            if ($line =~ /^>/) {
                $line =~ s/^>/>$locus\./; # fasta name.
                $id   = $line;
                next;
            }
            next if $flag && !($ctgs->{$id}||$ctgs->{substr($id,1)}); # white list for contig.
            push @{$ref->{$locus}}, $id;
            $seqs{$id} .= $line;
        }
    }
    
    return (\%seqs, $ref)
}

sub parse_fasta {
    # parse fasta.
    my $in_file = shift;
    my (%seqs, $buf, $id, $locus, $ref);
    open(my $in, "$in_file") or die "$!";
    $buf = "";
    $locus = "";
    while (<$in>) {
        chomp;
        $_ =~ s/\R//g;
        if (/^>/) {
            if (length($buf) > 0) {
                $seqs{$id}  = $buf;
                $buf = "";
            }
            $id = $_;
            $locus = (split(/\./, $_, 2))[0]; # Only for RADassembler.
            $locus =~ s/>//; # Locus name.
        } else {
            next if $flag && !($ctgs->{$id}||$ctgs->{$locus}||$ctgs->{substr($id,1)}); # white list for locus or contig name 
            push @{$ref->{$locus}}, $id;
            $buf .= $_;
        }
    }
    #The last one
    if ($flag && !($ctgs->{$id}||$ctgs->{$locus}||$ctgs->{substr($id,1)})) { 
    
    } else {
        $seqs{$id}  = $buf;
    }
    close $in;
    return (\%seqs, $ref);
}

