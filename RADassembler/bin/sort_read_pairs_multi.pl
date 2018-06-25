#!/usr/bin/env perl
#
# Copyright 2011-2015, Julian Catchen <jcatchen@illinois.edu>
#
# This file is part of Stacks.
#
# Stacks is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Stacks is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
#

#
# Sort paired-end sequences according to the stacks the non-paired-end 
# was found in. 
#
# By Julian Catchen <jcatchen@illinois.edu>
#
#
# Modified by Yulong Li
#

use strict;
use warnings;
use constant stacks_version => "180528";
use threads;
use Storable;
use Array::Shuffle 'shuffle_array'; # faster.
use constant true  => 1;
use constant false => 0;
#use Data::Dumper;

my $debug          = 0;
my $white_list     = "";
my $cat_white_list = "";
my $in_path        = "";
my $out_path       = "";
my $samp_path      = "";
my $out_type       = "fasta";
my $in_type        = "fastq";
my $gzipped        = false;
my $threads        = 1;
my $verify         = 1;
my ($depth, $mincov);
parse_command_line();
if (!-d $out_path) {
    mkdir $out_path;
} 

my (@files, %matches, %stacks, %marker_wl, $reads, %stacks_all);

my ($mindepth, $maxdepth);

if ($depth && $depth =~ /(\d+)\:(\d+)/) {
    $mindepth = $1;
    $maxdepth = $2;
} else {
    $mindepth = $depth;
    $maxdepth = 0;
}

build_file_list(\@files);

my ($file, $num_files, $i, $key);

if (length($cat_white_list) > 0) {
    load_white_list($cat_white_list, \%marker_wl);
    print STDERR "Loaded ", scalar(keys %marker_wl), " catalog IDs from '$cat_white_list'\n";
}

$num_files = scalar(@files);
$i         = 1;
foreach $file (@files) {
    printf(STDERR "Loading catalog matches, file % 2s of % 2s [%s]\n", $i, $num_files, $file->{'prefix'});
    #
    # Load the sstacks file, listing matches between the stacks and the catalog
    #
    load_matches($in_path, $file, \%matches, \%marker_wl, \%stacks_all);

    $i++;
}

#
# Determine which catalog loci have more than a single match from each sample and blacklist them.
#
my %multiple_matches;
if ($verify) {
print STDERR "Identifying catalog loci that have more than one match from a single sample...";
check_mult_catalog_matches(\%matches, \%multiple_matches);
}
    
print STDERR 
    "done\n",
    "  These loci will be excluded when collating paired-end reads;\n",
    "  A list of them has been recorded: $out_path/sort_read_pairs.log\n",
    scalar(keys %matches), " total catalog loci; ", scalar(keys %multiple_matches), 
    " will be excluded, processing ", scalar(keys %matches) - scalar(keys %multiple_matches), " loci.\n";


#
# Check if files already exist, if so, exit to prevent adding data to existing files.
#
# not needed.

#
# reads->sample_id->stack_id->[]
#
my ($cat_id, $path);
$i = 1;

foreach $file (@files) {
    
    printf(STDERR "Processing file % 2s of % 2s [%s]\n", $i, $num_files, $file->{'prefix'});

    my $t = threads->create(\&parse_tag_and_reads, $file);
    $i++;
    my $j = scalar(threads->list(threads::all));
    if ($j < $threads && $i < scalar(@files)) { next; }
    my @pool = threads->list(threads::all);
    while (scalar(@pool) >= $threads) {
        foreach(@pool) {
            if ($_->is_joinable()) {
                $_->join();
            }
        }
        @pool = threads->list(threads::all);
    }
}
foreach (threads->list(threads::all)) {$_->join();}

#
# Combine reads data structure.
#
print STDERR "Combining reads...";
cmb_reads(\@files);
print STDERR "done.\n";

#
# Output files (multi-individuals) for each loci.
#
print STDERR "Printing results...";
print_results($out_path, \%matches, \%stacks, $reads, \%multiple_matches, $depth);

#
# Clean up memory usage.
#



sub parse_tag_and_reads {
    
    my $file  = shift;
    #
    # Load the ustacks tag file for each sample, containing the read-to-stack mappings
    #
    $stacks{$file->{'prefix'}} = {};

    load_stacks($in_path, $file, $stacks{$file->{'prefix'}});
    
    #
    # Map the read-pairs to the stack/catalog match they correspond to.
    #
    
    
    $in_type eq "fastq" ?
    process_fastq_read_pairs($samp_path, $file, \%stacks, $reads) :
    process_fasta_read_pairs($samp_path, $file, \%stacks, $reads);
    
    undef(%{$stacks{$file->{'prefix'}}});
}

sub cmb_reads {
    my $files = shift;
    foreach my $file (@$files) {
        my $name = $file->{'prefix'};
        $reads->{$name} = retrieve("$out_path/$name");
        `rm $out_path/$name`;
    }
    
}

sub load_matches {
    my ($in_path, $in_file, $matches, $marker_wl) = @_;

    my ($file, $in_fh, $line, @parts, $key, $v);

    if ($gzipped == true) {
    $file  = $in_path . "/" . $in_file->{'prefix'} . $in_file->{'suffix'} . ".matches.tsv.gz";
    open($in_fh, "gunzip -c $file |") or die("Unable to open catalog matches file '$file', $!\n");
    } else {
    $file  = $in_path . "/" . $in_file->{'prefix'} . $in_file->{'suffix'} . ".matches.tsv";
    open($in_fh, "<$file") or die("Unable to open catalog matches file '$file', $!\n");
    }

    while ($line = <$in_fh>) {
        next if $line =~ /^#/;
        chomp $line;
        @parts = split(/\t/, $line);
        my ($cat_id, $loc_id);
        if (!$v) {
            # only the use the first line.
            my $npart = @parts;
            if ($npart == 8) {$v = 1;}
            elsif ($npart == 6) {$v = 2;}
            else {die "Stacks tags files error!";}
        }
        if ($v == 1) {
            $cat_id = $parts[2];
            $loc_id = $parts[4];
        } elsif ($v == 2) {
            $cat_id = $parts[0];
            $loc_id = $parts[2];
        } else {
            die "Stacks files error!";
        }
        if (length($cat_white_list) > 0) {
        next if (!defined($marker_wl->{$cat_id}));
    }

        if (!defined($matches->{$cat_id})) {
            $matches->{$cat_id} = {};
        }

        #
        # Index by catalog_ID -> sample_ID|stack_ID
        #
        $key = $in_file->{'prefix'} . "|" . $loc_id;
        $matches->{$cat_id}->{$key}++;
        $stacks_all{$key}++;
    }

    close($in_fh);
}

sub load_stacks {
    my ($in_path, $in_file, $stacks) = @_;

    my ($file, $in_fh, $line, @parts, $v);

    if ($gzipped == true) {
    $file = $in_path . "/" . $in_file->{'prefix'} . $in_file->{'suffix'} . ".tags.tsv.gz";
    open($in_fh, "gunzip -c $file |") or die("Unable to open '$file', $!\n");
    } else {
    $file = $in_path . "/" . $in_file->{'prefix'} . $in_file->{'suffix'} . ".tags.tsv";
    open($in_fh, "<$file") or die("Unable to open '$file', $!\n");
    }

    while ($line = <$in_fh>) {
        next if $line =~ /^#/;
        chomp $line;
        @parts = split(/\t/, $line);
        my ($model, $seq_id, $key, $loc_id);
        if (!$v) {
            # only the use the first line.
            my $npart = @parts;
            if ($npart == 14) {$v = 1;}
            elsif ($npart == 9) {$v = 2;}
            else {die "Stacks tags files error!";}
        }
        if ($v == 1) {
            $model  = $parts[6];
            $seq_id = $parts[8];
            $loc_id = $parts[2];
        } elsif ($v == 2) {
            $model  = $parts[2];
            $seq_id = $parts[4];
            $loc_id = $parts[1];
        } else {
            die "Stacks tags files error!";
        }
        next if ($model ne 'primary');
        $key = $in_file->{'prefix'} . "|" . $loc_id;
        next if not defined $stacks_all{$key};
        #
        # Index by sequence ID -> stack ID
        #
        if ($seq_id =~ /(.+)\s+(.+)$/) {
            #
            # Type I: 
            # HWI-ST1276:71:C1162ACXX:1:1101:1208:2458 1:N:0:CGATGT
            #
            
            $stacks->{$1} = $loc_id;
        } elsif ($seq_id =~ /(.+)[12]$/) {
            # Type II:
            # CTACAG_8_1103_15496_190439_1|1
            #
            # Type III:
            # 4_1101_13393_1801_1
            # 4_1101_13393_1801/1
            #
            $stacks->{$1} = $loc_id;
        } else {
            $stacks->{$seq_id} = $loc_id;
        }
    }

    close($in_fh);
}

sub process_fastq_read_pairs {
    my ($in_path, $in_file, $stacks, $a) = @_;
    my $local_reads;
    my ($file, $in_fh, $line, $seq, $qual, $key, $read_id, $read_dd, %size);

    if ($gzipped == true) {
        $file = $in_path . "/" . $in_file->{'prefix'} . "_1.fq.gz";
        if (-e $file) {
            open($in_fh, "gunzip -c $file |") or die("Unable to open paired-end input file $file\n");
        } else {
            $file = $in_path . "/" . $in_file->{'prefix'} . '_2.fq.gz';
            open($in_fh, "gunzip -c $file |") or die("Unable to open paired-end input file $file\n");
        }
    } else {
        $file = $in_path . "/" . $in_file->{'prefix'} . "_1.fq";
        if (-e $file) {
            open($in_fh, "<$file") or die("Unable to open paired-end input file $file\n");
        } else {
            $file = $in_path . "/" . $in_file->{'prefix'} . '_2.fq';
            open($in_fh, "<$file") or die("Unable to open paired-end input file $file\n");
        }
    }
    
    while ($line = <$in_fh>) {
        next if ($line !~ /^@/);
        chomp $line;
        if ($line =~ /^@(.+)\s+(.+)$/) {
            #
            # Type I: 
            # @HWI-ST1276:71:C1162ACXX:1:1101:1208:2458 1:N:0:CGATGT
            #
            $read_id = $1;      #Modified
            $read_dd = $line; #Original read id.
        } elsif ($line =~ /^@(.+)[12]$/) {
            # Type II:
            # @CTACAG_8_1103_15496_190439_1|1
            #
            # Type III:
            # @4_1101_13393_1801_1
            # @4_1101_13393_1801/1
            # 
            $read_id = $1;
            $read_dd = $line;
        } else {
            $read_id = substr($1, 1);
            $read_dd = $read_id;
        }
        $seq = <$in_fh>;
        chomp $seq;
        #
        # Read the repeated ID and the quality scores.
        #
        <$in_fh>;
        $qual = <$in_fh>;
        chomp $qual;
        $key = $stacks->{$in_file->{'prefix'}}->{$read_id};
        next if (!defined($key));
        
        if (!defined($local_reads->{$key})) {
            $local_reads->{$key} = [];
        }
        
        #
        # sample->stack_id->seqs.
        #
        #$size{$key}++;
        #next if ($maxdepth && $size{$key} > $maxdepth); # reduce the hash size.
        my @tmp = ($read_dd, $seq, $qual);
        push(@{$local_reads->{$key}}, \@tmp);
        
    }
    my $name = $in_file->{'prefix'};
    store $local_reads, "$out_path/$name";
    undef $local_reads;
}

sub process_fasta_read_pairs {
    my ($in_path, $in_file, $stacks, $reads) = @_;
    my $local_reads;
    my ($file, $in_fh, $line, $seq, $qual, $key, $read_id, $read_dd);

    if ($gzipped == true) {
        $file = $in_path . "/" . $in_file->{'prefix'} . "_1.fa.gz";
       if (-e $file) {
            open($in_fh, "gunzip -c $file |") or die("Unable to open paired-end input file $file\n");
        } else {
            $file = $in_path . "/" . $in_file->{'prefix'} . '_2.fa.gz';
            open($in_fh, "gunzip -c $file |") or die("Unable to open paired-end input file $file\n");
        }
    } else {
        $file = $in_path . "/" . $in_file->{'prefix'} . "_1.fa";
        if (-e $file) {
            open($in_fh, "$file") or die("Unable to open paired-end input file $file\n");
        } else {
            $file = $in_path . "/" . $in_file->{'prefix'} . '_2.fa';
            open($in_fh, "$file") or die("Unable to open paired-end input file $file\n");
        }
    }
    while ($line = <$in_fh>) {
    next if (substr($line, 0, 1) ne ">");
    chomp $line;
    if ($line =~ /^>(.+)\s+(.+)$/) {
        #
        # Type I: 
        # >HWI-ST1276:71:C1162ACXX:1:1101:1208:2458 1:N:0:CGATGT
        #
        $read_id = $1;      #Modified
        $read_dd = $line; #Original read id.
    } elsif ($line =~ /^>(.+)[12]$/) {
        # Type II:
        # >CTACAG_8_1103_15496_190439_1|1
        #
        #
        # Type III:
        # >4_1101_13393_1801_1
        #
        $read_id = $1;
        $read_dd = $line;
    } else {
        $read_id = substr($1, 1);
        $read_dd = $read_id;
    }

    $seq = <$in_fh>;
    chomp $seq;
    $key = $stacks->{$in_file->{'prefix'}}->{$read_id};

    next if (!defined($key));
    
    if (!defined($local_reads->{$key})) {
        $local_reads->{$key} = [];
    }
    my @tmp = ($read_dd, $seq, '');
    push(@{$local_reads->{$key}}, \@tmp); 
    }
    my $name = $in_file->{'prefix'};
    store $local_reads, "$out_path/$name";
    undef $local_reads;
}

sub print_results {
    my ($out_path, $matches, $stacks, $reads, $multiple_matches, $depth) = @_;

    my ($path, $cat_id, $sample, $stack_id, $read, $out_fh, $i, @keys, $count, $key, $mult_hits);
    
    my $tot = 0; # total num of loci.
    my $tod = 0; # total depth.
    my $out_seqs = {};
    # 
    # If a catalog ID matches stacks from multiple samples, print them out together.
    #
    foreach $cat_id (keys %{$matches}) {
        #
        # Check that this catalog ID only has a single match from each sample.
        #
        next if (defined($multiple_matches->{$cat_id}));
        
        my $cnt = 0; # Depth for a loci(single-end).
        my @out_seq; # Array holds output seqs for a catalog loci.
        if ($mincov) {
            # a safe method in case of no multiple hit verify.
            my %sam;
            foreach my $key (keys %{$matches->{$cat_id}}) {
                ($sample, $stack_id) = split(/\|/, $key);
                $sam{$sample}++;
            }
            next if scalar(keys %sam) < $mincov;
        }
        if ($mindepth) {
            if ($maxdepth) {
                #
                # randomly downscale to the maximum depth 
                # only applicable to single-end
                #
                foreach my $key (keys %{$matches->{$cat_id}}) {
                    
                    ($sample, $stack_id) = split(/\|/, $key);
                    foreach $read (@{$reads->{$sample}->{$stack_id}}) {
                        $cnt++;
                        my $out = '';
                        if ($out_type eq "fasta") {
                            $out .= ">". $cat_id. "|". $sample. "|". $stack_id. "|". $$read[0]. "\n";
                            $out .= $$read[1] . "\n";
                        } else {
                            $out .= "@". $cat_id. "|". $sample. "|". $stack_id. "|". $$read[0]. "\n";
                            $out .= $$read[1]. "\n";
                            $out .= "+\n";
                            $out .= $$read[2]. "\n";
                        }
                        push @out_seq, $out;
                    }
                }
                next if $cnt < $mindepth;
                $tot++;
                if ($cnt > $maxdepth) {
                    shuffle_array @out_seq; # need rewrite here.
                    $out_seqs->{$cat_id} = join('', @out_seq[1..$maxdepth]);
                    undef @out_seq;
                    $tod += $maxdepth;
                    next;
                } else {
                    $out_seqs->{$cat_id} = join('', @out_seq);
                    undef @out_seq;
                    $tod += $cnt;
                    next;
                }
                
            } else {
                foreach my $key (keys %{$matches->{$cat_id}}) {
                    ($sample, $stack_id) = split(/\|/, $key);
                    map {$cnt++} @{$reads->{$sample}->{$stack_id}};
                }
                next if $cnt < $mindepth;
            }
        }
        
        $tot++;
        my $out = '';
        foreach $key (keys %{$matches->{$cat_id}}) {
                
            ($sample, $stack_id) = split(/\|/, $key);
            $tod += scalar(@{$reads->{$sample}->{$stack_id}}); # read depth if no threshold is set.
            foreach $read (@{$reads->{$sample}->{$stack_id}}) {
                if ($out_type eq "fasta") {
                    $out .= ">". $cat_id. "|". $sample. "|". $stack_id. "|". $$read[0]. "\n";
                    $out .= $$read[1] . "\n";
                } else {
                    $out .= "@". $cat_id. "|". $sample. "|". $stack_id. "|". $$read[0]. "\n";
                    $out .= $$read[1]. "\n";
                    $out .= "+\n";
                    $out .= $$read[2]. "\n";
                }
            }

        }

        $out_seqs->{$cat_id} = $out;
    }
    store $out_seqs, "$out_path/exported_reads.dat";
    print STDERR "done.\n";
    print STDERR "A total of $tot loci were exported, with a mean depth of ". sprintf "%.2f.\n", $tod/$tot;
}

sub check_mult_catalog_matches {
    my ($matches, $multiple_matches) = @_;

    my ($fh, $key, $sample, $stack_id);

    #
    # Find catalog loci that have more than a single match from one or more samples
    # and log those loci that will be excluded.
    #
    open($fh, ">$out_path/sort_read_pairs.log") or die("Unable to open log file, $!\n");
    print $fh 
    "# The catalog loci listed below have more than a single match from one or more individuals, indicating undermerged or repetitive loci.\n",
    "# CatalogLocus <tab> Sample1:Locus1,Locus2;Sample2:Locus1,Locus2\n";

    foreach $cat_id (keys %{$matches}) {
    my %samples;

    foreach $key (keys %{$matches->{$cat_id}}) {
        ($sample, $stack_id) = split(/\|/, $key);
        push(@{$samples{$sample}}, $stack_id);
    }

    my $mult_hits = 0;
    my $str = "";

    foreach $sample (keys %samples) {
        if (scalar(@{$samples{$sample}}) > 1) {
        $mult_hits++;
        $str .= $sample . ":" . join(",", @{$samples{$sample}}) . "; ";
        }
    }
    
    if ($mult_hits > 0) {
        print $fh $cat_id, "\t", substr($str, 0, -1), "\n";
        $multiple_matches{$cat_id}++;
    }
    }

    close($fh);
}

sub count_reads {
    my ($catalog, $reads) = @_;

    my ($count, $key, $sample, $stack_id);

    $count = 0;

    foreach $key (keys %{$catalog}) {
        ($sample, $stack_id) = split(/\|/, $key);

    if (defined($reads->{$sample}->{$stack_id})) {
        $count += scalar(@{$reads->{$sample}->{$stack_id}});
    }
    }

    return $count;
}

sub build_file_list {
    my ($files) = @_;

    my (@ls, $line, $file, $prefix, $suffix);

    # Load a white list of files to process if it is supplied.
    my %wl;
    if (length($white_list) > 0) {
    load_white_list($white_list, \%wl);
    print STDERR "Loaded ", scalar(keys %wl), " filenames from '$white_list'\n";
    }

    @ls = glob("$in_path/*.tags.tsv");
    if (scalar @ls == 0) {
    @ls = glob("$in_path/*.tags.tsv.gz");
    $gzipped = true if (scalar @ls > 0);
    }

    foreach $line (@ls) {
    chomp $line;

    next if (length($line) == 0);    
    next if ($line =~ /.*catalog\..+/);

    ($file) = ($line =~ /$in_path\/(.+)\.tags\.tsv\.?g?z?$/); 

        if (length($white_list) > 0) {
        next if (!defined($wl{$file}));
    }

    if ($file =~ /\_1$/) {
        ($prefix, $suffix) = ($file =~ /^(.+)(\_1)$/);
    } else {
        $prefix = $file;
        $suffix = "";
    }

    push(@{$files}, {'prefix' => $prefix, 'suffix' => $suffix});
    }
}

sub load_white_list {
    my ($list, $wl) = @_;

    open(WHITE, "<" . $list) 
    or die("Unable to open white list file '$white_list': $!\n");

    my $line   = "";

    while ($line = <WHITE>) {
    chomp $line;

    next if (length($line) == 0);
    next if ($line =~ /^\s*#/);

    $wl->{$line}++;
    }

    close(WHITE);
}

sub parse_command_line {
    while (@ARGV) {
    $_ = shift @ARGV;
    if    ($_ =~ /^-p$/) { $in_path    = shift @ARGV; }
    elsif ($_ =~ /^-o$/) { $out_path   = shift @ARGV; }
    elsif ($_ =~ /^-s$/) { $samp_path  = shift @ARGV; }
    elsif ($_ =~ /^-t$/) { $out_type   = shift @ARGV; }
    elsif ($_ =~ /^-i$/) { $in_type    = shift @ARGV; }
    elsif ($_ =~ /^-W$/) { $white_list = shift @ARGV; }
    elsif ($_ =~ /^-w$/) { $cat_white_list = shift @ARGV; }
    elsif ($_ =~ /^-d$/) { $debug++;                  }
    elsif ($_ =~ /^-m$/) { $depth      = shift @ARGV; }
    elsif ($_ =~ /^-c$/) { $mincov     = shift @ARGV; }
    elsif ($_ =~ /^-T$/) { $threads    = shift @ARGV; }
    elsif ($_ =~ /^-x$/) { $verify     = 0; }
    elsif ($_ =~ /^-v$/) { version(); exit(); }
    elsif ($_ =~ /^-h$/) { usage(); }
    else {
        print STDERR "Unknown command line option: '$_'\n";
        usage();
    }
    }

    if ($out_type ne "fasta" && $out_type ne "fastq") {
    print STDERR "Output type must be either 'fasta' or 'fastq'.\n";
    usage();
    }

    if ($in_type ne "fasta" && $in_type ne "fastq") {
    print STDERR "Input type must be either 'fasta' or 'fastq'.\n";
    usage();
    }

    if (length($in_path) == 0) {
    print STDERR "You must specify a path to the Stacks output files.\n";
    usage();
    }

    if (length($out_path) == 0) {
    print STDERR "You must specify a path to write the collated output files.\n";
    usage();
    }

    if (length($samp_path) == 0) {
    print STDERR "You must specify a path to the paired-end reads.\n";
    usage();
    }

    $in_path   = substr($in_path, 0, -1)   if (substr($in_path, -1)   eq "/");
    $out_path  = substr($out_path, 0, -1)  if (substr($out_path, -1)  eq "/");
    $samp_path = substr($samp_path, 0, -1) if (substr($samp_path, -1) eq "/");
}

sub version {
    print STDERR "sort_read_pairs.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ; 
sort_read_pairs.pl -p path -s path -o path [-t type] [-W white_list] [-w white_list] [-d] [-h]
    p: path to the Stacks output files.
    s: path to paired-end sample files.
    o: path to output the collated FASTA files.
    i: input type, either 'fasta' or 'fastq' (default).
    t: output type, either 'fasta' (default) or 'fastq'.
    m: minimum depth for a locus (can also be mindepth:maxdepth, say 10:400, the minimum depth is 10, and
    the reads of a locus will also be randomly downscale to the maxdepth 400).
    c: coverage (num of individuals) for a locus
    W: a white list of files to process in the input path.
    w: a white list of catalog IDs to include.
    T: number of threads.
    x: do not verify multiple matches.
    h: display this help message.
    d: turn on debug output.

EOQ

exit(0);
}
