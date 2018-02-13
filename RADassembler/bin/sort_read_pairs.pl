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
# Modified by Yulong Li, new version
#

use strict;
use constant stacks_version => "1.32.2";
use Array::Shuffle 'shuffle_array'; # faster.
use constant true  => 1;
use constant false => 0;

my $debug          = 0;
my $white_list     = "";
my $cat_white_list = "";
my $in_path        = "";
my $out_path       = "";
my $samp_path      = "";
my $out_type       = "fasta";
my $in_type        = "fastq";
my $gzipped        = false;
my $depth;
my $mincov;
parse_command_line();
if (!-d $out_path) {
	mkdir $out_path;
} 
=cut
else {
    print STDERR "Error: output directory: $out_path already exist, \ndo you want to delete it ? y/n\n";
    my $info = <STDIN>;
    $info =~ s/[\r\n]//g;
    if ($info eq 'y' || $info eq 'yes') {
        `rm -rf $out_path`;
        mkdir $out_path;
    } else {
        print STDERR "Please rename your directory.\n";
        exit(1);
    }
    
}
=cut

my (@files, %matches, %stacks, %reads, %marker_wl);

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
    load_matches($in_path, $file, \%matches, \%marker_wl);

    $i++;
}

#
# Determine which catalog loci have more than a single match from each sample and blacklist them.
#
print STDERR "Identifying catalog loci that have more than one match from a single sample...";
my %multiple_matches;
check_mult_catalog_matches(\%matches, \%multiple_matches);
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
my ($cat_id, $path);
=cut
foreach $cat_id (keys %matches) {
    #
    # Check that this catalog ID only has a single match from each sample.
    #
    next if (defined($multiple_matches{$cat_id}));

    $path  = $out_path . "/" . $cat_id;
    $path .= $out_type eq "fasta" ? ".fa" : ".fq";

    if (-e $path) {
	die("Error: output files already exist. This program will append data to files if\n" .
	    "they already exist. Please delete these files and re-execute sort_read_pairs.pl.\n");
    }
}
=cut
$i = 1;
foreach $file (@files) {
    printf(STDERR "Processing file % 2s of % 2s [%s]\n", $i, $num_files, $file->{'prefix'});
    #
    # Load the ustacks tag file for each sample, containing the read-to-stack mappings
    #
    $stacks{$file->{'prefix'}} = {};

    print STDERR "  Loading tag file...";
    load_stacks($in_path, $file, $stacks{$file->{'prefix'}});
    print STDERR "done.\n";

    #
    # Map the read-pairs to the stack/catalog match they correspond to.
    #
    $reads{$file->{'prefix'}} = {};

    print STDERR "  Loading sample file...";
    $in_type eq "fastq" ?
	process_fastq_read_pairs($samp_path, $file, \%stacks, $reads{$file->{'prefix'}}) :
	process_fasta_read_pairs($samp_path, $file, \%stacks, $reads{$file->{'prefix'}});
    print STDERR "done.\n";
	
	print STDERR "  Clearing memory...";
	undef(%{$stacks{$file->{'prefix'}}});
	print STDERR "done.\n";
	
    $i++;
}

#
# Output files (multi-individuals) for each loci.
#
print STDERR "Printing results...";
print_results($out_path, \%matches, \%stacks, \%reads, \%multiple_matches, $depth);
print STDERR "done.\n";

#
# Clean up memory usage.
#

#print STDERR "  Clearing memory...";
#undef(%{$stacks{$file->{'prefix'}}});
#undef(%{$reads{$file->{'prefix'}}});
#print STDERR "  done.\n";
	
sub load_matches {
    my ($in_path, $in_file, $matches, $marker_wl) = @_;

    my ($file, $in_fh, $line, @parts, $key);

    if ($gzipped == true) {
	$file  = $in_path . "/" . $in_file->{'prefix'} . $in_file->{'suffix'} . ".matches.tsv.gz";
	open($in_fh, "gunzip -c $file |") or die("Unable to open catalog matches file '$file', $!\n");
    } else {
	$file  = $in_path . "/" . $in_file->{'prefix'} . $in_file->{'suffix'} . ".matches.tsv";
	open($in_fh, "<$file") or die("Unable to open catalog matches file '$file', $!\n");
    }

    while ($line = <$in_fh>) {
        chomp $line;
        @parts = split(/\t/, $line);

        if (length($cat_white_list) > 0) {
	    next if (!defined($marker_wl->{$parts[2]}));
	}

        if (!defined($matches->{$parts[2]})) {
            $matches->{$parts[2]} = {};
        }

        #
        # Index by catalog_ID -> sample_ID|stack_ID
        #
        $key = $in_file->{'prefix'} . "|" . $parts[4];
        $matches->{$parts[2]}->{$key}++;
    }

    close($in_fh);
}

sub load_stacks {
    my ($in_path, $in_file, $stacks) = @_;

    my ($file, $in_fh, $line, @parts);

    if ($gzipped == true) {
	$file = $in_path . "/" . $in_file->{'prefix'} . $in_file->{'suffix'} . ".tags.tsv.gz";
	open($in_fh, "gunzip -c $file |") or die("Unable to open '$file', $!\n");
    } else {
	$file = $in_path . "/" . $in_file->{'prefix'} . $in_file->{'suffix'} . ".tags.tsv";
	open($in_fh, "<$file") or die("Unable to open '$file', $!\n");
    }

    while ($line = <$in_fh>) {
        chomp $line;
        @parts = split(/\t/, $line);

        next if ($parts[6] eq "consensus" || $parts[6] eq "model" || $parts[6] eq "secondary");

        #
        # Index by sequence ID -> stack ID
        #
		if ($parts[8] =~ /(.+)\s+(.+)$/) {
			#
			# Type I: 
			# HWI-ST1276:71:C1162ACXX:1:1101:1208:2458 1:N:0:CGATGT
			#
			
			$stacks->{$1} = $parts[2];
		} elsif ($parts[8] =~ /(.+)[12]$/) {
            # Type II:
			# CTACAG_8_1103_15496_190439_1|1
			#
			# Type III:
			# 4_1101_13393_1801_1
			# 4_1101_13393_1801/1
            # ...
            # print $1, "\n";
			$stacks->{$1} = $parts[2];
		} else {
			$stacks->{$parts[8]} = $parts[2];
		}
    }

    close($in_fh);
}

sub process_fastq_read_pairs {
    my ($in_path, $in_file, $stacks, $reads) = @_;

    my ($file, $in_fh, $line, $seq, $qual, $key, $read_id, $read_dd);

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
	next if (substr($line, 0, 1) ne "@");
	chomp $line;
	#
	#
	#
	if ($line =~ /^@(.+)\s+(.+)$/) {
		#
		# Type I: 
		# @HWI-ST1276:71:C1162ACXX:1:1101:1208:2458 1:N:0:CGATGT
		#
		$read_id = $1;	  #Modified
		$read_dd = $line; #Original read id.
	} elsif ($line =~ /^@(.+)[12]$/) {
        # Type II:
		# @CTACAG_8_1103_15496_190439_1|1
		#
		# Type III:
		# @4_1101_13393_1801_1
        # @4_1101_13393_1801/1
		# ...
        # print $1, "\n";
		$read_id = $1;
		$read_dd = $line;
	} else {
		$read_id = substr($1, 1);
		$read_dd = $read_id;
	}
	#
	#
	#
	
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

    if (!defined($reads->{$key})) {
        $reads->{$key} = [];
    }

    push(@{$reads->{$key}}, {'id' => $read_dd, 'seq' => $seq, 'qual' => $qual}); #
    }
}

sub process_fasta_read_pairs {
    my ($in_path, $in_file, $stacks, $reads) = @_;

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
	
    #
	#
	#
	if ($line =~ /^>(.+)\s+(.+)$/) {
		#
		# Type I: 
		# >HWI-ST1276:71:C1162ACXX:1:1101:1208:2458 1:N:0:CGATGT
		#
		$read_id = $1;	  #Modified
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
	#
	#
	#
	$seq = <$in_fh>;
	chomp $seq;

        $key = $stacks->{$in_file->{'prefix'}}->{$read_id};

        next if (!defined($key));

        if (!defined($reads->{$key})) {
            $reads->{$key} = [];
        }

        push(@{$reads->{$key}}, {'id' => $read_dd, 'seq' => $seq, 'qual' => ""}); #
    }
}

sub print_results {
    my ($out_path, $matches, $stacks, $reads, $multiple_matches, $depth) = @_;

    my ($path, $cat_id, $sample, $stack_id, $read, $out_fh, $i, @keys, $count, $key, $mult_hits, $mindepth, $maxdepth);
	
	if ($depth =~ /(\d+)\:(\d+)/) {
		$mindepth = $1;
		$maxdepth = $2;
	} else {
		$mindepth = $depth;
		$maxdepth = 0;
	}
	my $tot = 0; # total num of loci.
    my $tod = 0; # total depth.
    # 
    # If a catalog ID matches stacks from multiple samples, print them out together.
    #
    foreach $cat_id (keys %{$matches}) {
        #
        # Check that this catalog ID only has a single match from each sample.
        #
        next if (defined($multiple_matches->{$cat_id}));
		
		my $cnt = 0; # Depth for a loci(single-end).
		my $sam = scalar (keys %{$matches->{$cat_id}}); # sample num for a locus.
		my @out_seq; # Array holds output seqs for a catalog loci.
		
        next if $mincov && $sam < $mincov;
		if ($mindepth) {
			#
			# a minimum depth for paired-end required...
			# 
			#
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
							$out .= ">". $cat_id. "|". $sample. "|". $stack_id. "|". $read->{'id'}. "\n";
							$out .= $read->{'seq'} . "\n";
						} else {
							$out .= "@". $cat_id. "|". $sample. "|". $stack_id. "|". $read->{'id'}. "\n";
							$out .= $read->{'seq'}. "\n";
							$out .= "+\n";
							$out .= $read->{'qual'}. "\n";
						}
						push @out_seq, $out;
					}
				}
				next if $cnt < $mindepth;
                $tot++;
				$path  = $out_path . "/" . $cat_id;
				$path .= $out_type eq "fasta" ? ".fa" : ".fq";
				open($out_fh, ">$path") or die("Unable to open $path; '$!'\n"); # all samples together.
					
				if ($cnt > $maxdepth) {
					shuffle_array @out_seq; # need rewrite here.
					map {print $out_fh $_} @out_seq[1..$maxdepth];
					close $out_fh;
					undef @out_seq;
                    $tod += $maxdepth;
					next;
				} else {
					map {print $out_fh $_} @out_seq;
					close $out_fh;
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
		$path  = $out_path . "/" . $cat_id;
		$path .= $out_type eq "fasta" ? ".fa" : ".fq";
        open($out_fh, ">>$path") or die("Unable to open $path; '$!'\n"); #single sample.	
        foreach $key (keys %{$matches->{$cat_id}}) {
				
			($sample, $stack_id) = split(/\|/, $key);
            $tod += scalar(@{$reads->{$sample}->{$stack_id}}); # read depth if no threshold is set.
            foreach $read (@{$reads->{$sample}->{$stack_id}}) {
				if ($out_type eq "fasta") {
					print $out_fh
					">", $cat_id, "|", $sample, "|", $stack_id, "|", $read->{'id'}, "\n",
					$read->{'seq'}, "\n";
				} else {
					print $out_fh
					"@", $cat_id, "|", $sample, "|", $stack_id, "|", $read->{'id'}, "\n",
					$read->{'seq'}, "\n",
					"+\n",
					$read->{'qual'}, "\n";
				}
            }
        }

        close($out_fh);
    }
    print STDERR "\nA total of $tot loci were exported, with a mean depth of ". sprintf "%.2f.\n", $tod/$tot;
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
	next if ($line =~ /batch_\d+\.catalog/);

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
	elsif ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-m$/) { $depth      = shift @ARGV; }
	elsif ($_ =~ /^-c$/) { $mincov     = shift @ARGV; }
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
    h: display this help message.
    d: turn on debug output.

EOQ

exit(0);
}
