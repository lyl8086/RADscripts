#!/usr/bin/env perl
#@ Yulong Li
use strict; 
use warnings; 
use Parallel::ForkManager;
use Getopt::Long;
use Storable;
use constant V => 180731;

my (@files, $file, $in_tags, $sub_lst, $sec_asmb);
my ($in_path, $out_path, $loci, $out_seqs, $time);
my $num_threads = 1;
my $assembler   = 'cap3';
my $kmer        = 27;
my $cmd         = '';    ##Parameters parsed to CAP3 for the fist assembly.
my $collect     = 0;
my $help        = 0;
my $dat         = 0;
my $usage = 
'    
    Options [default values]:
    
    -i input path of fasta files
    -o output path for assembly files
    -l catalog file from stacks [*catalog.tags.tsv.gz]
    -t number of threads [1]
    -c collect the assembly results for the second reads [0]
    -p parameters parsed to CAP3
    -f sub file list for assembly
    -2 continue to run the second assembly
    -d use storable data structure.
    -a assembler, either "cap3" or "velvet".
    -k hash length for velvet, default 27.
    -h help
';

GetOptions ("par=s"         => \$cmd,
            "in=s"          => \$in_path,
            "out=s"         => \$out_path,
            "loci=s"        => \$loci,
            "threads=i"     => \$num_threads,
            "collect:1"     => \$collect,
            "file:s"        => \$sub_lst,
            "2:1"           => \$sec_asmb,
            "d:1"           => \$dat,
            "assembler:s"   => \$assembler,
            "kmer:i"        => \$kmer,
            "help:1"        => \$help
            ) 
or die ("Error in command line arguments\n");
    
if (!($in_path && $out_path && $loci) or $help) {print "$usage\n";exit(0)};

## If the outpath does not exist, create it;
if (! -e "$out_path/assembly_1st") {
    
    `mkdir -p $out_path/assembly_1st`;
    `mkdir -p $out_path/log`
    
} else {
    print STDERR "Out path exists, will delete fasta files if any...\n";
    `find $out_path/assembly_1st/ -name "*.fa" | xargs -P $num_threads rm -rf`;
}

##delete any tmp files;

`find $in_path/ -name "*cap*" | xargs -P $num_threads rm -rf`;

$assembler = ($assembler eq 'velvet') ? 'velvet' : 'cap3';
$kmer      = $kmer ? $kmer : 27;
unless($cmd) {
    #Parameters for short reads.
    if ($assembler eq 'cap3') {
        $cmd = '-r 0 -i 30 -j 31 -o 18 -s 300 -p 85';
    } else {
        $cmd = '-i 30 -j 31 -o 18 -s 300 -p 85';
    }
}

`echo "cap3 $cmd" > $out_path/log/assembly.1par`;


## if from storable data structure.
if ($dat) {
    $out_seqs = retrieve("$in_path/exported_reads.dat");
}

## if provided sub file list...   
if ($sub_lst) {
    open(my $in_fh, "$sub_lst") or die "$!";
    while(<$in_fh>) {
        next if /^$/;
        chomp;
        $_ =~ s/\.fa|\.fasta//g if $dat;
        push @files, $_;
    }
    close $in_fh;
} elsif ($dat) {
    @files = keys %{$out_seqs};
    open(my $out_fh, ">$out_path/log/reads_1st.log");
    map {print $out_fh $_.'.fa',"\n";} @files;
    close $out_fh;
} else {
    get_flist($in_path,$out_path,"reads_1st","fa");
}
die "No fasta files in $in_path.\n" if scalar(@files) == 0;

## parse tags.
my (%seqs, $v);
parse_tags();

## multi-threading part for the 1st assembly;
my $thread_m = new Parallel::ForkManager($num_threads);
unless ($sec_asmb) {
    
    $time = localtime;
    print "[$time] Starting the first run...\n";
    my $total = @files;
    my $i;
    foreach $file (@files){
    
        $i++;
        $thread_m->start and next;
        if ($dat) {
            my $seq  = $out_seqs->{$file};
            my $path = "$in_path/". $file . '.fa';
            open(my $out, ">$path");
            print $out $seq;
            close $out;
            $file .= '.fa';
        }
        if ($assembler eq 'cap3') {
            run_assembly($file, $in_path, "$out_path/assembly_1st", 1, $cmd);
        } else {
            run_hybrid_assembly($file, $in_path, "$out_path/assembly_1st", 1, $cmd);
        }
        print STDERR "  Assembling $i of $total \r";
        $thread_m->finish;
        
    }
    $thread_m->wait_all_children;
    `find $in_path -name "*.fa"|xargs -P $num_threads rm -rf` if $dat;
}

## 2nd assembly;
$time = localtime;
print "\n[$time] Starting the second run...\n";
`find $out_path/assembly_1st/ -name "*cap*" | xargs -P $num_threads rm -rf` if $sec_asmb;
`echo "cap3 $cmd" > $out_path/log/assembly.2par`;

if (! -e "$out_path/assembly_2nd") {
    
    `mkdir -p $out_path/assembly_2nd`;
    
} else {
    print STDERR "Out path exists, will delete fasta files if any...\n";
    `find $out_path/assembly_2nd/ -name "*.fa" | xargs -P $num_threads rm -rf`;
}

if ($dat && $sec_asmb) {
    # continue 2nd assembly with storable dat.
    $out_seqs = retrieve("$out_path/assembly_1st/reads_1st.dat");
    @files = keys %{$out_seqs};
    open(my $out_fh, ">$out_path/log/reads_2nd.log");
    map {print $out_fh $_.'.fa',"\n";} @files;
    close $out_fh;
} else {
    get_flist("$out_path/assembly_1st/", $out_path, "reads_2nd", "fa");
}
die "No fasta files in $out_path/assembly_1st/.\n" if scalar(@files) == 0;

$thread_m = new Parallel::ForkManager($num_threads*2); # the second run could be 2X faster.
my $total = @files;
my $i     = 0;
foreach $file (@files){
    
    $i++;
    $thread_m->start and next;
    if ($dat && $sec_asmb) {
        my $seq  = $out_seqs->{$file};
        my $path = "$out_path/assembly_1st/" . $file . '.fa';
        open(my $out, ">$path");
        print $out $seq;
        close $out;
        $file .= '.fa';
    }
    
    run_assembly($file, "$out_path/assembly_1st", "$out_path/assembly_2nd", 2, $cmd);
    
    print STDERR "  Assembling $i of $total \r";
    $thread_m->finish;
    
}

$thread_m->wait_all_children;

## store fasta files of the 1st assembly.
$time = localtime;
print "\n[$time] Storing the 1st assembled results...";
store_dat("$out_path/assembly_1st/", "$out_path/assembly_1st/", \@files, 'reads_1st.dat') if ($dat && !$sec_asmb);
print "done.\n";

## get the list of second assembly files;
get_flist("$out_path/assembly_2nd", $out_path, "2nd_assembled","fa");
die "No fasta files in $out_path/assembly_2nd/.\n" if scalar(@files) == 0;

## store fasta files of the 2nd assembly.
$time = localtime;
print "\n[$time] Storing the 2nd assembled results...";
store_dat("$out_path/assembly_2nd/", "$out_path/assembly_2nd/", \@files, 'reads_2nd.dat') if $dat;
print "done.\n";

## collect the final results, and delet useless files;
$time = localtime;
print "\n[$time] Starting to collect the final contigs...";
if (1) {
    my $in_path = $out_path . "/" . "assembly_2nd";
    my $name_fa = "collected_final.fa";
    if ($dat) {
        my $in_fa = "$in_path/reads_2nd.dat";
        retrieve_dat($in_fa, $out_path, $name_fa, 0);
    } else {
        parse_fasta(\@files, $in_path, $out_path, $name_fa);
    }
    print "done.\n";
}
## if collect contigs for reads2;
if ($collect) {
    
    my $in_path = $out_path . "/" . "assembly_1st";
    my $name_fa = "collected_reads1.fa";
    if ($dat) {
        my $in_fa = "$in_path/reads_1st.dat";
        retrieve_dat($in_fa, $out_path, $name_fa, 1);
    } else {
        get_flist($in_path, $out_path,"collect_reads1","fa");
    
        parse_fasta(\@files, $in_path, $out_path, $name_fa, 1);
    }
    
    `sed -i /Consensus/d  $out_path/$name_fa`;
    
}


sub parse_fasta {

    my ($files, $in_path, $out_path, $name_fa, $flag) =@_;
    my ($in_fh, $out_fh, $file);
    
    open ($out_fh, ">$out_path/$name_fa") or die "$!";
    
    foreach $file (@$files) {
        open ($in_fh, "<$in_path/$file") or die "$!";
        while (<$in_fh>) {
            if (defined $flag && $_ =~ /Consensus/) {
                # if collect contigs only for read2.
                chomp;
            }
            my $locus = substr($file, 0, -3);
            if (/^>/) {
                print $out_fh "\n" if $. > 1; # a new start of fasta.
                $_ =~ s/^>/>$locus\./ if $_ !~ /Consensus/; # if not consensus.
            } else {
                $_ =~ s/\n$//; # seq to one line.
            }
            print $out_fh $_;
        }
        print $out_fh "\n";
        close $in_fh;
    }
    
    close $out_fh;
    
}

sub parse_tags {

    # extract consensus from stacks catalogs;
    if (substr($loci, -2) eq 'gz') {
        open($in_tags, "gunzip -c $loci|") or die "$!";
    } elsif (substr($loci, -2) eq '7z') {
        open($in_tags, "7zcat.sh $loci|") or die "$!";
    } else {
        open($in_tags, "$loci") or die "$!";
    }
    
    while (<$in_tags>) {
    
        next if /^#/;
        my ($id, $seq);
        my @parts = split(/\t/, $_);
        if (!$v) {
            # only the use the first line.
            my $npart = @parts;
            if ($npart == 14) {$v = 1;}
            elsif ($npart == 9) {$v = 2;}
            else {die "Stacks tags files error!";}
        }
   
        if ($v == 1) {
            $id  = $parts[2];
            $seq = $parts[9];
        } elsif ($v == 2) {
            $id  = $parts[1];
            $seq = $parts[5];
        } else {
            die "Stacks tags files error!";
        } 
        # die "catalog file is wrong!!\n" if ($id+1 ne $.); # this is not needed for stacks 2.0 
        $seq = uc(reverse $seq);
        $seq =~ tr/ATCGN/TAGCN/;
        $seqs{$id} = $seq;
        
    }

    close $in_tags;

}

sub store_dat {
    my ($in_path, $out_path, $files, $name) = @_;
    my $out_seqs = {};
    foreach (@$files) {
        my $file = $_;
        my $seq  = '';
        open(my $in_fh, "$in_path/$file");
        sysread($in_fh, $seq, 1e10);
        close $in_fh;
        $out_seqs->{substr($file, 0, -3)} = $seq;
    }
    store $out_seqs, "$out_path/$name";
    `find $in_path -name "*.fa"|xargs -P $num_threads rm -rf`;
}

sub retrieve_dat {
    
    my ($in_fa, $out_path, $name_fa, $flag) = @_;
    my $out_seqs = retrieve($in_fa);
    open(my $out_fh, ">$out_path/$name_fa") or die "$!";
    foreach my $locus (sort keys %{$out_seqs}) {
        my $seq = $out_seqs->{$locus};
        $seq    =~ s/>/>$locus\./g;
        $seq    =~ s/[^\d]\n[^>]//ig; # seq to one line.
        $seq    =~ s/Consensus.*[\r\n]+/Consensus/g if $flag;
        print $out_fh $seq;
    }
    close $out_fh;
}

sub run_assembly {

    
    my ($file, $in_path, $out_path, $run, $cmd) = @_;
    
    `cap3 $in_path/$file $cmd 1>/dev/null 2>>$out_path/error.log`; # Changed, maybe some bugs here.
    
    if (!-z "$in_path/$file.cap.contigs") {
        # reads are assembled.
        `mv $in_path/$file.cap.contigs $out_path/$file`;
        
        if ($run == 2) {
            my $read1 = `grep '_Consensus' $in_path/$file.cap.singlets 2>/dev/null`;
            if ($read1) {
                # read1 is not assembled, link it.
                my $locus = substr($file, 0, -3);
                my $seq   = '>'.$locus.'_Consensus'."\n";
                   $seq  .= $seqs{$locus};
                open(my $out_fh, ">>$out_path/$file") or die "$!";
                print $out_fh $seq;
                close $out_fh;
                link_fasta($out_path, $out_path, $file);
            }
        }
    } elsif ($run == 2) {
        # connect read1 and contig of read2 by 10xN.
        link_fasta($in_path, $out_path, $file);
    }
    
    if ($run == 1) {
        my $locus = substr($file, 0, -3);
        my $seq = $seqs{$locus};
        open (my $out_fh, ">>$out_path/$file");
        print $out_fh 
        ">" . $locus . "_Consensus" . "\n" . $seq . "\n";
        
        close $out_fh;
    }
        
    my $tmp .= "$in_path/$file.cap.ace $in_path/$file.cap.contigs.links ";
       $tmp .= "$in_path/$file.cap.contigs.qual $in_path/$file.cap.info ";
       $tmp .= "$in_path/$file.cap.singlets ";
       $tmp .= "$in_path/$file.cap.contigs ";
    `rm -rf $tmp`;
    
}

sub run_hybrid_assembly {
    
    my ($file, $in_path, $out_path, $run, $cmd) = @_;
    # Run Velvet assembly;
    `velveth $out_path/$file\_out $kmer -fasta -short $in_path/$file && \
    velvetg $out_path/$file\_out -cov_cutoff auto -exp_cov auto -clean yes && \
    mv $out_path/$file\_out/contigs.fa $out_path/$file`;
    
    my $locus = substr($file, 0, -3);
    my $seq = $seqs{$locus};
    open (my $out_fh, ">>$out_path/$file");
    print $out_fh 
    ">" . $locus . "_Consensus" . "\n" . $seq . "\n";
    close $out_fh;
    `rm -rf $out_path/$file\_out`;
}

sub get_flist {
    my ($in_path, $out_path, $name, $type) = @_;
    
    `mkdir -p $out_path/log/`;
    
    undef (@files);
    opendir (D, $in_path) or die "$!";
    open(my $out_list, ">$out_path/log/$name\.log") or die "$!";

    # get the file list;
    
    while (($file = readdir(D))) {
        next if $file !~ /.+\.$type$/;
        print $out_list $file . "\n";
        push (@files, $file);
    }
    close $out_list;
    close D;
}

sub link_fasta {
    my ($in_path, $out_path, $file) = @_;
    open(IN, "$in_path/$file") or die "$!";
    my ($ctg, $seq, $read1, $num);
    my $buf  = "";
    my $plen = 0;
    my $id   = "";
    while(<IN>) {
        chomp;
        $_ =~ s/\r//g;
        if (/^>/) {
            my $len = length($buf);
            if($id =~ /Consensus/) {
                $read1 = $buf; # read1 consensus.
            } elsif ($len > $plen) {
                $ctg   = $buf; # longest contig.
                $plen  = $len;
            }
            $buf = "";
            $id = $_;
            $num++;
        } else {
            $buf .= $_;
        }
    }
    if($id =~ /Consensus/) {
        $read1 = $buf; # read1 consensus.
    } elsif (length($buf) > $plen) {
        $ctg   = $buf; # longest contig.
    }
    close IN;
    
    if ($num == 1) {
        $seq = $buf;
    } else {
        $id = '>ctg_10N_rd1';
        if ($ctg =~ /NNNNNNNNNN$/) {
            $seq = $ctg . $read1;
        } else {
            $seq = $ctg . 'N'x10 . $read1;
        }
    }
    open (my $out_fh, ">$out_path/$file");
    print $out_fh $id, "\n",
    $seq, "\n";
    close $out_fh;
}