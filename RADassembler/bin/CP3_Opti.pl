#!/usr/bin/env perl
#@ Yu-Long Li
use strict; 
use warnings; 
use Parallel::ForkManager;
use Getopt::Long;
use constant V => 180502;

my (@files, $file, $out_fh, $in_fh, $in_tags, $sub_lst, $sec_asmb, $assembler, $kmer);
my $in_path;
my $out_path;
my $loci;
my $num_threads = 1;
my $cmd = '';    ##Parameters parsed to CAP3 for the fist assembly.
my $collect;
my $help = 0;
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
            "assembler:s"   => \$assembler,
            "kmer:i"        => \$kmer,
            "help:1"        => \$help
            ) 
or die ("Error in command line arguments\n");
    
die "$usage\n" if !($in_path && $out_path && $loci) or ($help);

## If the outpath does not exist, create it;

if (! -e "$out_path/assembly_1st") {
    
    `mkdir -p $out_path/assembly_1st`;
    `mkdir -p $out_path/log`
    
} else {
    print STDERR "Out path exists, will delete fasta files if any...\n";
    `find $out_path/assembly_1st/ -name "*.fa" | xargs rm -rf`;
}

##delete any tmp files;

`find $in_path/ -name "*cap*" |xargs rm -rf`;

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

## extract consensus from stacks catalogs;

if (substr($loci, -2) eq 'gz') {
    open($in_tags, "gunzip -c $loci|") or die "$!";
    } 
elsif (substr($loci, -2) eq '7z') {
    open($in_tags, "7zcat.sh $loci|") or die "$!";
    }
else {
    open($in_tags, "$loci") or die "$!";
}
my (%seqs, $v);
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

## if provided sub file list...
if ($sub_lst) {
    open(my $in_fh, "$sub_lst") or die "$!";
    while(<$in_fh>) {
        next if /^$/;
        chomp;
        push @files, $_;
    }
    close $in_fh;
} else {
    get_flist($in_path,$out_path,"reads_1st","fa");
}

## multi-threading part for the 1st assembly;
my $thread_m = new Parallel::ForkManager($num_threads);
unless ($sec_asmb) {
    
    print "Starting the first run...\n";
    my $total = @files;
    my $i;
    foreach $file (@files){
    
        $i++;
        $thread_m->start and next;
        if ($assembler eq 'cap3') {
            run_assembly($file, $in_path, "$out_path/assembly_1st", 1, $cmd);
        } else {
            run_hybrid_assembly($file, $in_path, "$out_path/assembly_1st", 1, $cmd);
        }
        print STDERR "  Assembling $i of $total \r";
        $thread_m->finish;
        
    }
    $thread_m->wait_all_children;
    `find $in_path/ -name "*cap*" |xargs rm -rf`;
}

## 2nd assembly;
print "\nStarting the second run...\n";

`echo "cap3 $cmd" > $out_path/log/assembly.2par`;

if (! -e "$out_path/assembly_2nd") {
    
    `mkdir -p $out_path/assembly_2nd`;
    
} else {
    print STDERR "Out path exists, will delete fasta files if any...\n";
    `find $out_path/assembly_2nd/ -name "*.fa" | xargs rm -rf`;
}

get_flist("$out_path/assembly_1st/", $out_path, "reads_2nd", "fa");

$thread_m = new Parallel::ForkManager($num_threads*2); # the second run could be 2X faster.
my $total = @files;
my $i     = 0;
foreach $file (@files){
    
    $i++;
    $thread_m->start and next;
    
    run_assembly($file, "$out_path/assembly_1st", "$out_path/assembly_2nd", 2, $cmd);
    
    print STDERR "  Assembling $i of $total \r";
    $thread_m->finish;
    
}

$thread_m->wait_all_children;

`find $out_path/assembly_1st/ -name "*cap*" |xargs rm -rf`;

##collect the final results, and delet useless files;
print "\nStarting to collect the final contigs...\n";

get_flist("$out_path/assembly_2nd", $out_path, "2nd_assembled","fa");    # get the list of second assembly files;

parse_fasta(\@files, "$out_path/assembly_2nd", "$out_path", "collected_final.fa");

##    if collect contigs for reads2;
if ($collect) {
    
    my $in_fa = $out_path . "/" . "assembly_1st";
    my $out_fa = $out_path;
    my $name_fa = "collected_reads1.fa";
    
    get_flist("$out_path/assembly_1st",$out_path,"collect_reads1","fa");
    
    parse_fasta(\@files, $in_fa, $out_fa, $name_fa, 1);
    
    `sed -i /Consensus/d  $out_fa/$name_fa`;
    
}


sub parse_fasta {

    my ($files, $in_fa, $out_fa, $name_fa, $flag) =@_;
    my ($in_fh, $out_fh, $file);
    
    open ($out_fh, ">$out_fa/$name_fa") or die "$!";
    
    foreach $file (@$files) {
        open ($in_fh, "<$in_fa/$file") or die "$!";
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

sub run_assembly {

    
    my ($file, $in_path, $out_path, $run, $cmd) = @_;
    
    `cap3 $in_path/$file $cmd 1>/dev/null 2>>$out_path/error.log`; #Changed, maybe some bugs here.
    if (!-z "$in_path/$file.cap.contigs") {
        # reads are assembled.
        `mv $in_path/$file.cap.contigs $out_path/$file`;
        
        if ($run == 2) {
            my $read1 = `grep '_Consensus' $in_path/$file.cap.singlets`;
            if ($read1) {
                # read1 is not assembled, link it.
                my $locus = substr($file, 0, -3);
                my $seq   = '>'.$locus.'_Consensus'."\n";
                   $seq  .= $seqs{$locus};
                open(my $out_fh, ">>$out_path/$file") or die "$!";
                print $out_fh $seq;
                close $out_fh;
                #`echo "$seq" >> $out_path/$file`;
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
        
    my $cm = "rm -rf $in_path/$file.cap.ace $in_path/$file.cap.links $in_path/$file.cap.qual $in_path/$file.cap.info $in_path/$file.cap.singlets";
    `$cm`;
    
}

sub run_hybrid_assembly {
    
    my ($file, $in_path, $out_path, $run, $cmd) = @_;
    ##Run Velvet assembly;
    `velveth $out_path/$file\_out $kmer -fasta -short $in_path/$file`;
    `velvetg $out_path/$file\_out -cov_cutoff auto -exp_cov auto -clean yes`;
    `mv $out_path/$file\_out/contigs.fa $out_path/$file`;
    
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
    open(my $out_list, ">$out_path/log/$name\.txt") or die "$!";

    ##    get the file list;
    
    while (($file = readdir(D))) {
        next if $file !~ /.+\.$type$/;
        print $out_list $file . "\n";
        push (@files, $file);
    }
    
    return @files;
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