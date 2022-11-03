# POD documentation - main docs before the code

=head1 NAME

sam_disambiguate - parses sequencing read alignments and filters categories of reads

=head1 SYNOPSIS

This script takes in a sorted (by read name) SAM file of reads aligned vs many different references
in our case small non-coding RNA references. Each read 'block' is disambiguated, checking if
the various alignments on different databases are generally overlapping on the read.
Also, reads containing specific lengths of softclip are outputed on a 'potentially chimeric' set

=head1 DESCRIPTION

The input must be a SAM file sorted by read name

Each block of reads starting with the same name is first disambiguated:

For blocks with more than one read:
- the matching part of each read (M and I characters) is identified.
- the union and intersection length of matches is calculated
- if the intersection length is 80%+ of the union length, the reads are considered similar enough
- if that is not the case, the whole block is discarded (these are reads that have different parts matching)
- for reads that were similar enough, the read with the longest S part is selected as representative of the block
- the reference name of the representative read has the references of the other members of the block appended with |
    (e.g. miR-432|snorna-1)

Blocks with one read, and representatives of multi-read blocks follow here:
- Reads with less than 15 matching nts are discarded
- Reads with 15 or more S on one side, and 6 or less S on the other are sent to 'potentially chimeric'
- Other reads (i.e. more than 15M and less than 6S on each side) are printed as SINGLE reads (non-chimeric)

Potentially chimeric reads:
- the whole read is printed in 'ofile_chim' file (SAM)
- the long (>=15) Softclipped part is extracted and a fastq file made 'ofile_chim_softclip'
- the mapped (>=15) part is also extracted and a ??? file is made 'ofile_chim_mapped'

=cut

# Let the code begin...


#!/usr/bin/env perl

use warnings;
use strict;
use autodie;
use Getopt::Long::Descriptive;
use Pod::Usage;

my ($opt, $usage) = describe_options(
	"\nUsage: %c %o",
	[],
    ["sam_disambiguate - parses sequencing read alignments and filters categories of reads"],
	[],
	['ifile=s',
		'#infile must be SAM file SORTED by read name. Use - for STDIN',
		{required => 1}],
	['ofile_full=s',
		'#outfile for fully aligned reads (SAM)',
		{required => 1}],
    ['ofile_chim=s',
		'#outfile for potentially chimeric reads (SAM)',
		{required => 1}],
    ['ofile_chim_softclip=s',
		'#outfile for potentially chimeric reads (FASTQ of softclipped part)',
		{required => 1}],
    ['ofile_chim_mapped=s',
		'#outfile for potentially chimeric reads (SAM of mapped part)',
		{required => 1}],
    ['ofile_chim_refref=s',
		'#outfile for potentially chimeric reads multiple parts mapping on reference files',
		{required => 1}],
	['verbose|v', 'Print progress'],
	['help|h', 'Print usage and exit',
		{shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;


my $IN = filehandle_for($opt->ifile);
my $ofile_full = $opt->ofile_full;
my $ofile_chim = $opt->ofile_chim;
my $ofile_chim_softclip = $opt->ofile_chim_softclip;
my $ofile_chim_mapped = $opt->ofile_chim_mapped;
my $ofile_chim_refref = $opt->ofile_chim_refref;

open (my $OUTSINGLE, ">", $ofile_full) or die "cannot open $ofile_full for output\n";
open (my $OUTFULLCHIM, ">", $ofile_chim) or die "cannot open $ofile_chim for output\n";
open (my $OUTSOFTCHIM, ">", $ofile_chim_softclip) or die "cannot open $ofile_chim_softclip for output\n";
open (my $OUTMAPPEDCHIM, ">", $ofile_chim_mapped) or die "cannot open $ofile_chim_mapped for output\n";
open (my $OUTREFREFCHIM, ">", $ofile_chim_refref) or die "cannot open $ofile_chim_refref for output\n";

my $startime = time();

my @block_of_reads = ();
my $current_name = undef;
my $linecount = 0;
my $we_are_not_done_yet = 1;

while ($we_are_not_done_yet == 1) {
    my $line = <$IN>;
    if (!defined $line){ #this little loop allows the last 'block' to be resolved
        $line = "EOF";
        $we_are_not_done_yet = 0;
    }
    chomp $line;
    
    $linecount++;
    if ($linecount % 100000 == 0){
        my $current_time = time();
        if ($opt->verbose){
            warn $linecount."\t".($current_time-$startime)."sec\n";
        }
    }

    if ($line =~ /^@/){next;} #skip headers
    my $name = (split(/\t/, $line))[0];
    if (!defined($current_name) || ($name eq $current_name)){ #we are still within the same block
        push @block_of_reads, $line; #add the current line in the disambiguation block
        if (!defined($current_name)){
            $current_name = $name;
        }
    }
    else{ #we are in the start of a new block
        
        my ($outread, $is_ref_ref_chimera) = disambiguate_block(@block_of_reads);
        if (defined $is_ref_ref_chimera){
            foreach my $oread (@block_of_reads){
                print $OUTREFREFCHIM $oread."\n";
            }
        }
        @block_of_reads = (); #clean block
        push @block_of_reads, $line; #add the current line in the disambiguation block
        $current_name = $name; # set the current block name

        if (defined $is_ref_ref_chimera){
            next;
        }

        my ($toprint_flag, $read_type, $Slength, $matched_start, $matched_end) = filter_read($outread);
        if ($toprint_flag == 1){
            print $OUTSINGLE $outread."\n";
        }
        elsif ($toprint_flag == 2){
            print $OUTFULLCHIM $outread."\n";
            
            my @split_outread = split("\t", $outread);
            my $fastq_readname = $split_outread[0]."-".$split_outread[2];
                     
            if ($read_type eq 'S_start'){
                my $fastq_seq = substr($split_outread[9], 0, $Slength);
                my $qual = substr($split_outread[10], 0, $Slength);  
                print $OUTSOFTCHIM '@'.$split_outread[0]."\n$fastq_seq\n+\n$qual\n";
                my $matched_seq = substr($split_outread[9], $matched_start, ($matched_end-$matched_start));
                my $matched_qual = substr($split_outread[10], $matched_start, ($matched_end-$matched_start));
                $split_outread[9] = $matched_seq;
                $split_outread[10] = $matched_qual;
                $split_outread[5] = length($matched_seq)."M";
                print $OUTMAPPEDCHIM join("\t", @split_outread)."\n";

            }
            if ($read_type eq 'S_end'){
                my $fastq_seq = substr($split_outread[9], (-1 * $Slength));
                my $qual = substr($split_outread[10], (-1 * $Slength));
                print $OUTSOFTCHIM '@'.$split_outread[0]."\n$fastq_seq\n+\n$qual\n";
                my $matched_seq = substr($split_outread[9], $matched_start, ($matched_end-$matched_start));
                my $matched_qual = substr($split_outread[10], $matched_start, ($matched_end-$matched_start));
                $split_outread[9] = $matched_seq;
                $split_outread[10] = $matched_qual;
                $split_outread[5] = length($matched_seq)."M";
                print $OUTMAPPEDCHIM join("\t", @split_outread)."\n";
            }
        }
    }
}

close $IN;
close $OUTSINGLE;
close $OUTSOFTCHIM;
close $OUTMAPPEDCHIM;
close $OUTFULLCHIM;
close $OUTREFREFCHIM;

exit;


sub disambiguate_block{
    
    my @block_of_reads = @_;
    my $outread = undef;
    
    if ($#block_of_reads == 0){ #if there is only 1 read to disambiguate (i.e highest index is 0)
        $outread = $block_of_reads[0];
    }
    else{
        my $min_start = undef;
        my $min_stop = undef;
        my $max_start = undef;
        my $max_stop = undef;
        my $longest_S_read = $block_of_reads[0];
        my $longest_S_number = 0;
        my $outreadname = undef;
        foreach my $read (@block_of_reads){
            
            my @splitline = split(/\t/, $read);
            my $cigar = $splitline[5];
            my $cigar_expanded = '';
            if (!defined($outreadname)){
                $outreadname = $splitline[2];
            }
            else{
                $outreadname .= "|".$splitline[2];
            }
            my $reference = $splitline[2];
            my $readseq = $splitline[9];
            while($cigar =~ /(\d+)([A-Z])/g){
                my $num = $1;
                my $letter = $2;
                if ($letter eq 'D'){next;}
                if ($letter eq 'N'){next;}
                if ($letter eq 'S'){
                    if ($num > $longest_S_number){
                        $longest_S_number = $num;
                        $longest_S_read = $read;
                    }
                }
                $cigar_expanded .= $letter x $num;
            }
            $cigar_expanded =~ /([MI]+)/;
            #warn join("\t", ($reference, $cigar_expanded, $readseq, $-[0], $+[0]))."\n";

            if (!defined($min_start) || ($-[0] < $min_start)){$min_start = $-[0];}
            if (!defined($min_stop) || ($+[0] < $min_stop)){$min_stop = $+[0];}
            if (!defined($max_start) || ($-[0] > $max_start)){$max_start = $-[0];}
            if (!defined($max_stop) || ($+[0] > $max_stop)){$max_stop = $+[0];}
        }
        #warn join("\t",($min_start, $max_start, $min_stop, $max_stop)). "\n";
        my $union_length = $max_stop - $min_start;
        my $inters_length = $min_stop - $max_start;
        if ($inters_length / $union_length >= 0.8){
            my @out_sp = split("\t", $longest_S_read);
            $out_sp[2] = $outreadname;
            $outread = join("\t", @out_sp);
        }
        else{
            #less 80% overlapping set
            #warn join("\n", (@block_of_reads, '####'))
            return (undef, 1);
        }
    }
    return ($outread, undef);
}

sub filter_read{
    my ($read) = @_;
    if (!defined($read)){
        return (0, 'unmapped_read', undef, undef, undef);
    }
    
    my @splitread = split("\t", $read);
    my $cigar = $splitread[5];
    my $cigar_expanded = '';
    while($cigar =~ /(\d+)([A-Z])/g){
        my $num = $1;
        my $letter = $2;
        if ($letter eq 'D'){next;}
        if ($letter eq 'N'){next;}
        $cigar_expanded .= $letter x $num;
    }
    $cigar_expanded =~ /([MI]+)/;
    my $matched_part = $1;
    my $matched_start = $-[0];
    my $matched_end = $+[0];
    #warn join("\t", ($cigar_expanded, $matched_start, $matched_end, substr($cigar_expanded, $matched_start, ($matched_end-$matched_start))))."\n";
    if (!defined($matched_part)){
        return (0, 'unmapped_read', undef, undef, undef);
    }
    if (length($matched_part) < 15){
        return (0, 'match_under_15', undef, $matched_start, $matched_end);
    }

    my $startS = 0;
    my $endS = 0;
    #warn $cigar_expanded."\n";
    if ($cigar_expanded =~ /^(S+)/){
        $startS = length($1);
    }
    if ($cigar_expanded =~ /(S+)$/){
        $endS = length($1);
    }
    #warn $startS, "\t", $endS."\n";;
    if ($startS >= 15 && $endS <= 6){
        return (2, 'S_start', $startS, $matched_start, $matched_end);
    }
    if ($startS <= 6 && $endS >= 15){
        return (2, 'S_end', $endS, $matched_start, $matched_end);
    }

    if ($startS >= 6 && $endS >= 6){
        return (0, 'over_6_S_bothsides', undef, $matched_start, $matched_end);
    }
    
    return (1, 'normal_read', undef, $matched_start, $matched_end );
}

sub filehandle_for {
#   this opens a filehande by filename or opens STDIN if the filename is "-"
	my ($file) = @_;

	if ($file eq '-'){
		open(my $IN, "<-");
		return $IN;
	}
	else {
		open(my $IN, "<", $file);
		return $IN
	}
}
