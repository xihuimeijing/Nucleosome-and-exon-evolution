#!/bin/usr/perl -w
use strict;
use 5.010;
use Getopt::Long;
my ($inFile,$fragment,%hash);
GetOptions(
            'i|bam=s'          => \$inFile,
	    'f|size=s'         => \$fragment,
            'h|help'           => sub{usage()}
	  ) || usage();
open IN,"samtools view $inFile|" or die "Can't open file $inFile:$!";
my($min_length,$max_length)=split /,/,$fragment;
while(<IN>){
    chomp;
    my @fields=split;
    my($readName,$chrom,$pos,$Mchrom,$Mpos,$inSize,$seq)=@fields[0,2,3,6,7,8,9];
    next if($Mchrom ne "=" || (abs $inSize)<$min_length || (abs $inSize)>$max_length);
    my $readLen=length $seq;
    if(exists $hash{$readName}){
        next;
    }else{
        my $start;
	$inSize=abs($inSize);
	if($readLen<$inSize){
	    $pos>$Mpos?$start=$Mpos:$start=$pos;
	}else{
	    $pos>$Mpos?$start=$pos:$start=$Mpos;
	}
        say join "\t",$chrom,$start+int($inSize/2)-1,$start+int($inSize/2),$inSize;
        $hash{$readName}='';    
    }
}
sub usage{
print <<HELP;
Usage:perl $0 -i file.bam -f min_size,max_size >out.bed
Author:Yumei Li,2016-3-2
       Revised by Yumei Li,2016-5-9,add comparison between read length and fragment range.
       Revised by Yumei Li,2017-10-17, add auto detection of read length.
Description:This script prints fragment midpoints(The second point if the length is even) from input bam file based on size range [min_size,max_size].The fourth column of the output file is the fragment size.
Options:
    -i|--bam   FILE    The input bam format file.(Pair-end reads)
    -f|--size  STRING  The range of fragment length that seperated by comma.
    -h|--help   Print this help information.
HELP
    exit(-1);
}