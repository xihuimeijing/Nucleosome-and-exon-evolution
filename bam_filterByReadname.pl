#!/bin/usr/perl -w
use strict;
use 5.010;
use Getopt::Long;
my ($inFile,$readName);
GetOptions(
            'i|bam=s'          => \$inFile,
            'r|read=s'         => \$readName,
            'h|help'           => sub{usage()}
	  ) || usage();
my ($IN,%hash,$READ);
if(defined $inFile){
    open $IN,"samtools view -h $inFile|" or die "Can't open file $inFile:$!";
}else{
    $IN=\*STDIN;
}
if(defined $readName){
    open $READ,"$readName" or die "Can't open file $readName:$!";
}else{
    $READ=\*STDIN;    
}
while(<$READ>){
    chomp;
    if(! defined $hash{$_}){
	$hash{$_}="";
    }
}
while(<$IN>){
    chomp;
    if(/^@/){say $_;next;}
    my @split=split /\t/;
    if(! defined $hash{$split[0]}){
	say $_;
    }
}
sub usage{
print <<HELP;
Usage:perl $0 -i <BAM> 
Author:Yumei Li,2017-1-5
Description:This script will filter reads by raed names.
Options:
    -i|--bam   FILE    The input bam file.[Can be read from STDIN with header]
    -r|--read  FILE    The file contains read names to be filtered.[Can be read from STDIN]
    -h|--help          Print this help information. 
HELP
    exit(-1);
}