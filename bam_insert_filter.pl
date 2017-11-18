#!/bin/usr/perl -w
use strict;
use 5.010;
use Getopt::Long;
my ($inFile,$size_range,$mapQ,$nm);
GetOptions(
            'i|bam=s'          => \$inFile,
            's|size=s'         => \$size_range,
            'q|mapQ=i'         => \$mapQ,
            'n|nm=i'           => \$nm,
            'h|help'           => sub{usage()}
	  ) || usage();
my($IN,$minS,$maxS);
if(defined $inFile){
    open $IN,"samtools view -h $inFile|" or die "Can't open file $inFile:$!";
}else{
    $IN=\*STDIN;
}
if(defined $size_range){
    ($minS,$maxS)=split /,/,$size_range;
}else{
    $minS=50;
    $maxS=300;
}
if(! defined $mapQ){$mapQ=30;}
if(! defined $nm){$nm=10;}

my ($tag,$read,$insert_first,$mapq_first,$nm_first);
while(<$IN>){
    chomp;
    if(/^@/){say $_;next;}
    my @split=split /\t/;
    my($name,$mapQ_now,$M_name_now,$insertS_now)=@split[0,4,6,8];
    $_=~/.*NM:i:(\d+)\s+.*/;
    my $nm_now=$1;
    if(! defined $tag){
        $read=$_;
        $tag=$name;
        $insert_first=$insertS_now;
        $mapq_first=$mapQ_now;
        $nm_first=$nm_now;
    }else{
        if($name eq $tag){
            if($M_name_now eq "=" && (abs($insertS_now)>=$minS && abs($insertS_now)<=$maxS) && (abs($insert_first)>=$minS && abs($insert_first)<=$maxS) && ($mapQ_now>=$mapQ && $mapq_first>=$mapQ) && ($nm_now<=$nm && $nm_first<=$nm) && ($insert_first+$insertS_now==0)){
                say $read;
                say $_;
            }
        }else{
            $read=$_;
            $tag=$name;
            $insert_first=$insertS_now;
            $mapq_first=$mapQ_now;
            $nm_first=$nm_now;
        }
    }
}
sub usage{
print <<HELP;
Usage:perl $0 -i <BAM> 
Author:Yumei Li,2016-5-8
Description:This script will print reads in proper pair based on position and insert size.It will output a SAM file.
Options:
    -i|--bam   FILE    The input bam file which should be sorted by read name and uniquely mapped reads.[Can be read from STDIN with header]
    -s|--size  STRING  The insert size range separated by comma.[default:50,300]
    -q|--mapQ  INT     The minimum mapping quality for both reads in pair.[default:30]
    -n|--nm    INT     The minimum edit distance for both reads in pair.[Require NM tag,default:10]
    -h|--help          Print this help information. 
HELP
    exit(-1);
}