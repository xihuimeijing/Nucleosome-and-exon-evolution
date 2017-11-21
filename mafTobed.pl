#!/bin/env perl
use strict;
use 5.010;
use Getopt::Long;
my($maf,$suffix);
my $opt=GetOptions(
                'm|maf=s'       => \$maf,
                'h|help'	=> sub{&usage;exit(-1);}
                  );
my $IN;
if(defined $maf){
    open $IN,"$maf" or die "Can't open file $maf:$!";
}else{
    $IN=\*STDIN;
}

my $blockNum=1;
my $tag=0;
while(<$IN>){
    chomp;
    next if /^#/;
    next if /^[aieq] /;
    if(/^$/){print "\n";$blockNum+=1;}
    if(/^s /){
          my($name, $Chr, $Start, $Length, $Strand, $SrcSize)
                            = /s\s+(\S+?)                    #genome
                                \.(\S+)                     #Chr
                                \s+(\d+)                    #Start
                                \s+(\d+)                    #Length
                                \s+([+-])                   #tStrand
                                \s+(\d+)                    #SrcSize
                                \s+[ATCGatcgNn-]+     #targetSeq
                              /x;
        my($outS,$outE);
        if($Strand eq '+'){
            $outS=$Start;
            $outE=$Start+$Length;
        }else{
            $outS=$SrcSize-($Start+$Length);
            $outE=$outS+$Length;
        }
        if($tag == $blockNum){
            print "\t";
            print join "\t",$Chr,$outS,$outE,"$name"."_"."$blockNum",$Strand;
        }else{
            print join "\t",$Chr,$outS,$outE,"$name"."_"."$blockNum",$Strand;
            $tag=$blockNum;
        }
    }
}
sub usage{
print <<HELP;
Usage: perl $0 -m in.maf >out.bed 
Author: Yumei Li, 2017-4-7
This script output block regions of each species in maf into the same out line.
Options:
    -m|--maf       The maf file to be prosecced.[Can be read from STDIN.]
    -h --help      Print this help information.
HELP
    exit(-1);
}