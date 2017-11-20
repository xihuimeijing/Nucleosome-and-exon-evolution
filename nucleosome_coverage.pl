#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
my ($bed,$symbol,$mid_e,%hash);
my ($inFile,$column,$nucleosomeDyad);
GetOptions(
            'b|bed=s'          => \$inFile,
            'i|ib=i'           => \$column,
	    'n|nc=s'           => \$nucleosomeDyad,
            'h|help'           => sub{usage()}
	  ) || usage();
open IN,"intersectBed -wo -a $inFile -b $nucleosomeDyad|" or die "Can't open file:$!";
my $total=`wc -l $nucleosomeDyad|cut -f1 -d ' '`;
while(<IN>){
    chomp;
    my @split=split;
    my $dyad_e=$split[$column+2];
    $mid_e=$split[2]-($split[2]-$split[1]-1)/2;
    if($split[5] ne "+" && $split[5] ne "-" || $column<=4){
	$hash{$dyad_e-$mid_e}+=1;
    }else{
	if($split[5] eq "+"){
	    $hash{$dyad_e-$mid_e}+=1;
	}else{
	    $hash{$mid_e-$dyad_e}+=1;
	}
    }
}
for my $i(sort keys %hash){
    say join "\t",$i,$hash{$i},$hash{$i}/$total;
}
sub usage{
print <<HELP; 
Usage:	perl $0 -b interest.bed6 -i 7 -n nucleosome_dyad.bed >result.file 2>log
        Output the nucleosome dyad fraction across the input region's midpoint.This script will also generate a file called closest
Author: Yumei Li, 2015-8-21
        Revised by Yumei Li in 2015-12-23, including the intersectBed process in the script.
Output: relative_position count frequency
Options:
        -b  FILE           The input region in bed6/bed6+ format or bed4 format(Without considering strand information).(Must have unique name column)
        -i  INT            The total column number of the input BED file.
        -n  FILE           The nucleosome dyad positions in bed format.
        'help|h'            Print this help message    
HELP
    exit(-1);
}