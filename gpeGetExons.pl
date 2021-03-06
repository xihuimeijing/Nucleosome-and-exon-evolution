#!/bin/usr/perl -w
use strict;
use 5.010;
use Getopt::Long;
my ($inGpe,$bin);
GetOptions(
            'i|gpe=s'          => \$inGpe,
            'b|bin'            => \$bin,
            'h|help'           => sub{usage()}
	  ) || usage();
my $IN;
if(defined $inGpe){
    open $IN,"$inGpe" or die "Can't open file $inGpe:$!";
}else{
    $IN=\*STDIN;
}
while(<$IN>){
    chomp;
    my @splits=split /\t/;
    if(defined $bin){shift @splits;}
    my($name,$chr,$strand)=@splits[0,1,2];
    my(@starts,@ends);
    @starts=split /,/,$splits[8];
    @ends=split /,/,$splits[9];   
    for(my $i=0;$i<=$#starts;$i++){
        if($strand eq "+"){
            my $tag=$i+1;
            say join "\t",$chr,$starts[$i],$ends[$i],"$name"."|"."$tag","0",$strand;
        }else{
            my $tag=$#starts-$i+1;
            say join "\t",$chr,$starts[$i],$ends[$i],"$name"."|"."$tag","0",$strand;
        }   
    }
    
}
         
sub usage{
print <<HELP;
Usage:perl $0 -i <IN.gpe> -b >exon.bed6
Author:Yumei Li,2016-9-23
Description:This script will extract exons from GPE file.
Options:
    -i|--gpe   FILE     The input gpe file.[Can be read from STDIN]
    -b|--bin   LOGIC    With bin column in GPE file.
    -h|--help           Print this help information. 
HELP
    exit(-1);
}