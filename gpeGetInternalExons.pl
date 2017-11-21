#!/bin/usr/perl -w
use strict;
use 5.010;
use Getopt::Long;
my ($inGpe,$bin,$add,$gene);
GetOptions(
            'i|gpe=s'          => \$inGpe,
            'b|bin'            => \$bin,
	    'a|add'            => \$add,
	    'g|gene'           => \$gene,
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
    next if $#starts<=1;
    for(my $i=1;$i<$#starts;$i++){
	if(defined $add){
	    if($strand eq "+"){
		my $tag=$i+1;
		say join "\t",$chr,$starts[$i],$ends[$i],"$name"."|"."$tag","0",$strand,$ends[$i-1],$starts[$i],$ends[$i],$starts[$i+1];
	    }else{
		my $tag=$#starts-$i+1;
		say join "\t",$chr,$starts[$i],$ends[$i],"$name"."|"."$tag","0",$strand,$ends[$i],$starts[$i+1],$ends[$i-1],$starts[$i];
	    }
	}elsif(defined $gene){
	    if($strand eq "+"){
		my $tag=$i+1;
		say join "\t",$chr,$starts[$i],$ends[$i],$splits[11]."|"."$name"."|"."$tag","0",$strand;
	    }else{
		my $tag=$#starts-$i+1;
		say join "\t",$chr,$starts[$i],$ends[$i],$splits[11]."|"."$name"."|"."$tag","0",$strand;
	    }
	}else{
	    if($strand eq "+"){
		my $tag=$i+1;
		say join "\t",$chr,$starts[$i],$ends[$i],"$name"."|"."$tag","0",$strand;
	    }else{
		my $tag=$#starts-$i+1;
		say join "\t",$chr,$starts[$i],$ends[$i],"$name"."|"."$tag","0",$strand;
	    }
	}
    }
    
}
         
sub usage{
print <<HELP;
Usage:perl $0 -i <IN.gpe> -b >exon.bed6
Author:Yumei Li,2016-9-23
Description:This script will extract internal exons from GPE file.
Options:
    -i|--gpe   FILE     The input gpe file.[Can be read from STDIN]
    -b|--bin   LOGIC    With bin column in GPE file.
    -a|--add   LOGIC    Add the flanking intron of the internal exon(The last four columns are upstream intron and downstream intron coordinate).
    -g|--gene  LOGIC    Add gene name in the name field(4th column) of the output file.
    -h|--help           Print this help information. 
HELP
    exit(-1);
}