#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use Bio::DB::Fasta;
my ($inFile,$refGene,$upstream,$downstream,$chrom,$dn,$quiet);
GetOptions(
            'b|bed=s'          => \$inFile,
	    'n|dn=s'           => \$dn,
            'f|fa=s'           => \$refGene,
	    'u|up=i'           => \$upstream,
	    'd|down=i'         => \$downstream,
	    'c|chrom'          => \$chrom,
	    'q|quiet'          => \$quiet,
            'h|help'           => sub{usage()}
	  ) || usage();
my ($BED,%hash);
if(defined $inFile){
    open $BED,"$inFile" or die "Can't open file $inFile:$!";
}else{
    $BED=\*STDIN;
}
my $db = Bio::DB::Fasta->new($refGene);
my $process=0;
my @dinucleo=split /,/,$dn;
for(my $i=$upstream;$i<$downstream;$i++){
    $hash{$i}{"total"}=0;
    for(my $n=0;$n<=$#dinucleo;$n++){
	$hash{$i}{'dinucleo'}{$dinucleo[$n]}=0;
    }
}
while(<$BED>){
    chomp;
    $process+=1;
    my @split=split;
    my $chr=$split[0];
    my $start=$split[1]+$upstream+1; #1-base corrdinate
    my $end=$split[2]+$downstream+1; #The last dinucleotide
    if(! defined $chrom){$chr=~s/chr//;}
    my @ref_seq=split //,uc($db->seq($chr,$start,$end));
    my $tag=$upstream;
    for(my $i=0;$i<$#ref_seq;$i++){
	if($ref_seq[$i]=~/A|T|G|C/ && $ref_seq[$i+1]=~/A|T|G|C/){
	    $hash{$tag}{"total"}+=1;
	}
	for(my $j=0;$j<=$#dinucleo;$j++){
	    my @nucleotide=split //,$dinucleo[$j];
	    if($ref_seq[$i] eq $nucleotide[0] && $ref_seq[$i+1] eq $nucleotide[1]){
		$hash{$tag}{'dinucleo'}{$dinucleo[$j]}+=1;
	    }
	}
	$tag++;
    }
    if(!defined $quiet){
	say STDERR "$process lines have beed processed...";
    }
}
print "pos";
for my $base(sort keys %{$hash{'0'}{'dinucleo'}}){
    print "\t$base";
}
print "\n";
foreach my $key(keys %hash){
    print "$key";
    foreach my $base(sort keys %{$hash{$key}{'dinucleo'}}){
	my $freq=$hash{$key}{'dinucleo'}{$base}/$hash{$key}{'total'};
	print "\t$freq";
    }
    print "\n";
}
sub usage{
print <<HELP;
Usage: perl $0 -b nucleosome_dyad.bed -n AA -u -500 -d 500 -f hg19.fa >result
Author: Yumei Li,2015-9-6
        Revise by Yumei Li in 2015-12-13. Including the handle of mutiple dinucleotides.
Description: Output the dinucleotide frequency for the input point and its flanking region.
Options:
    -b  FILE           The input point in bed format.[If not given, it can be read from STDIN]
    -n  STRING         The dinucleotide you want to calculate frequency.
			If you have mutiple dinucleotides to sum up,you can seperate them by ','; eg "AA,AT,TT".
    -u  INT            The upstream region, negative value
    -d  INT            The downstream region, positive value
    -f  FILE           Reference genome for specific species in fa format.
    -c  LOGIC          The header of reference genome is ">chr*" format.[default:">*"]
    -q  LOGIC          Do not print processing log to STDERR.
    -h  --help         Print this help information.
HELP
    exit(-1);
}