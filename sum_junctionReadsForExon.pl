#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
my($exon,$junc);
my $opt=GetOptions(
                        'e|exon=s'     => \$exon,
                        'j|junc=s'     => \$junc,
                        'h|help'       => sub{&usage;exit(-1);}
                  );
open EXON,"$exon" or die "Can't open file $exon:$!";
open JUNC,"$junc" or die "Can't open file $junc:$!";
my %junction;
while(<JUNC>){
    chomp;
    my @split=split;
    my($chr,$corr)=split /:/,$split[3];
    my($start,$end)=split /-/,$corr;
    my $strand=$split[5];
    my $reads=$split[4];
    $junction{$chr}{$strand}{$start}=$reads;
    $junction{$chr}{$strand}{$end}=$reads;
}
while(<EXON>){
    chomp;
    my @split=split;
    my($chr,$start,$end,$name,$number,$strand)=@split[0..5];
    if(exists $junction{$chr}{$strand}{$start}){
        print "$_\t$junction{$chr}{$strand}{$start}\t";   
    }else{
        print "$_\t0\t";
    }
    if(exists $junction{$chr}{$strand}{$end}){
        say $junction{$chr}{$strand}{$end};
    }else{
        say "0";
    }
}
sub usage{
print STDERR <<HELP 
Usage:	perl $0 -e exon.bed6 -j junc.bed12 >file.rst 
Output: Output all the junction reads for the input exons.
        'e|exon'    FILE     Exon in BED6/BED6+ format.
        'j|junc'    FILE     Junction file from tophat in BED12 format.(The fifth column is junction reads number.)
         'help|h'            Print this help message    
HELP
}