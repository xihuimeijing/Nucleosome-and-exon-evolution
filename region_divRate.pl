#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use Bio::DB::Fasta;
my ($inFile,$refGene,$ancestor);
GetOptions(
            'b|bed=s'          => \$inFile,
            'a|ans=s'          => \$ancestor,
            'f|fa=s'           => \$refGene,
            'h|help'           => sub{usage()}
	  ) || usage();
my ($BED,$total,$change,%hash);
if(defined $inFile){
    open $BED,"$inFile" or die "Can't open file $inFile:$!";
}else{
    $BED=\*STDIN;
}
my $db = Bio::DB::Fasta->new($refGene);
my $ancestor_db = Bio::DB::Fasta->new($ancestor);
while(<$BED>){
    chomp;
    my @split=split;
    next if($split[0] eq 'chr');
    my($chr,$start,$end)=@split[0..2];
    my $ref_seq=$db->seq($chr,$start+1,$end);
    my $anc_seq=$ancestor_db->seq($chr,$start+1,$end);
        if((!defined $ref_seq) || (!defined $anc_seq)){
	say "$_\tNA";
	next;
    }
    my @ref_seq=split //,uc($ref_seq);
    #$chr=~s/chr//; #the header in ancestor fa has no 'chr'
    my @anc_seq=split //,uc($anc_seq);
    $hash{'A-T'}{'change'}=0;
    $hash{'C-T'}{'change'}=0;
    $hash{'G-C'}{'change'}=0;
    $hash{'G-T'}{'change'}=0;
    $hash{'T-C'}{'change'}=0;
    $hash{'T-G'}{'change'}=0;
    $hash{'A-T'}{'total'}=0;
    $hash{'C-T'}{'total'}=0;
    $hash{'G-C'}{'total'}=0;
    $hash{'G-T'}{'total'}=0;
    $hash{'T-C'}{'total'}=0;
    $hash{'T-G'}{'total'}=0;
    for(my $i=0;$i<=$#ref_seq;$i++){
	next if($anc_seq[$i]!~/A|T|G|C/);
	if($anc_seq[$i] eq 'C'){
            $hash{'C-T'}{'total'}+=1;
            $hash{'G-T'}{'total'}+=1;#C->A
            $hash{'G-C'}{'total'}+=1;#C->G
            if($ref_seq[$i] eq 'G'){
                $hash{'G-C'}{'change'}+=1;
            }elsif($ref_seq[$i] eq 'T'){
                $hash{'C-T'}{'change'}+=1;
            }elsif($ref_seq[$i] eq 'A'){
                $hash{'G-T'}{'change'}+=1;
            }
        }elsif($anc_seq[$i] eq 'G'){
	    $hash{'C-T'}{'total'}+=1;#G->A
            $hash{'G-T'}{'total'}+=1;
            $hash{'G-C'}{'total'}+=1;
            if($ref_seq[$i] eq 'C'){
                $hash{'G-C'}{'change'}+=1;
            }elsif($ref_seq[$i] eq 'T'){
                $hash{'G-T'}{'change'}+=1;
            }elsif($ref_seq[$i] eq 'A'){
                $hash{'C-T'}{'change'}+=1;
            }
        }elsif($anc_seq[$i] eq 'A'){
            $hash{'A-T'}{'total'}+=1;
            $hash{'T-C'}{'total'}+=1;#A->G
            $hash{'T-G'}{'total'}+=1;#A->C
            if($ref_seq[$i] eq 'T'){
                $hash{'A-T'}{'change'}+=1;
            }elsif($ref_seq[$i] eq 'G'){
                $hash{'T-C'}{'change'}+=1;
            }elsif($ref_seq[$i] eq 'C'){
                $hash{'T-G'}{'change'}+=1;
            }
        }elsif($anc_seq[$i] eq 'T'){
            $hash{'A-T'}{'total'}+=1;#T->A
            $hash{'T-C'}{'total'}+=1;
            $hash{'T-G'}{'total'}+=1;
            if($ref_seq[$i] eq 'A'){
                $hash{'A-T'}{'change'}+=1;
            }elsif($ref_seq[$i] eq 'G'){
                $hash{'T-G'}{'change'}+=1;
            }elsif($ref_seq[$i] eq 'C'){
                $hash{'T-C'}{'change'}+=1;
            }
        }
    }
    say join "\t",@split,$hash{'A-T'}{'total'},$hash{'C-T'}{'total'},$hash{'A-T'}{'change'},$hash{'C-T'}{'change'},$hash{'G-C'}{'change'},$hash{'G-T'}{'change'},$hash{'T-C'}{'change'},$hash{'T-G'}{'change'};
    undef %hash;
}
sub usage{
print <<HELP;
Usage: perl $0 -b exon.bed  -a ancestor_mm9.fa -f mm9.fa >result
Author: Yumei Li,2015-8-17
Revision: Recised by Yumei Li adding different substitution type on 2017-6-29.
Description: Output the divergence counts(rate) from ancestor nucleotide to derived nucleotide.
Output: Inupt_fields ancestor_A/T ancestor_G/C A->T C->T G->C G->T T->C T->G    
Options:
    -b  FILE           The input region in bed format.[If not given, it can be read from STDIN]
    -a  FILE           The ancestor sequence for this species in fa format.
    -f  FILE           Reference genome for specific species in fa format.
    -h  --help         Print this help information.
HELP
    exit(-1);
}