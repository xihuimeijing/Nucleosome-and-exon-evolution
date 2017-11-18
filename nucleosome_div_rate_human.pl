#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use Bio::Perl;
use Bio::DB::Fasta;
my($bed,$refGene,$ancestor,$upstream,$downstream,$remove);
my $opt=GetOptions(
                        'b|bed=s'       => \$bed,
                        'u|up=i'        => \$upstream,
                        'd|down=i'      => \$downstream,
                        'a|anc=s'       => \$ancestor,
                        'f|fa=s'        => \$refGene,
                        'r'             => \$remove,
                        'h|help'	=> sub{&usage;exit(-1);}
                  );
my ($BED,%hash,$process);
if(defined $bed){
    open $BED,"$bed" or die "Can't open file $bed:$!";
}else{
    $BED=\*STDIN;
}
my $db = Bio::DB::Fasta->new($refGene);
my $ancestor_db = Bio::DB::Fasta->new($ancestor);
#Initialize hash value
for(my $i=$upstream;$i<=$downstream;$i++){
    $hash{$i}{"ancestor"}{"A"}=0;
    $hash{$i}{"ancestor"}{"T"}=0;
    $hash{$i}{"ancestor"}{"G"}=0;
    $hash{$i}{"ancestor"}{"C"}=0;
    $hash{$i}{"change"}{"A-T"}=0;
    $hash{$i}{"change"}{"A-G"}=0;
    $hash{$i}{"change"}{"A-C"}=0;
    $hash{$i}{"change"}{"T-A"}=0;
    $hash{$i}{"change"}{"T-C"}=0;
    $hash{$i}{"change"}{"T-G"}=0;
    $hash{$i}{"change"}{"G-T"}=0;
    $hash{$i}{"change"}{"G-A"}=0;
    $hash{$i}{"change"}{"G-C"}=0;
    $hash{$i}{"change"}{"C-T"}=0;
    $hash{$i}{"change"}{"C-A"}=0;
    $hash{$i}{"change"}{"C-G"}=0;
}
while(<$BED>){
    chomp;
    $process+=1;
    my @split=split;
    my $chr=$split[0];
    my $point_s=$split[1]+int((split[2]-split[1])/2);
    my $point_e=$point_s+1;
    my $start=$point_s+$upstream+1; #1-base corrdinate
    my $end=$point_e+$downstream;
    my $ref_seq=uc($db->seq($chr,$start,$end));
    my @ref_seq=split //,$ref_seq;
    #$chr=~s/chr//;
    my $anc_seq=uc($ancestor_db->seq($chr,$start,$end));
    my @anc_seq=split //,$anc_seq;
    if(defined $remove){
        next if($ref_seq=~/N|-|\./ || $anc_seq=~/N|-|\./);
    }
    my $tag=$upstream;
    for(my $i=0;$i<=$#ref_seq;$i++){
        if($anc_seq[$i] eq 'A'){
            $hash{$tag}{"ancestor"}{"A"}+=1;
            if($ref_seq[$i] eq "T"){
                $hash{$tag}{"change"}{"A-T"}+=1;
            }elsif($ref_seq[$i] eq "G"){
                $hash{$tag}{"change"}{"A-G"}+=1;
            }elsif($ref_seq[$i] eq "C"){
                $hash{$tag}{"change"}{"A-C"}+=1;
            }
        }elsif($anc_seq[$i] eq 'T'){
            $hash{$tag}{"ancestor"}{"T"}+=1;
            if($ref_seq[$i] eq "A"){
                $hash{$tag}{"change"}{"T-A"}+=1;
            }elsif($ref_seq[$i] eq "G"){
                $hash{$tag}{"change"}{"T-G"}+=1;
            }elsif($ref_seq[$i] eq "C"){
                $hash{$tag}{"change"}{"T-C"}+=1;
            }
        }elsif($anc_seq[$i] eq 'G'){
            $hash{$tag}{"ancestor"}{"G"}+=1;
            if($ref_seq[$i] eq "T"){
                $hash{$tag}{"change"}{"G-T"}+=1;
            }elsif($ref_seq[$i] eq "A"){
                $hash{$tag}{"change"}{"G-A"}+=1;
            }elsif($ref_seq[$i] eq "C"){
                $hash{$tag}{"change"}{"G-C"}+=1;
            }
        }elsif($anc_seq[$i] eq 'C'){
            $hash{$tag}{"ancestor"}{"C"}+=1;
            if($ref_seq[$i] eq "T"){
                $hash{$tag}{"change"}{"C-T"}+=1;
            }elsif($ref_seq[$i] eq "G"){
                $hash{$tag}{"change"}{"C-G"}+=1;
            }elsif($ref_seq[$i] eq "A"){
                $hash{$tag}{"change"}{"C-A"}+=1;
            }
        }
        $tag++;
    }
    say STDERR "$process lines have beed processed...";
}

print "Position\t";
foreach my $key1(sort keys %{$hash{$upstream}{"ancestor"}}){
    print "$key1\t";
}
foreach my $key2(sort keys %{$hash{$upstream}{"change"}}){
    my $rate="$key2"."_rate";
    print "$key2\t$rate\t";
}
print "Total_rate\n";
foreach my $pos(sort keys %hash){
    my($total,$total_change);
    print "$pos\t";
    foreach my $ancestor(sort keys %{$hash{$pos}{"ancestor"}}){
        $total+=$hash{$pos}{"ancestor"}{$ancestor};
        print "$hash{$pos}{'ancestor'}{$ancestor}\t";
    }
    foreach my $ref(sort keys %{$hash{$pos}{"change"}}){
        $total_change+=$hash{$pos}{"change"}{$ref};
        my $rate=$hash{$pos}{"change"}{$ref}/$total;
        printf "%s\t%.6f\t",$hash{$pos}{'change'}{$ref},$rate;
    }
    my $all_rate=$total_change/$total;
    printf "%.6f\n",$all_rate;
}
sub usage{
print STDERR <<HELP 
Usage:	perl $0 -b point.bed -a human_ancestor.fa -f hg19.fa >result.file 2>log
        Statistics the interspecies divergence rate for the input point and upstream and downstream regions. 
Output: relative_position ancestor_base_count(4 columns) base_change_count/rate total_div_rate 
        'b|bed'    FILE     The input file in bed3 or bed3+ format.(If not given, it can be read from STDIN)
        'u|up'     INT      The upstream region, negative value
        'd|down'   INT      The downstream region, positive value
        'f|fa'     FILE     Reference genome, fasta format
        'a|anc'    FILE     The ancestor sequence, fasta format
        'r'        LOGIC    Filter out the lines as long as ref/ancestor sequence have one N/-/.
        'help|h'            Print this help message    
HELP
}
