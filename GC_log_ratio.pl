#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use Bio::DB::Fasta;
my($bed,$refGene,$equal,$seperate);
my $opt=GetOptions(
                        'b|bed=s'       => \$bed,
                        'f|fa=s'        => \$refGene,
                        'e|eq'          => \$equal,
                        's|sep'         => \$seperate,
                        'h|help'	=> sub{&usage;exit(-1);}
                  );
my $BED;
if(defined $bed){
    open $BED,"$bed" or die "Can't open file $bed:$!";
}else{
    $BED=\*STDIN;
}
my $db = Bio::DB::Fasta->new($refGene);
if(defined $equal){
    while(<$BED>){
        chomp;
        my @split=split;
        my($chr,$start,$end)=@split[0..2];
        my $length=$end-$start;
        #numerator
        my @up=&GC_count($chr,$start+1,$end);
        my $up;
        if($up[1]==0){
           $up="NA";     
        }else{
            $up=$up[0]/$up[1];
        }
        my @down1=&GC_count($chr,$start-$length+1,$start);#upstream
        my @down2=&GC_count($chr,$end+1,$end+$length);#downstream
        if(defined $seperate){
            my($upstream,$downstream);
            if($down1[1]==0){
                $upstream="NA";
            }else{
                $upstream=$down1[0]/$down1[1];
            }
            if($down2[1]==0){
                $downstream="NA";
            }else{
                $downstream=$down2[0]/$down2[1];
            }
            my($ratio_up,$ratio_down);
            if($up eq "NA" || $upstream eq "NA" || $upstream==0 || $up==0){
                $ratio_up="NA";
            }else{
                $ratio_up=log($up/$upstream)/log(2);
            }
            if($up eq "NA" || $downstream eq "NA" || $up==0 || $downstream==0){
                $ratio_down="NA";
            }else{
                $ratio_down=log($up/$downstream)/log(2);
            }
            say join "\t",@split,$up,$upstream,$downstream,$ratio_up,$ratio_down;
        }else{
            my $down;
            if($down1[1]+$down2[1]==0){
               $down="NA"; 
            }else{
                $down=($down1[0]+$down2[0])/($down1[1]+$down2[1]);
            }
            my $ratio;
            if($up eq "NA" || $down eq "NA" || $up==0 || $down==0){
                $ratio="NA";
            }else{
                $ratio=log($up/$down)/log(2);
            }
            say join "\t",@split,$up,$down,$ratio;
        }
    }
}else{
    while(<$BED>){
        chomp;
        my @split=split;
        my @up=&GC_count($split[0],$split[1]+1,$split[2]);
        my $up;
        if($up[1]==0){
           $up="NA";     
        }else{
            $up=$up[0]/$up[1];
        }
        my @down1=&GC_count($split[0],$split[-4]+1,$split[-3]);#upstreaml
        my @down2=&GC_count($split[0],$split[-2]+1,$split[-1]);#downstream
        if(defined $seperate){
           my($upstream,$downstream);
            if($down1[1]==0){
                $upstream="NA";
            }else{
                $upstream=$down1[0]/$down1[1];
            }
            if($down2[1]==0){
                $downstream="NA";
            }else{
                $downstream=$down2[0]/$down2[1];
            }
            my($ratio_up,$ratio_down);
            if($up eq "NA" || $upstream eq "NA" || $upstream==0 || $up==0){
                $ratio_up="NA";
            }else{
                $ratio_up=log($up/$upstream)/log(2);
            }
            if($up eq "NA" || $downstream eq "NA" || $up==0 || $downstream==0){
                $ratio_down="NA";
            }else{
                $ratio_down=log($up/$downstream)/log(2);
            }
            say join "\t",@split,$up,$upstream,$downstream,$ratio_up,$ratio_down;
        }else{
            my $down;
            if($down1[1]+$down2[1]==0){
               $down="NA"; 
            }else{
                $down=($down1[0]+$down2[0])/($down1[1]+$down2[1]);
            }
            my $ratio;
            if($up eq "NA" || $down eq "NA" || $up==0 || $down==0){
                $ratio="NA";
            }else{
                $ratio=log($up/$down)/log(2);
            }
            say join "\t",@split,$up,$down,$ratio;
        }
    }
}
sub GC_count{
    my($chr,$start,$end)=@_;
    my $seq=$db->seq($chr,$start,$end);
    $seq=uc $seq;
    my @base=split //,$seq;
    my $total=0;
    my $gc=0;
    for(my $i=0;$i<=$#base;$i++){
        if($base[$i]=~/A|T|G|C/){
            $total+=1;
        }
        if($base[$i]=~/G|C/){
            $gc+=1;
        }
    }
    my @gc_ratio=($gc,$total);
    return @gc_ratio;
}
sub usage{
print STDERR <<HELP 
Usage:	perl $0 -b hg19_region -f hg19.fa >> result.file
        Statistics the base-2-log-ratio between between GC content of exon and flanking intron regions(The last four columns are two intron regions).
        Notice: The value will be set to 'NA' if GC content of the exon or the flanking intron region is 0.
Output: input_fields exonGC upGC downGC ratio(upstream_ratio downstraem_ratio)
        'b|bed'    FILE     The input file in bed6 or bed6+ format.(If not given, it can be read from STDIN)
        'f|fa'     FILE     Reference genome, fasta format
        's|sep'    LOGIC    Compute GC log ratio for the upstream and downstream intron regions seperately.
        'e|eq'     LOGIC    Use length of flanking intron the same as exon length.     
        'help|h'            Print this help message    
HELP
}
