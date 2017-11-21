#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
my($bed,$bw,$strandNum,$upstream,$downstream,$equal);
my $opt=GetOptions(
                        'b|bed=s'       => \$bed,
                        'w|bw=s'        => \$bw,
                        's|sd=i'        => \$strandNum,
                        'u|up=i'        => \$upstream,
                        'd|down=i'      => \$downstream,
                        'e|eq'          => \$equal,
                        'h|help'	=> sub{&usage;exit(-1);}
                  );
my $BED;
if(defined $bed){
    open $BED,"$bed" or die "Can't open file $bed:$!";
}else{
    $BED=\*STDIN;
}
if(! defined $strandNum){
    while(<$BED>){
        chomp;
        my @split=split;
        my($chr,$start,$end)=@split[0..2];
        if(defined $equal){
            $upstream=$end-$start;
            $downstream=$end-$start;
        }
        my($exon,$up_intron,$down_intron,$up,$down);
        chomp($exon=`bigWigSummary "$bw" $chr $start $end 1`);
        my $up_s=$start-$upstream;
        my $do_e=$end+$downstream;
        chomp($up_intron=`bigWigSummary "$bw" $chr $up_s $start 1`);
        chomp($down_intron=`bigWigSummary "$bw" $chr $end $do_e 1`);
        $up=$exon-$up_intron;
        $down=$exon-$down_intron;
        say join "\t",@split,$exon,$up_intron,$down_intron,$up,$down;
    }
}else{
    while(<$BED>){
        chomp;
        my @split=split;
        my($chr,$start,$end,$strand)=@split[0..2,$strandNum-1];
        if(defined $equal){
            $upstream=$end-$start;
            $downstream=$end-$start;
        }
        my($exon,$up_intron,$down_intron,$up,$down);
        chomp($exon=`bigWigSummary "$bw" $chr $start $end 1`);
        if($strand eq "+"){
            my $up_s=$start-$upstream;
            my $do_e=$end+$downstream;
            chomp($up_intron=`bigWigSummary "$bw" $chr $up_s $start 1`);
            chomp($down_intron=`bigWigSummary "$bw" $chr $end $do_e 1`);
        }else{
            my $up_e=$end+$upstream;
            my $do_s=$start-$downstream;
            chomp($up_intron=`bigWigSummary "$bw" $chr $end $up_e 1`);
            chomp($down_intron=`bigWigSummary "$bw" $chr $do_s $start 1`);
        }
        $up=$exon-$up_intron;
        $down=$exon-$down_intron;
        say join "\t",@split,$exon,$up_intron,$down_intron,$up,$down;
    }
}

sub usage{
print STDERR <<HELP 
Usage:	perl $0 -b Human_exon_age_hg19.bed6+ -w GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bigWig
>> result.file
        Statistics the difference between nucleosome occupancy of exon and equal(fixed) length flanking intron regions.
Output: input_fields exon_NC upIntron_NC downIntron_NC exon-upIntron exon-downIntron
        'b|bed'    FILE     The input file separated by tab or blank.[can be bed4 bed6 or bed4+](If not given, it can be read from STDIN)
        'w|bw'     FILE     The nucleosome occupancy data in bigWig format.
        's|sd'     INT      The strand column number.(If undifined,ignore the strand information.)
        'u|up'     INT      The upstream length of intron region.
        'd|down'   INT      The downstream length of intron region.
        'e|eq'     LOGIC    Using equal length of flanking intron regions. 
        'help|h'            print this help message    
HELP
}