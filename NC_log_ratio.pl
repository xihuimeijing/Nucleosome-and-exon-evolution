#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
my($bed,$bw,$strandNum,$upstream,$downstream,$equal,$plus);
my $opt=GetOptions(
                        'b|bed=s'       => \$bed,
                        'w|bw=s'        => \$bw,
                        's|sd=i'        => \$strandNum,
                        'p|plus=f'      => \$plus,
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
if(! defined $strandNum){$strandNum=6;}
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
    if(! defined $exon){
        $up="NA";
        $down="NA";
    }elsif($exon == 0){
        if(defined $plus){
            if($up_intron == 0){$up_intron+=$plus;}
            if($down_intron == 0){$down_intron+=$plus;}
            $up=log(($exon+$plus)/$up_intron)/log(2);
            $down=log(($exon+$plus)/$down_intron)/log(2);
        }else{
            $up="-Inf"; #Not reasonable when up_intron or down_intron is zero
            $down="-Inf";
        }
    }else{
        if(! defined $up_intron){
           $up="NA";
        }elsif($up_intron==0){
            if(defined $plus){
                $up=log($exon/($up_intron+$plus))/log(2);  
            }else{
                $up="Inf";
            }
        }else{
           $up=log($exon/$up_intron)/log(2);
        }
        if(! defined $down_intron){
           $down="NA";
        }elsif($down_intron==0){
            if(defined $plus){
                $down=log($exon/($down_intron+$plus))/log(2);   
            }else{
                $down="Inf";
            }
        }else{
           $down=log($exon/$down_intron)/log(2);
        }
    }
    say join "\t",@split,$up,$down;
}
sub usage{
print STDERR <<HELP 
Usage:	perl $0 -b Human_exon_age_hg19.bed6+ -w GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bigWig
>> result.file
        Statistics the base-2-log-ratio(Set to NA if there is no data avaiable) between nucleosome occupancy of exon and equal(fixed) length flanking intron regions.
Output: input_fields exon/upstream_intron exon/downstream_intron
        'b|bed'    FILE     The input file separated by tab or blank.[can be bed4 bed6 or bed4+](If not given, it can be read from STDIN)
        'w|bw'     FILE     The nucleosome occupancy data in bigWig format.
        'p|plus'   FLOAT    The pseudo value to be assigned when the NC is zero.IF not define, the result will be assigned as Inf/-Inf. 
        's|sd'     INT      The strand column number.[default:6]
        'u|up'     INT      The upstream length of intron region.
        'd|down'   INT      The downstream length of intron region.
        'e|eq'     LOGIC    Using equal length of flanking intron regions. 
        'help|h'            print this help message    
HELP
}