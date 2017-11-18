#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use List::Util qw(min max);
my ($bed,$upstream,$downstream);
my $opt=GetOptions(
                        'b|bed=s'       => \$bed,
                        'u|up=i'        => \$upstream,
                        'd|down=i'      => \$downstream,
                        'h|help'	=> sub{&usage;exit(-1);}
                  );
my ($BED,%hash);
if(defined $bed){
    open $BED,"$bed" or die "Can't open file $bed:$!";
}else{
    $BED=\*STDIN;
}
while(<$BED>){
    chomp;
    my @split=split;
    my $length=$split[4];
    $hash{'0'}+=1;
    if($length%2 == 1){ #odd number
        my $range=int($length/2);
        for(my $i=1;$i<=(min $range,$downstream);$i++){
            $hash{$i}+=1;
        }
        for(my $j=-1;$j>=(max -$range,$upstream);$j--){
            $hash{$j}+=1;
        }
    }else{ # even number
        my $range=$length/2;
        for(my $i=1;$i<=(min $range,$downstream);$i++){
            $hash{$i}+=1;
        }
        for(my $j=-1;$j>=(max 1-$range,$upstream);$j--){
            $hash{$j}+=1;
        }
    }
}
foreach my $pos(keys %hash){
    say join "\t",$pos,$hash{$pos},$hash{$pos}/$hash{'0'};
}
sub usage{
print STDERR <<HELP 
Usage:	perl $0 -b <*.bed> -u -500 -d 500 >result.file 2>log
        Statistics the coverage across th input point. 
Author: Yumei Li, 2015-8-21
Output: relative_position count frequency
        'b|bed'    FILE     The midpoint of a region in bed format,the length of the region in 5th column.(If not given, it can be read from STDIN)
        'u|up'     INT      The upstream region from midpoint, negative value
        'd|down'   INT      The downstream region from midpoint, positive value
        'help|h'            Print this help message    
HELP
}