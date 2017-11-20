#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
my($bed,$length,$discard);
my $opt=GetOptions(
                        'b|bed=s'       => \$bed,
                        'l|length=i'    => \$length,
                        'd|dis=i'       => \$discard,
                        'h|help'	=> sub{&usage;exit(-1);}
                  );
my $IN;
if(defined $bed){
    open $IN,"$bed" or die "Can't open file $bed:$!";
}else{
    $IN=\*STDIN;
}
if(! defined $length){$length=100;}
if(! defined $discard){$discard=$length/2;}
while(<$IN>){
    chomp;
    my @split=split /\t/;
    my($chr,$start,$end)=@split[0..2];
    my $out_s;
    for($out_s=$start;$out_s<=$end-$length;$out_s+=$length){
        say join "\t",$chr,$out_s,$out_s+$length;
    }
    if($end-$out_s>=$discard){say join "\t",$chr,$out_s,$end;}
}
sub usage{
print STDERR <<HELP 
Usage:	perl $0 -b <*.bed> -l 100 -d 50 > output.bed
        Split the input bed regions to fixed length of intervals.
        'b|bed'    FILE     The bed format file for the region.[can be bed4 bed6 or bed4+](If not given, it can be read from STDIN)
        'l|length' INT      The length of the output intervals.[100]
        'd|dis'    INT      The mininum interval length to output(Set for the last interval of input regions).[half of interval length]
        'help|h'            print this help message    
HELP
}
