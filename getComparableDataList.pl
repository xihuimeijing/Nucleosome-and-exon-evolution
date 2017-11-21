#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
my($in1,$in2,$field1,$field2,$diff);
my $opt=GetOptions(
                        'i1|in1=s'        => \$in1,
                        'i2|in2=s'        => \$in2,
                        'f1|field1=i'     => \$field1,
                        'f2|field2=i'     => \$field2,
                        'd|diff=f'        => \$diff,
                        'h|help'	  => sub{&usage;exit(-1);}
                  );
open IN1,"$in1" or die "Can't open file $in1:$!";
open IN2,"$in2" or die "Can't open file $in2:$!";
my %tag;
while(my $line1=<IN1>){
    chomp($line1);
    my @split=split /\t/,$line1;
    my $num1=$split[$field1-1];
    seek(IN2,0,0);
    my $lineNum=0;
    while(my $line2=<IN2>){
        chomp($line2);
        $lineNum++;
        my @split=split /\t/,$line2;
        my $num2=$split[$field2-1];
        if(abs($num2-$num1)<$diff){
            if(! exists($tag{$lineNum})){
                say "$line1\t$line2";
                $tag{$lineNum}="";
                last;
            }
        }
    }
}
sub usage{
print STDERR <<HELP 
Usage:	perl $0 -i1 file1.tsv -i2 file2.tsv -f1 1 -f2 1 -d 5 >result.txt
Description: Get subset from two files without significant difference.
Author: Yumei Li,2016/11/10        
        'i1|in1'     FILE     The first input file separated by tab.
        'i2|in2'     FILE     The second input file separated by tab.
        'f1|field1'  INT      The column number of file1 to extract subset.
        'f2|filed2'  INT      The column number of file2 to extract subset.
        'd|diff'     INT      The maximum difference between file1 adn file2.
        'help|h'              Print this help message    
HELP
}