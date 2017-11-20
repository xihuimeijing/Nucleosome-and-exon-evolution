#!/usr/bin/perl 
use 5.010;
use strict;
use warnings;
use Getopt::Long;   
use File::Basename;
use lib dirname $0;     
use pm::fqParser;

my ($input1File, $input2File, $minLength, $baseFrac, $qualCutoff, $qualFrac, $minAveQual, $tail,
    $output1File, $output2File, $help);
my ($base, $platform) = ('N', 'Sanger'); 

my $opt=GetOptions(           
                        'i1=s'            => \$input1File,
                        'i2=s'            => \$input2File,
                        'f|baseFrac=s'    => \$baseFrac,
                        'b|base=s'        => \$base,
                        'l|minLength=i'   => \$minLength,
                        'q|qualFrac=s'    => \$qualFrac,
                        'c|cutoff=i'      => \$qualCutoff,
                        'a|aveQual=s'     => \$minAveQual,
                        't|tail=s'        => \$tail,
                        'p|platform=s'    => \$platform,
                        'o1=s'            => \$output1File,
                        'o2=s'            => \$output2File,
                        'h|help'          => \$help
                        );
usage() if(defined $help);

open IN1,"$input1File" or die "Can't open $input1File: $!";
open IN2,"$input2File" or die "Can't open $input2File: $!";
open OUT1,">$output1File" or die "Can't open $output1File: $!";
open OUT2,">$output2File" or die "Can't open $output2File: $!";

my $read1Name;
my ($passedPeN, $badTailN, $totalPeN) = (0, 0, 0);

while($read1Name=<IN1>){
   chomp(my $seq1=<IN1>);
   my $plus1=<IN1>;
   chomp(my $qualSeq1=<IN1>);
   my $read2Name=<IN2>;
   chomp(my $seq2=<IN2>);
   my $plus2=<IN2>;
   chomp(my $qualSeq2=<IN2>);
   $totalPeN++;
   if( defined $minLength){
      next if(length $seq1 < $minLength || length $seq2 < $minLength);
   } 
   if(defined $baseFrac){
      next if( &fqParser::getBaseFraction($seq1, $base) > $baseFrac ||     
               &fqParser::getBaseFraction($seq2, $base) > $baseFrac
               );
   }
   if( defined $qualCutoff && defined $qualFrac){
      next if( &fqParser::getLowQualBaseFraction($qualSeq1, $platform, $qualCutoff) > $qualFrac ||
               &fqParser::getLowQualBaseFraction($qualSeq2, $platform, $qualCutoff) > $qualFrac
             );
   }
   if( defined $minAveQual){
      next if (&fqParser::countAveQual($qualSeq1, $platform) < $minAveQual ||
               &fqParser::countAveQual($qualSeq2, $platform) < $minAveQual
              );
   }
   if( defined $tail ){
      if( &fqParser::isBadTail($tail, $qualSeq1, $platform) || &fqParser::isBadTail($tail, $qualSeq2, $platform)){
            $badTailN++;
            next;
      }
   }
   print OUT1 $read1Name;   
   say OUT1 $seq1;
   print OUT1 $plus1;
   say OUT1 $qualSeq1;
   print OUT2 $read2Name;
   say OUT2 $seq2;
   print OUT2 $plus2;
   say OUT2 $qualSeq2;
   $passedPeN++;
}
say STDERR "Total pairs = $totalPeN";
say STDERR "Passed pairs = $passedPeN (". sprintf ("%.2f", ($passedPeN/$totalPeN*100)). "%)"; 
say STDERR "Bad tail pairs = $badTailN" if defined $tail;


sub usage{
   my $scriptName = basename $0;    
print <<HELP;
Usage: perl $scriptName --i1 input1.fastq --i2 input2.fastq -f 0.1 -q 0.5 -c 5 --o1 output1.fq --o2 output2.fq 2>report.log
      Statistic report is outputted to STDERR
   
   Input files
         --i1           FILE        The fastq file with the 1st read
         --i2           FILE        The fastq file with the 2nd read
         
   Sequence filter                        
      -f|--baseFrac     FLOAT       Discard read pairs with fraction of ambiguous bases > FLOAT. 0.1 is recommended
      -b|--base         CHAR        Ambiguous base[N]
      -l|--minLength    INT         Discard reads with length < INT
      
   Quality filter
      -q|--qualFrac     FLOAT       Discard read pairs with fraction of low-quality bases > FLOAT. 0.5 is recommended
      -c|--cutoff       INT         Low-quality base when its quality < INT. 5 is recommended
      -a|--minAveQual   FLOAT       Discard reads with average quality < FLOAT
      -t|--tail         INT1,INT2   Pair is discarded when any base quality of tail INT1 ones in either read < INT2
                                    E.g.: '2,5' means discarding pair when any base quality of tail 2 ones in either read is less than 5
      -p|--platform     CHAR|STR    The quality scoring platform. It can be
                                    S|Sanger (default),
                                    X|Solexa,
                                    I|Illumina1.3+,
                                    J|Illumina1.5+,
                                    L|Illumina1.8+
   Output files
         --o1           FILE        The fastq file with the filtered 1st read
         --o2           FILE        The fastq file with the filtered 2nd read
      -h|--help                     This help information
HELP
    exit(-1);
}
