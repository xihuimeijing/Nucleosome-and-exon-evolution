#!/usr/bin/perl
use strict;
use 5.010;
use Getopt::Long;
###All the lines about output fasta file of query speceis in the script are discarded. Because of maf is target species centered and there exists indel in query species.
my($mafFile,$output,$target,$query1,$query_out,$suffix);
GetOptions(
            'm|maf=s'       => \$mafFile,
            't|target=s'    => \$target,
            'q|query=s'     => \$query1,
            'o|out=s'       => \$output,
            '--om=s'        => \$query_out,
            '--suffix=s'    => \$suffix,
            'h|help'        => sub{usage()}
        ) || usage();
my($MAF,$targetLine,$seqNum,%hash,$OUT,%chrSize,%query_chrSize,%query_hash);
if(defined $mafFile){
    open $MAF,"$mafFile" or die "Can't open file $mafFile:$!";
}else{
    $MAF=\*STDIN;;
}
if(defined $output){
    open $OUT,">$output" or die "Can't open file $output:$!";
}else{
    $OUT=\*STDOUT;
}
my @queries=split /,/,$query1;
#variable store information
my($targetChr, $targetStart, $targetLength, $targetStrand, $targetChrS, $targetSeq, $targetEnd, @query_seqs, %query_infors);
while(<$MAF>){
    chomp;
    next if /^#/;
    next if /^[aieq] /;
    if (/^$/){
        $targetLine = "";
        #Store ancestor information of last block
        if($seqNum == $#queries+1){
            my @targetSeq=split //,uc($targetSeq);
            #defined a 2-dim array
            my @q_Seqs;
            for(my $i=0;$i<=$#queries;$i++){
                $q_Seqs[$i]=[split //,uc($query_seqs[$i])];
            }
            my $ancestor;
            for(my $i=0;$i<=$#targetSeq;$i++){
                next if($targetSeq[$i]!~/A|T|G|C|N/);
                my $finalSeq=$targetSeq[$i];
                for(my $j=0;$j<=$#query_seqs;$j++){
                    if($q_Seqs[$j][$i]=~/A|T|G|C|N/){
                        $finalSeq.=$q_Seqs[$j][$i];
                    }
                }
                if((length $finalSeq)<($seqNum+1)){
                    $ancestor.="N";
                    next;
                }
                my $A_count=($finalSeq=~tr/A/A/);
                my $T_count=($finalSeq=~tr/T/T/);
                my $G_count=($finalSeq=~tr/G/G/);
                my $C_count=($finalSeq=~tr/C/C/);
                my $support=$seqNum/2+1;
                if($A_count>=$support){
                    $ancestor.="A";
                }elsif($T_count>=$support){
                    $ancestor.="T";
                }elsif($G_count>=$support){
                    $ancestor.="G";
                }elsif($C_count>=$support){
                    $ancestor.="C";
                }else{
                    $ancestor.="N";
                }
            }
            $hash{$targetChr}{$targetStart}{$targetEnd}=$ancestor;
            ##Output query fasta files
            #if(defined $query_out){
            #    my @query_outs=split /,/,$query_out;
            #    for(my $i=0;$i<=$#query_outs;$i++){
            #        my @infors=split /_/,$query_infors{$query_outs[$i]};
            #        $query_hash{$query_outs[$i]}{$infors[0]}{$infors[1]}{$infors[2]}{$infors[3]}=$ancestor; 
            #    }
            #}
        }
        $seqNum=0;
        next;
    }
    if(/^s /){
        if($targetLine eq ""){
            $targetLine = $_;
            ($targetChr, $targetStart, $targetLength, $targetStrand, $targetChrS, $targetSeq)
                            = /s\s+.+\.(\S+)                #targetChr
                                \s+(\d+)                    #targetStart
                                \s+(\d+)                    #targetLength
                                \s+([+-])                   #targetStrand
                                \s+(\d+)                    #chrSize
                                \s+([ATCGatcgNn-]+)         #targetSeq
                                /x;
            $targetEnd = $targetStart + $targetLength;
            if(! defined $chrSize{$targetChr}){$chrSize{$targetChr}=$targetChrS;}
        }else{
            my ($queryName, $queryChr, $queryStart, $queryLength, $queryStrand, $queryChrSize, $querySeq)
                        = /s\s+(\S+?)               #queryName
                            \.(\S+)                 #queryChr
                            \s+(\d+)                  #queryStart
                            \s+(\d+)                  #queryLength
                            \s+([+-])                 #queryStrand
                            \s+(\d+)                  #queryChrSize 
                            \s+([ATCGatcgNn-]+)     #querySeq
                        /x;
            #$queryStart = $queryChrSize - ($queryStart+$queryLength) if $queryStrand eq '-';
            #my $queryEnd = $queryStart + $queryLength;
            for(my $i=0;$i<=$#queries;$i++){
                if($queryName eq $queries[$i]){
                    $query_seqs[$i]=$querySeq;
                    #$query_infors{$queryName}=join "_",$queryChr,$queryStrand,$queryStart,$queryEnd;
                    #if(! defined $query_chrSize{$queryName}{$queryChr}){$query_chrSize{$queryName}{$queryChr}=$queryChrSize}
                    $seqNum++;
                }
            }
        }
    }
}
#Output ancestor sequence in fasta format
my $faStart=0;
foreach my $chr(sort keys %hash){
    say $OUT ">$chr";
    $faStart=0;
    my $outSeq;
    foreach my $start(sort {$a <=> $b} keys %{$hash{$chr}}){
        foreach my $end(keys %{$hash{$chr}{$start}}){
            my $seq=$hash{$chr}{$start}{$end};
            if($start>$faStart){
                my $tmpSeq="."x($start-$faStart);
                $outSeq.=$tmpSeq;
            }
            $outSeq.=$seq;
            $faStart=$end;
        }
    }
    if($faStart<$chrSize{$chr}){$outSeq.="."x($chrSize{$chr}-$faStart);}
    &myPrint($outSeq,$OUT);
}
close($OUT);
#if(defined $query_out){
#    my @query_outs=split /,/,$query_out;
#        for(my $i=0;$i<=$#query_outs;$i++){
#            my $outname=$query_outs[$i]."_".$suffix;
#            my $OUT_q;
#            open $OUT_q,">$outname" or die "Can't open output file:$!";
#            my $faStart=0;
#            foreach my $chr(sort keys %{$query_hash{$query_outs[$i]}}){
#                say $OUT_q ">$chr";
#                $faStart=0;
#                my $outSeq;
#                my @strand=keys %{$query_hash{$query_outs[$i]}{$chr}};
#                foreach my $start(sort {$a <=> $b} keys %{$query_hash{$query_outs[$i]}{$chr}{$strand[0]}}){
#                    foreach my $end(keys %{$query_hash{$query_outs[$i]}{$chr}{$strand[0]}{$start}}){
#                        my $seq;
#                        if($strand[0] eq "+"){
#                            $seq=$query_hash{$query_outs[$i]}{$chr}{$strand[0]}{$start}{$end};
#                        }else{
#                            $seq=reverse $query_hash{$query_outs[$i]}{$chr}{$strand[0]}{$start}{$end};
#                            $seq=~tr/ACGT/TGCA/;
#                        }
#                        if($start>$faStart){
#                            my $tmpSeq="."x($start-$faStart);
#                            $outSeq.=$tmpSeq;
#                        }
#                        $outSeq.=$seq;
#                        $faStart=$end;  
#                    }
#                }
#                if($faStart<$query_chrSize{$query_outs[$i]}{$chr}){$outSeq.="."x($query_chrSize{$query_outs[$i]}{$chr}-$faStart);}
#                &myPrint($outSeq,$OUT_q);
#            }
#            close($OUT_q);
#        }
#}
sub myPrint{
    my($a,$b)=@_;
    my $length=length $a;
    for(my $i=0;$i<=$length;$i+=50){
        say $b substr($a,$i,50);
    }
}
sub usage{
print <<HELP;
Usage: perl $0 -m <in.maf> -t hg19 -q rheMac2,calJac1 -o out.fa 
    The script will output the most common ancestor of target species (inferred by parsimony) and query species.The total species number must be odd 'i', 'e' and 'q' lines are skipped.
Author: Yumei Li, 2017/7/25
Revision: Adding processing more than three species at 2017/8/3.
    -m|--maf            The input maf file.[Can be read from STDIN]
    -t|--target         The target species name in maf file.
    -q|--query          The query species name in maf file separated by comma.
    -o|--out            The output fasta file.
    -h|--help           Print this help information.
HELP
    exit(-1);
}