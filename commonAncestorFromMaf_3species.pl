#!/usr/bin/perl
use strict;
use 5.010;
use Getopt::Long;
my($mafFile,$target,$query,$outgroup,$output1,$output2);
GetOptions(
            'm|maf=s'       => \$mafFile,
            't|target=s'    => \$target,
            'q|query=s'     => \$query,
            '--og=s'        => \$outgroup,
            '--o1=s'        => \$output1,
            '--o2=s'        => \$output2,
            'h|help'        => sub{usage()}
        ) || usage();
my($MAF,$targetLine,%hash,$OUT1,$OUT2,%chrSize,%query_chrSize,%query_hash);
if(defined $mafFile){
    open $MAF,"$mafFile" or die "Can't open file $mafFile:$!";
}else{
    $MAF=\*STDIN;;
}
open $OUT1,">$output1" or die "Can't open file $output1:$!";
open $OUT2,">$output2" or die "Can't open file $output2:$!";
#variable store information
my($targetChr, $targetStart, $targetLength, $targetStrand, $targetChrS, $targetSeq, $targetEnd, $queryChr, $queryStart, $queryStrand, $queryEnd, $querySeq, $outgroupSeq);
while(<$MAF>){
    chomp;
    next if /^#/;
    next if /^[aieq] /;
    if (/^$/){
        $targetLine = "";
        #Store ancestor information of last block
        if((defined $querySeq) && (defined $outgroupSeq)){
            my @targetSeq=split //,uc($targetSeq);
            #defined a 2-dim array
            my @querySeqs=split //,$querySeq;
            my @outgroupSeqs=split //,$outgroupSeq;
            my $ancestor;
            my $ancestor_q;
            for(my $i=0;$i<=$#targetSeq;$i++){
                next if($targetSeq[$i]!~/A|T|G|C|N/);
                my $finalSeq=$targetSeq[$i];
                $finalSeq.=$querySeqs[$i];
                $finalSeq.=$outgroupSeqs[$i];
                if($finalSeq=~/N|-/){
                    $ancestor.="N";
                    next;
                }
                my $A_count=($finalSeq=~tr/A/A/);
                my $T_count=($finalSeq=~tr/T/T/);
                my $G_count=($finalSeq=~tr/G/G/);
                my $C_count=($finalSeq=~tr/C/C/);
                if($A_count>=2){
                    $ancestor.="A";
                }elsif($T_count>=2){
                    $ancestor.="T";
                }elsif($G_count>=2){
                    $ancestor.="G";
                }elsif($C_count>=2){
                    $ancestor.="C";
                }else{
                    $ancestor.="N";
                }
            }
            $hash{$targetChr}{$targetStart}{$targetEnd}=$ancestor;
            for(my $i=0;$i<=$#querySeqs;$i++){
                next if($querySeqs[$i]!~/A|T|G|C|N/);
                my $finalSeq=$targetSeq[$i];
                $finalSeq.=$querySeqs[$i];
                $finalSeq.=$outgroupSeqs[$i];
                if($finalSeq=~/N|-/){
                    $ancestor_q.="N";
                    next;
                }
                my $A_count=($finalSeq=~tr/A/A/);
                my $T_count=($finalSeq=~tr/T/T/);
                my $G_count=($finalSeq=~tr/G/G/);
                my $C_count=($finalSeq=~tr/C/C/);
                if($A_count>=2){
                    $ancestor_q.="A";
                }elsif($T_count>=2){
                    $ancestor_q.="T";
                }elsif($G_count>=2){
                    $ancestor_q.="G";
                }elsif($C_count>=2){
                    $ancestor_q.="C";
                }else{
                    $ancestor_q.="N";
                }
            }
            $query_hash{$queryChr}{$queryStart}{$queryEnd}{$queryStrand}=$ancestor_q;
        }
        undef($querySeq);
        undef($outgroupSeq);
        undef($targetSeq);
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
            my ($queryName, $queryChr_tmp, $queryStart_tmp, $queryLength, $queryStrand_tmp, $queryChrSize, $tmpSeq)
                        = /s\s+(\S+?)                 #queryName
                            \.(\S+)                   #queryChr
                            \s+(\d+)                  #queryStart
                            \s+(\d+)                  #queryLength
                            \s+([+-])                 #queryStrand
                            \s+(\d+)                  #queryChrSize 
                            \s+([ATCGatcgNn-]+)     #querySeq
                        /x;
                if($queryName eq $query){
                    $querySeq=uc($tmpSeq);
                    if($queryStrand_tmp eq "-"){
                        $queryStart = $queryChrSize - ($queryStart_tmp+$queryLength);
                    }else{
                        $queryStart = $queryStart_tmp;
                    }
                    $queryEnd = $queryStart + $queryLength;
                    ($queryChr,$queryStrand)=($queryChr_tmp,$queryStrand_tmp);
                    if(! defined $query_chrSize{$queryChr}){$query_chrSize{$queryChr}=$queryChrSize;}
                }
                if($queryName eq $outgroup){
                    $outgroupSeq=uc($tmpSeq);
                }
        }
    }
}
#Output ancestor sequence in fasta format
my $faStart=0;
foreach my $chr(sort keys %hash){
    say $OUT1 ">$chr";
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
    &myPrint($outSeq,$OUT1);
}
close($OUT1);
my $faStart_q=0;
foreach my $chr(sort keys %query_hash){
    say $OUT2 ">$chr";
    $faStart_q=0;
    my $outSeq;
    my $formerSeq;
    foreach my $start(sort {$a <=> $b} keys %{$query_hash{$chr}}){
        foreach my $end(sort {$a <=> $b} keys %{$query_hash{$chr}{$start}}){
            my @strand=keys %{$query_hash{$chr}{$start}{$end}};
            my $seq;
            if($strand[0] eq "+"){
                $seq=$query_hash{$chr}{$start}{$end}{$strand[0]};
            }else{
                $seq=reverse $query_hash{$chr}{$start}{$end}{$strand[0]};
                $seq=~tr/ACGT/TGCA/;
            }
            #There may be duplications in target or query species.
            if(defined $formerSeq){
                my($former_start,$former_end,$former_seq)=split /;/,$formerSeq;
                if($former_start>$faStart_q){
                    my $tmpSeq="."x($former_start-$faStart_q);
                    $outSeq.=$tmpSeq;
                }
                my @seq1=split //,$former_seq;
                my @seq2=split //,$seq;
                if($start==$former_start){ #The same start
                    my $nowSeq;
                    for(my $i=0;$i<=$#seq1;$i++){
                        if($seq1[$i] eq $seq2[$i]){
                            $nowSeq.=$seq1[$i];
                        }else{
                            $nowSeq.="N";
                        }
                    }
                    if($#seq2>$#seq1){$nowSeq.=join "",@seq2[($#seq1+1)..$#seq2];}
                    $faStart_q=$start;
                    #End may be not greater than former end
                    if($end>$former_end){
                        $formerSeq=join ";",$start,$end,$nowSeq;
                    }else{
                        $formerSeq=join ";",$start,$former_end,$nowSeq;
                    }
                }elsif($start>$former_start && $start<$former_end){ #partial overlap
                    $outSeq.=join "",@seq1[0..($start-$former_start-1)];
                    my $nowSeq;
                    if($end<$former_end){
                        #Seq2 is shorter
                        my $i;
                        for($i=$start-$former_start;$i<=$start-$former_start+$#seq2;$i++){
                            if($seq1[$i] eq $seq2[$i-($start-$former_start)]){
                                $nowSeq.=$seq1[$i];
                            }else{
                                $nowSeq.="N";
                            }      
                        }
                        $nowSeq.=join "",@seq1[$i..$#seq1];
                        $formerSeq=join ";",$start,$former_end,$nowSeq;
                    }elsif($end=$former_end){
                        for(my $i=$start-$former_start;$i<=$#seq1;$i++){
                            if($seq1[$i] eq $seq2[$i-($start-$former_start)]){
                                $nowSeq.=$seq1[$i];
                            }else{
                                $nowSeq.="N";
                            }       
                        }
                        $formerSeq=join ";",$start,$former_end,$nowSeq;
                    }else{
                        my $i;
                        for($i=$start-$former_start;$i<=$#seq1;$i++){
                            if($seq1[$i] eq $seq2[$i-($start-$former_start)]){
                                $nowSeq.=$seq1[$i];
                            }else{
                                $nowSeq.="N";
                            }
                        }
                        $nowSeq.=join "",@seq2[($i-($start-$former_start))..$#seq2];
                        $formerSeq=join ";",$start,$end,$nowSeq;
                    }
                    $faStart_q=$start;
                }elsif($start>=$former_end){ #No overlap
                    $outSeq.=$former_seq;
                    $faStart_q=$former_end;
                    $formerSeq=join ";",$start,$end,$seq;
                }
            }else{ #The first line
                $formerSeq=join ";",$start,$end,$seq;
            }
        }
    }
    my($former_start,$former_end,$former_seq)=split /;/,$formerSeq;
    #Output the last line
    if($former_start>$faStart_q){
        my $tmpSeq="."x($former_start-$faStart_q);
        $outSeq.=$tmpSeq;
    }
    $faStart_q=$former_end;
    $outSeq.=$former_seq;
    if($faStart_q<$query_chrSize{$chr}){$outSeq.="."x($query_chrSize{$chr}-$faStart_q);}
    &myPrint($outSeq,$OUT2);
    undef($outSeq);
    undef($formerSeq);
}
close($OUT2);
sub myPrint{
    my($a,$b)=@_;
    my $length=length $a;
    for(my $i=0;$i<=$length;$i+=50){
        say $b substr($a,$i,50);
    }
}
sub usage{
print <<HELP;
Usage: perl $0 -m <in.maf> -t hg19 -q rheMac2,calJac1 --o1 out.hg19.fa --o2 out.rheMac2.fa
    The script will output the most common ancestor of target species (inferred by parsimony) and query species with one outgroup species. 'i', 'e' and 'q' lines are skipped.
Author: Yumei Li, 2017/9/17,
        Final revision, 2017/10/09
    -m|--maf            The input maf file.[Can be read from STDIN]
    -t|--target         The target species name in maf file.
    -q|--query          The query species name in maf file.
    --og                The outgroup species name in maf file;
    --o1                The output for target species in fasta file.
    --o2                The output for query species in fasta file.
    -h|--help           Print this help information.
HELP
    exit(-1);
}