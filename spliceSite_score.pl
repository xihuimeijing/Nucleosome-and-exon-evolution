#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use Bio::DB::Fasta;
my ($bedFile,$faFile,$ss3_def,$ss5_def);
GetOptions(
            'b|bed=s'          => \$bedFile,
            'f|fa=s'           => \$faFile,
	    'ss3'              => \$ss3_def,
	    'ss5'              => \$ss5_def,
            'h|help'           => sub{usage()}
	  ) || usage();
my $BED;
if(defined $bedFile){
    open $BED,"$bedFile" or die "Can't open file $bedFile:$!";
}else{
    $BED=\*STDIN;
}
my $db = Bio::DB::Fasta->new($faFile);
while(<$BED>){
    chomp;
    $tag+=1;
    my @split=split;
    my($chr,$start,$end,$strand)=@split[0..2,5];
    my($ss3_score,$ss5_score);
    if(defined $ss3_def){
	if($strand eq "+"){
	    my $ss3_seq=$db->seq($chr,$start-19,$start+3);#20 bases in intron,3 bases in exon
	    $ss3_score=`/mnt/share/liym/tools/MaxEnd/score3.pl -s $ss3_seq|cut -f2`;
	}else{
	    my $ss3_seq=$db->seq($chr,$end+20,$end-2,);
	    $ss3_score=`/mnt/share/liym/tools/MaxEnd/score3.pl -s $ss3_seq|cut -f2`;
	}
	chomp($ss3_score);
	say join "\t",@split,$ss3_score;
    }elsif(defined $ss5_def){
	if($strand eq "+"){
	    my $ss5_seq=$db->seq($chr,$end-2,$end+6); #3 bases in exon,6 bases in intron
	    $ss5_score=`/mnt/share/liym/tools/MaxEnd/score5.pl -s $ss5_seq|cut -f2`;
	}else{
	    my $ss5_seq=$db->seq($chr,$start+3,$start-5,);
	    $ss5_score=`/mnt/share/liym/tools/MaxEnd/score5.pl -s $ss5_seq|cut -f2`;
	}
	chomp($ss5_score);
	say join "\t",@split,$ss5_score;
    }else{
	if($strand eq "+"){
	    my $ss3_seq=$db->seq($chr,$start-19,$start+3);#20 bases in intron,3 bases in exon
	    $ss3_score=`/mnt/share/liym/tools/MaxEnd/score3.pl -s $ss3_seq|cut -f2`;
	    my $ss5_seq=$db->seq($chr,$end-2,$end+6); #3 bases in exon,6 bases in intron
	    $ss5_score=`/mnt/share/liym/tools/MaxEnd/score5.pl -s $ss5_seq|cut -f2`;
	}else{
	    my $ss5_seq=$db->seq($chr,$start+3,$start-5,);
	    $ss5_score=`/mnt/share/liym/tools/MaxEnd/score5.pl -s $ss5_seq|cut -f2`;
	    my $ss3_seq=$db->seq($chr,$end+20,$end-2,);
	    $ss3_score=`/mnt/share/liym/tools/MaxEnd/score3.pl -s $ss3_seq|cut -f2`;
	}
	chomp($ss3_score);
	chomp($ss5_score);
	say join "\t",@split,$ss5_score,$ss3_score;
    }
}
sub usage{
print <<HELP; 
Usage:	perl $0 -b file.bed -f hg19.fa >result.file 2>log
Author: Yumei Li,2015-8-17
Revision: Adding option of calculating only 3'ss or 5'ss score.
Output: input_fields 5'ss_score 3'ss_score
        Statistics the 5' and 3' splice site score for the given region in BED format file.(~/ToolKit/MaxEnd/score3.pl;~/ToolKit/MaxEnd/score5.pl)
        'b|bed'    FILE         The bed format file for the region.[can be bed6 or bed6+](If not given, it can be read from STDIN)
        'f|fa'     FILE         reference genome, fasta format
	'--3s'     LOGIC        Only calculate 3'ss score.
	'--5s'     LOGIC        Only calculate 5'ss score.
        'help|h'                print this help message    
HELP
    exit(-1);
}
