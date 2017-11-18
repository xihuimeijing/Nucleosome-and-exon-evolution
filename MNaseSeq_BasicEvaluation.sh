#!/bin/sh
usage(){
  cat << EOF
Description: Perform basic evaluation for MNase-seq bam and called nucleosomes by DANPOS2.
Dependencies: samtools, picard, RSeQC, bwtool and inhouse scripts peReadsFragmentMid.pl, region_dinucleotideFreq.pl, dinucleotide_plot.R, nucleosome_coverage.pl, TSS_nucleosomeCov.R.
Options:
    i   FILE    The orginal BAM format file infull path.
    b   FILE    The nucleosome occupancy in bw file. 
    g   FILE    The gene structure file in gpe format.
    d   DIR     The DANPOS2 result directory name.
    f   FILE    The genome file in fasta format.
    e   FILE    The internal exons in bed6/bed6+ format.
    l   INT     Read length of the input bam file.
    o   DIR     The output directory.[defalut:./]
    h           Print this help information.
EOF
    exit 0
}
[ $1 ] || usage

outDir=$PWD
perlScript=/mnt/share/liym/nucleosome/scripts
while getopts "hi:b:g:d:f:e:l:o:" OPTION
do
    case $OPTION in
        h) usage;;
        i) inputBam=$OPTARG;;
        b) occupancy=$OPTARG;;
        g) gpeStructure=$OPTARG;;
        d) dirDanpos=$OPTARG;;
        f) refFa=$OPTARG;;
        e) exonBed=$OPTARG;;
        l) readLen=$OPTARG;;
        o) outDir=$OPTARG;;
        ?) usage;;
    esac
done

#Get full path of input file and directory
inputBam=$(realpath $inputBam)
occupancy=$(realpath $occupancy)
gpeStructure=$(realpath $gpeStructure)
dirDanpos=$(realpath $dirDanpos)
refFa=$(realpath $refFa)
exonBed=$(realpath $exonBed)

cd $outDir
if [ -d evaluation ];then
    [ -d evaluation_old ] && rm -r evaluation_old
    mv evaluation evaluation_old
fi
mkdir evaluation && cd evaluation

if [ -d log ];then
    [ -d log_old ] && rm -r log_old
    mv log log_old
fi
mkdir log

echo "Fragment size distribution-----"
java -jar /mnt/share/share/tool/picard-tools/picard.jar CollectInsertSizeMetrics I=$inputBam O=insert_size_metrics.txt H=insert_size_histogram.pdf > log/collectInsertS.log 2> log/collectInsertS.err
echo "Calculate mismatch profile-----"
/mnt/share/liym/tools/RSeQC-2.6.3/install/home/share/local/bin/mismatch_profile.py -i $inputBam -l $readLen -o mismatchProfile 2>log/mismatchProfile.log &
echo "Calculate dinucleotide frequency across nucleosomes-----"
perl $perlScript/peReadsFragmentMid.pl -i $inputBam -f 147,147 -l $readLen >n147_mid.txt 
perl $perlScript/region_dinucleotideFreq.pl -b n147_mid.txt -n AA,AT,TA,TT -u -150 -d 150 -f $refFa -c >AT_freq.tsv 2>log/AT_freq.log &
perl $perlScript/region_dinucleotideFreq.pl -b n147_mid.txt -n GG,GC,CG,CC -u -150 -d 150 -f $refFa -c >GC_freq.tsv 2>log/GC_freq.log
if [ -s AT_freq.tsv && -s GC_freq.tsv ];then
    $perlScript/dinucleotide_plot.R AT_freq.tsv GC_freq.tsv dinucleotide.pdf 2>log/dinucleotide_R.log
fi
echo "Calculate nucleosome profile across TSS----"
awk '{if($3 == "+"){print $2"\t"$4-1500"\t"$4+1501"\t"$2"_"$4"\t0\t"$3;}else{print $2"\t"$5-1501"\t"$5+1500"\t"$2"_"$5-1"\t0\t"$3}}' $gpeStructure |awk -v OFS="\t" '{if($2<0){print $1,"0",$3,$4,$5,$6}}'|sort|uniq >TSS_updn1.5k.bed6
perl $perlScript/nucleosome_coverage.pl -b TSS_updn1.5k.bed6 -i 6 -n $dirDanpos/pooled/nucleosome_dyad.bed3 |sort -k1,1n >TSS_nucleosomeCov.tsv
Rscript $perlScript/TSS_nucleosomeCov.R TSS_nucleosomeCov.tsv
rm TSS_updn1.5k.bed6
$perlScript/TSSprofile_bwtool.sh -u 1500 -d 1500 -g $gpeStructure -b $occupancy -o TSS_profile.tsv
Rscript /mnt/share/liym/bin/lines.R -i=TSS_profile.tsv -x="TSS relative position" -y="Normalized nucleosome occupancy" -c="blue" -o=TSS_profile.pdf 
echo "Calculate nucleosome profile across exon/intron boundary----"
bwtool agg 500:200:500 $exonBed $occupancy exon_meta_profile.tsv
Rscript /mnt/share/liym/bin/lines.R -i=exon_meta_profile.tsv -x="Exon relative position" -y="Normalized nucleosome occupancy" -c="blue" -o=exon_meta_profile.tsv.pdf
bwtool agg 500:100 $exonBed $occupancy exon_up_profile.tsv -starts
Rscript /mnt/share/liym/bin/lines.R -i=exon_up_profile.tsv -x="Exon relative position" -y="Normalized nucleosome occupancy" -c="blue" -o=exon_up_profile.tsv.pdf
bwtool agg 100:100 $exonBed $occupancy exon_down_profile.tsv -starts
Rscript /mnt/share/liym/bin/lines.R -i=exon_down_profile.tsv -x="Exon relative position" -y="Normalized nucleosome occupancy" -c="blue" -o=exon_down_profile.tsv.pdf

