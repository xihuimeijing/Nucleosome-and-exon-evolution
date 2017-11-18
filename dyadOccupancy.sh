#!/bin/sh
usage(){
  cat << EOF
Description: Generate dayd occupancy form MNase-seq bam.
Dependencies: bamtools, samtools, bedgraph_to_wig.pl, tooRunner.sh wigToBigWig. Please ensure the above softwares are in your PATH.
Options:
    i   FILE    The orginal BAM format file.
    g   FILE    The tab delimited chromosome size file.
    d   LOGIC   Use DANPOS2 to call nucleosomes and calculate nucleosome occupancy[default:Use toolRunner.sh to calculate nucleosome dyad occupancy].
    h           Print this help information.
EOF
    exit 0
}
[ $1 ] || usage

while getopts "hi:g:d" OPTION
do
    case $OPTION in
        h) usage;;
        i) inputBam=$OPTARG;;
        g) chrSize=$OPTARG;;
        d) danpos=1;;
        ?) usage;;
    esac
done

samtools view $inputBam |grep -E "SA|XA" |cut -f1|sort|uniq >multAligName.txt
perl /mnt/share/liym/bin/bam_filterByReadname.pl -i $inputBam -r multAligName.txt |samtools view -@ 10 -bSu - |samtools sort -@ 10 - out.sorted.uniq
bamtools stats -insert -in out.sorted.uniq.bam >out.sorted.uniq.bamStats &
samtools rmdup out.sorted.uniq.bam out.sorted.uniq.rmdup.bam >samtoolsRmdup.log 2>samtoolsRmdup.err
samtools sort -@ 10 -n out.sorted.uniq.rmdup.bam out.uniq.rmdup.nameSorted
perl /mnt/share/liym/bin/bam_insert_filter.pl -i out.uniq.rmdup.nameSorted.bam -s 100,250 -q 30 -n 10 | samtools view -@ 10 -bSu - | samtools sort -@ 10 -n - out.uniq.rmdup.filter.sortedName
bamtools stats -insert -in out.uniq.rmdup.filter.sortedName.bam >out.uniq.rmdup.filter.sortedName.bamStats &
if [ $danpos ];then
    danpos.py dpos out.uniq.rmdup.filter.sortedName.bam -jd 147 -m 1 --extend 74 -o danpos_rst >danpos.log 2>danpos.err
    awk -v OFS="\t" '{if($1!="chr"){print $1,$4-1,$4}}' danpos_rst/pooled/*positions.xls >danpos_rst/pooled/nucleosome_dyad.bed3
    perl /mnt/share/liym/bin/wigNormalize.pl -i danpos_rst/pooled/*.wig -c $chrSize >danpos_rst/pooled/NC_norm.wig
    wigToBigWig -clip  danpos_rst/pooled/NC_norm.wig $chrSize danpos_rst/pooled/NC_norm.bw 2> danpos_rst/pooled/wigToBigwig.log
else
    bamToBed -bedpe -i out.uniq.rmdup.filter.sortedName.bam > fragment.bed 2> bamToBedpe.log
    cut -f 1,2,6,7 fragment.bed | awk -v OFS="\t" '{d=int(($3-$2)/2);print $1,$2+d,$2+d+1;}' | sort -k1,1 | bedtools genomecov -bg -i stdin -g $chrSize | sort -k1,1 -k2,2n > dyad.occupancy.bg
    bedgraph_to_wig.pl --bedgraph dyad.occupancy.bg --wig dyad.occupancy.wig --step 1
    toolRunner.sh wigmath.Scale -i dyad.occupancy.wig -f -p 10 -o dyad.occupancy.scale.wig >wigmath.Scale.log
    toolRunner.sh wigmath.GaussianSmooth -f -i dyad.occupancy.scale.wig -o dyad.occupancy.scale.smooth.wig -p 10 >wigmath.GaussianSmooth.log
    wigToBigWig dyad.occupancy.scale.smooth.wig $chrSize dyad.occupancy.scale.smooth.bw
    rm dyad.occupancy.bg dyad.occupancy.wig dyad.occupancy.scale.wig
fi