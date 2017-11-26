#!/bin/sh
#Author: Yumei Li, 2017/11/20
#E-mail: liyumei@pku.edu.cn
#This is the whole pipeline for nucleosome project processing.(Data are located at pluto and jupiter)
scriptPath=/mnt/share/liym/nucleosome/scripts
#1 Data processing and evaluation
##1.1 MNase-seq data
cd /rd1/user/liym/nucleosome # @pluto
bwaDB=/data/bwa/hg19/bwtsw #For human
bwaDB=/data/bwa/rheMac2/bwtsw #For monkey
bwaDB=/mnt/share/liym/data/bwa/tupBel1/tupBel1 #Fro tree shrew
bwaDB=/data/bwa/mm9/bwtsw # For mouse
bwaDB=/mnt/share/liym/data/genome/susScr3/bwa/susScr3.fa # For pig
####Human as an example
ls human_tissue*/*/*_1.fq | while read file;do
    dir=$(dirname $file);
    cd $dir;
    prefix=$(basename $file|sed 's/_1.clean.fq//');
    fqPeFilter.pl -f 0.05 -b N -l 100 -q 20 -a 10 -t 5,5 --o1 ${prefix}_1.filter.fq --o2 ${prefix}_2.filter.fq ${prefix}_1.clean.fq ${prefix}_2.clean.fq 2>fqFilter.log;
    bwa mem -t 10 $bwaDB ${prefix}_1.filter.fq ${prefix}_2.filter.fq 2>bwa.err| samtools view -bSu - | samtools sort -@ 8 - out.sorted 2> samtoolsSort.err
    $scriptPath/dyadOccupancy.sh -d -i out.sorted.bam -g /data/chr.size/hg19.size >dyadOccupancy.log 2> dyadOccupancy.err
    $scriptPath/MNaseSeq_BasicEvaluation.sh -i out.uniq.rmdup.filter.sortedName.bam -b danpos_rst/pooled/NC_norm.bw -g /mnt/share/liym/data/structure/hg19/2016-09-24.refGene.noBin.gpe -d danpos_rst/ -f /mnt/share/liym/data/genome/hg19/fa/hg19.fa -e /mnt/share/liym/data/structure/hg19/hg19_refGene_internal_uniqExon.bed6+ -l 150 > MNaseSeq_BasicEvaluation.log 2> MNaseSeq_BasicEvaluation.err
    cd ../../;
done;

##1.2 RNA-seq data
cd /rd1/user/liym/nucleosome/RNA-seq # @pluto
ls hg19/*/*_1.fq |while read file;do
    dir=$(dirname $file);
    prefix=$(basename $file|sed 's/_1.fq//');
    cd $dir;
    fqTrimer.pl -l [ATGCN]{10} --i2 ${prefix}_2.fq --o2 readTrimer_2.fq ${prefix}_1.fq >readTrimer_1.fq 2> fqTrimer.log;
    fqPeFilter.pl -f 0.05 -b N -l 100 -q 0.5 -c 20 -a 20 --o1 readFilter_1.fq --o2 readFilter_2.fq readTrimer_1.fq readTrimer_2.fq 2>fqFilter.log;
    tophat2 -p 10 -N 8 --read-edit-dist 8 -a 8 --library-type fr-firststrand --no-coverage-search --segment-length 28 --segment-mismatches 2 /data/bowtie2/hg19/hg19 readFilter_1.fq readFilter_2.fq >tophat2.log 2>tophat2.err;
    samtools view -bu -q 50 -@ 5 tophat_out/accepted_hits.bam | samtools sort -@ 5 -m 10G - uniq.sorted
    samtools index uniq.sorted.bam
    ln -s *_1_fastqc final1_fastqc
    ln -s *_2_fastqc final2_fastqc
    $scriptPath/rnaSeqPe.sh -t 2 -e -f /data/fna/hg19/all.fa -g /mnt/share/liym/data/structure/hg19/2016-09-24.refGene.noBin.gpe -s 4 -l fr-firststrand -p /home/share/tool/picard-tools/ >runEva.log 2>runEva.err
    cd ../../;
done;
ls rheMac2/*/*_1.fq |while read file;do
    dir=$(dirname $file);
    prefix=$(basename $file|sed 's/_1.fq//');
    cd $dir;
    fqTrimer.pl -l [ATGCN]{10} --i2 ${prefix}_2.fq --o2 readTrimer_2.fq ${prefix}_1.fq >readTrimer_1.fq 2> fqTrimer.log;
    fqPeFilter.pl -f 0.05 -b N -l 100 -q 0.5 -c 20 -a 20 --o1 readFilter_1.fq --o2 readFilter_2.fq readTrimer_1.fq readTrimer_2.fq 2>fqFilter.log;
    tophat2 -p 15 -N 8 --read-edit-dist 8 -a 8 --library-type fr-firststrand --no-coverage-search --segment-length 28 --segment-mismatches 2 /mnt/share/liym/data/bowtie2/rheMac2/rheMac2 readFilter_1.fq readFilter_2.fq >tophat2.log 2>tophat2.err
    samtools view -bu -q 50 -@ 5 tophat_out/accepted_hits.bam | samtools sort -@ 5 -o uniq.sorted.bam -
    samtools index uniq.sorted.bam
    ln -s *_1_fastqc final1_fastqc
    ln -s *_2_fastqc final2_fastqc
    $scriptPath/rnaSeqPe.sh -t 2 -e -f /mnt/share/liym/data/genome/rheMac2/rheMac2.fa -g /mnt/share/share/data/structure/gpe/rheMac2.rb2Gene.gpe -s 4 -l fr-firststrand -p /home/share/tool/picard-tools/ >runEva.log 2>runEva.err
    cd ../../;
done;
ls tupBel1/*/*_1.fq|while read file;do
    dir=$(dirname $file);
    prefix=$(basename $file|sed 's/_1.fq//');
    cd $dir;
    fqTrimer.pl -l [ATGCN]{10} --i2 ${prefix}_2.fq --o2 readTrimer_2.fq ${prefix}_1.fq >readTrimer_1.fq 2> fqTrimer.log;
    fqPeFilter.pl -f 0.05 -b N -l 100 -q 0.5 -c 20 -a 20 --o1 readFilter_1.fq --o2 readFilter_2.fq readTrimer_1.fq readTrimer_2.fq 2>fqFilter.log;
    tophat2 -N 10 --read-edit-dist 10 --segment-mismatches 3 --segment-length 28 -a 8 --library-type fr-firststrand -p 10 /data/bowtie2/tupBel1/tupBel1 readFilter_1.fq readFilter_2.fq >tophat2.log 2>tophat2.err;
    samtools view -bu -q 50 -@ 5 tophat_out/accepted_hits.bam | samtools sort -@ 5 -o uniq.sorted.bam -
    samtools index uniq.sorted.bam
    ln -s *_1_fastqc final1_fastqc
    ln -s *_2_fastqc final2_fastqc
    /mnt/share/liym/nucleosome/scripts/rnaSeqPe.sh -t 2 -e -p /home/share/tool/picard-tools/ -f /data/fna/tupBel1/tupBel1.fa -g /mnt/share/liym/data/structure/tupBel1/tupBel1.ensGene.noBin.gpe -l fr-firststrand >runEva.log 2>runEva.err
    cd ../../;
done;
ls mm9/*/*_1.fq|while read file;do
    dir=$(dirname $file);
    prefix=$(basename $file|sed 's/_1.fq//');
    cd $dir;
    fqTrimer.pl -l [ATGCN]{10} --i2 ${prefix}_2.fq --o2 readTrimer_2.fq ${prefix}_1.fq >readTrimer_1.fq 2> fqTrimer.log;
    fqPeFilter.pl -f 0.05 -b N -l 100 -q 0.5 -c 20 -a 20 --o1 readFilter_1.fq --o2 readFilter_2.fq readTrimer_1.fq readTrimer_2.fq 2>fqFilter.log;
    tophat2 -p 16 -N 8 --read-edit-dist 8 -a 8 --library-type fr-firststrand --no-coverage-search --segment-length 28 --segment-mismatches 2 /mnt/share/share/data/bowtie2/mm9/mm9 readFilter_1.fq readFilter_2.fq>tophat2.log 2>tophat2.err
    samtools view -bu -q 50 -@ 5 tophat_out/accepted_hits.bam | samtools sort -@ 5 -o uniq.sorted.bam -
    samtools index uniq.sorted.bam
    ln -s *_1_fastqc final1_fastqc
    ln -s *_2_fastqc final2_fastqc
    /mnt/share/liym/nucleosome/scripts/rnaSeqPe.sh -t 2 -e -f /mnt/share/liym/data/genome/mm9/mm9.fa -g /mnt/share/share/data/structure/gpe/mm9.refGene.gpe -s 4 -l fr-firststrand -p /home/share/tool/picard-tools/ >runEva.log 2>runEva.err
    cd ../../;
done;
ls susScr3/*/*_1.fq|while read file;do
   dir=$(dirname $file);
    prefix=$(basename $file|sed 's/_1.fq//');
    cd $dir;
    fqTrimer.pl -l [ATGCN]{10} --i2 ${prefix}_2.fq --o2 readTrimer_2.fq ${prefix}_1.fq >readTrimer_1.fq 2> fqTrimer.log;
    fqPeFilter.pl -f 0.05 -b N -l 100 -q 0.5 -c 20 -a 20 --o1 readFilter_1.fq --o2 readFilter_2.fq readTrimer_1.fq readTrimer_2.fq 2>fqFilter.log;
    tophat2 --read-mismatches 10 --read-gap-length 3 --read-edit-dist 10 --min-anchor 8 --splice-mismatches 0 --num-threads 15 --mate-inner-dist 0 --no-coverage-search --segment-length 30 --segment-mismatches 2 --library-type fr-firststrand -o tophat_out /mnt/share/liym/data/bowtie2/susScr3/susScr3 readFilter_1.fq readFilter_2.fq >tophat2.log 2>tophat2.err
    samtools view -bu -q 50 -@ 5 tophat_out/accepted_hits.bam | samtools sort -@ 5 -o uniq.sorted.bam -
    samtools index uniq.sorted.bam
    ln -s *_1_fastqc final1_fastqc
    ln -s *_2_fastqc final2_fastqc
    /mnt/share/liym/nucleosome/scripts/rnaSeqPe.sh -t 2 -e -p /home/share/tool/picard-tools/ -f /mnt/share/liym/data/genome/susScr3/susScr3.fa -g /mnt/share/liym/data/structure/susScr3/susScr3.ensGene.noBin.gpe -l fr-firststrand >runEva.log 2>runEva.err
    cd ../../;
done;
###Generate coverage and junction files
ls -d */*/|while read dir;do
    species=$(echo $dir|cut -f1 -d '/');
    cd $dir;
    mkdir junc
    sam2junc.pl -l fr-firststrand -p uniq.sorted.bam >junc/junc.bed12 2>junc/junc.sam
    mkdir cov
    samtools view uniq.sorted.bam | wiggles /dev/stdin cov/uniq.bg;
    bedGraphToBigWig cov/uniq.bg /data/chr.size/${speceies}.size cov/coverage.bw;
    cd ../../;
done;

#2 Whole genome comparison
cd /rd1/user/liym/nucleosome/mutiSpecies/
##2.1 Prepare data
mkdir data && cd data
mkdir hg19 rheMac2 tupBel1 mm9 susScr3
ln -s /rd1/user/liym/nucleosome/human_tissues/brain/danpos_rst/pooled/brain_final.sorted.smooth.norm.bw hg19/brain.norm.bw
ln -s /rd1/user/liym/nucleosome/human_tissues/muscle/danpos_rst/pooled/muscle_final.sorted.smooth.norm.bw hg19/muscle.norm.bw

ln -s /rd1/user/liym/nucleosome/monkey_tissues/brain/danpos_rst/pooled/monkey_tissues_brain_final.sorted.smooth.norm.bw rheMac2/brain.norm.bw
ln -s /rd1/user/liym/nucleosome/monkey_tissues/heart/danpos_rst/pooled/monkey_tissues_heart_final.sorted.smooth.norm.bw rheMac2/heart.norm.bw
ln -s /rd1/user/liym/nucleosome/monkey_tissues/kidney/danpos_rst/pooled/monkey_tissues_kidney_final.sorted.smooth.norm.bw rheMac2/kidney.norm.bw
ln -s /rd1/user/liym/nucleosome/monkey_tissues/liver/danpos_rst/pooled/monkey_tissues_liver_final.sorted.smooth.norm.bw rheMac2/liver.norm.bw
ln -s /rd1/user/liym/nucleosome/monkey_tissues/muscle/danpos_rst/pooled/monkey_tissues_muscle_final.sorted.smooth.norm.bw rheMac2/muscle.norm.bw

ln -s /rd1/user/liym/nucleosome/treeShrew_tissues/brain/danpos_rst/pooled/treeShrew_tissues_brain_final.sorted.smooth.norm.bw tupBel1/brain.norm.bw
ln -s /rd1/user/liym/nucleosome/treeShrew_tissues/heart/danpos_rst/pooled/treeShrew_tissues_heart_final.sorted.smooth.norm.bw tupBel1/heart.norm.bw
ln -s /rd1/user/liym/nucleosome/treeShrew_tissues/kidney/danpos_rst/pooled/treeShrew_tissues_kidney_final.sorted.smooth.norm.bw tupBel1/kidney.norm.bw
ln -s /rd1/user/liym/nucleosome/treeShrew_tissues/muscle/danpos_rst/pooled/treeShrew_tissues_muscle_final.sorted.smooth.norm.bw tupBel1/muscle.norm.bw
##Mouse data moved from Jupiter
ln -s mm9/brain_final.sorted.smooth.norm.bw mm9/brain.norm.bw 
ln -s mm9/brain_final.sorted.smooth.norm.bw mm9/heart.norm.bw
ln -s mm9/brain_final.sorted.smooth.norm.bw mm9/liver.norm.bw
ln -s mm9/brain_final.sorted.smooth.norm.bw mm9/kidney.norm.bw
ln -s mm9/brain_final.sorted.smooth.norm.bw mm9/muscle.norm.bw

ln -s /rd1/user/liym/nucleosome/pig_tissue_data/brain/danpos_rst/pooled/pig_tissue_data_brain_final.sorted.smooth.norm.bw susScr3/brain.norm.bw
ln -s /rd1/user/liym/nucleosome/pig_tissue_data/heart/danpos_rst/pooled/pig_tissue_data_heart_final.sorted.smooth.norm.bw susScr3/heart.norm.bw
ln -s /rd1/user/liym/nucleosome/pig_tissue_data/kidney/danpos_rst/pooled/pig_tissue_data_kidney_final.sorted.smooth.norm.bw susScr3/kidney.norm.bw
ln -s /rd1/user/liym/nucleosome/pig_tissue_data/liver/danpos_rst/pooled/pig_tissue_data_liver_final.sorted.smooth.norm.bw susScr3/liver.norm.bw
ln -s /rd1/user/liym/nucleosome/pig_tissue_data/muscle/danpos_rst/pooled/pig_tissue_data_muscle_final.sorted.smooth.norm.bw susScr3/muscle.norm.bw
##2.2 Nucleosome Position
cd /rd1/user/liym/nucleosome/mutiSpecies/position
###Nucleosome dyad position from other species to human
ls /rd1/user/liym/nucleosome/monkey_tissues/*/danpos_rst/pooled/*final*xls|while read file;do tissue=$(echo $file|cut -f7 -d '/');grep -vE "start|M|random|Un|gl" $file|awk -v OFS="\t" '{print $1,$4-74,$4+73}' >monkey.tmp;liftOver monkey.tmp /data/liftover/rheMac2/rheMac2ToHg19.over.chain monkey_${tissue}.hg19.bed3 monkey_${tissue}.unmapped;done;
ls /rd1/user/liym/nucleosome/mouse_tissues/*/*final*xls|while read file;do tissue=$(echo $file|cut -f7 -d '/');grep -vE "start|M|random|Un|gl" $file|awk -v OFS="\t" '{print $1,$4-74,$4+73}' >mouse.tmp;liftOver mouse.tmp /data/liftover/mm9/mm9ToHg19.over.chain mouse_${tissue}.hg19.bed3 mouse_${tissue}.unmapped;done;
ls /rd1/user/liym/nucleosome/pig_tissue_data/*/danpos_rst/pooled/*final*xlswhile read file;do tissue=$(echo $file|cut -f7 -d '/');grep -vE "start|M|random|Un|gl" $file|awk -v OFS="\t" '{print $1,$4-74,$4+73}' >pig.tmp;liftOver pig.tmp /mnt/share/liym/data/liftOver/susScr3ToHg19.over.chain pig_${tissue}.hg19.bed3 pig_${tissue}.unmapped;done;
ls /rd1/user/liym/nucleosome/treeShrew_tissues/*/danpos_rst/pooled/*final*xls|while read file;do tissue=$(echo $file|cut -f7 -d '/');grep -vE "start|M|random|Un|gl" $file|awk -v OFS="\t" '{print $1,$4-74,$4+73}' >tree.tmp;liftOver tree.tmp /mnt/share/liym/data/liftOver/tupBel1ToHg19.over.chain tree_${tissue}.hg19.bed3 tree_${tissue}.unmapped;done;
####n147 homolog regions
cd /mnt/share/liym/data/conservation/5way/
cat /data/conservation/multiz/hg19/chr*.maf |mafSpeciesSubset stdin species.list 4species.maf 2>mafSpeciesSubset.err 
perl /mnt/share/liym/bin/mafTobed.pl -m 4species.maf >4species.bed6+ 
awk 'NF==20' 4species.bed6+ >4species.filtered.bed6+
cat ../hg19_susScr3/*.maf |perl /mnt/share/liym/bin/mafTobed.pl >hg19.susScr3.bed6+
cd /rd1/user/liym/nucleosome/mutiSpecies/occupancy
ln -s /mnt/share/liym/data/conservation/5way/hg19.susScr3.bed6+ hg19.susScr3.bed
ln -s /mnt/share/liym/data/conservation/5way/4species.filtered.bed6+ 4species.coor.bed
cut -f1-3 hg19.susScr3.bed |sort|uniq >hg19.susScr3.uniq.bed3
cut -f1-3 4species.coor.bed |bedtools intersect -wo -a stdin -b hg19.susScr3.uniq.bed3 >5species.intersect.bed
awk -v ORS="" '{print $1"\t";if($2>$5){print $2"\t"}else{print $5"\t"};if($3>$6){print $6"\t"}else{print $3}print "\n"}' 5species.intersect.bed |sort|uniq >5species.homo.regions.bed3 
sort -k1,1 -k2,2n 5species.homo.regions.bed3 |bedtools merge -i stdin >5species.homo.merged.bed
cd ../position
perl /mnt/share/liym/bin/bedToIntervals.pl -b ../occupancy/5species.homo.merged.bed -l 147 -d 147 >5species.n147.bed3
#Intersect
mkdir intersect
ls *.hg19.bed3 |while read file;do 
        prefix=$(echo $file|sed 's/.hg19.bed3//');
        awk -v OFS="\t" '{if($3-$2>=100 && $3-$2<=200){print $1,$2+int(($3-$2)/2),$2+int(($3-$2)/2)+1}}' $file|intersectBed -loj -a 5species.n147.bed3 -b stdin >intersect/${prefix}.bed 
done;
grep -vE "start|M|random|Un|gl" /rd1/user/liym/nucleosome/human_tissues/brain/danpos_rst/pooled/brain_final.sorted.smooth.positions.xls|awk -v OFS="\t" '{print $1,$4-1,$4}' |intersectBed -loj -a 5species.n147.bed3 -b stdin >intersect/human_brain.bed 
grep -vE "start|M|random|Un|gl" /rd1/user/liym/nucleosome/human_tissues/muscle/danpos_rst/pooled/muscle_final.sorted.smooth.positions.xls|awk -v OFS="\t" '{print $1,$4-1,$4}' |intersectBed -loj -a 5species.n147.bed3 -b stdin >intersect/human_muscle.bed 
cd intersect
ls *.bed|sed 's/.bed//'|while read file;do cut -f1-3 ${file}.bed |sort|uniq -c|awk -v OFS="\t" '{if($1==1){print $2,$3,$4}}'|bedtools intersect -wa -f 1.0 -r -a ${file}.bed -b stdin >${file}_filter.bed3+;done;
bedtools intersect -wo -f 1.0 -r -a human_brain_filter.bed3+ -b human_muscle_filter.bed3+|cut -f1-6,10-12 >summary.tsv
tag=9
ls *.bed3+|grep -vE "human_muscle|human_brain"|while read file;do
        start=$(($tag+4))
        end=$(($tag+6))
        bedtools intersect -wo -f 1.0 -r -a summary.tsv -b $file |cut -f1-${tag},${start}-${end} >tmp;
        mv tmp summary.tsv
        tag=$(($tag+3))
done;
##assign distance
perl -e '
        open IN,"$ARGV[0]";
        while(<IN>){
                chomp;
                my @split=split /\t/;
                print join "\t",@split[0..2];
                for(my $i=4;$i<=54;$i+=3){
                        if($split[$i]<0){
                                print "\tNA";
                        }else{
                                my $dis=$split[$i]-$split[1];
                                print "\t$dis";
                        }
                }
                print "\n";
        }
' summary.tsv >summary.distance.tsv 
awk '{na=0;for(i=4;i<=20;i++){if($i=="NA"){na+=1;}}if(na<=15){print $0}}' summary.distance.tsv >summary.distance.filter.tsv
#####Start R#####
data<-read.delim(file="summary.distance.filter.tsv",header=F)
data<-data[,-c(1,2,3)]
colnames(data)<-c("human_brain","human_muscle","monkey_brain","monkey_heart","monkey_kidney","monkey_liver","monkey_muscle","mouse_GSM1399400","mouse_GSM717558","mouse_brain","mouse_heart","mouse_kidney","mouse_liver","mouse_muscle","mouse_sperm","pig_brain","pig_heart","pig_kidney","pig_liver","pig_muscle","tree_brain","tree_heart","tree_kidney","tree_muscle")
cor<-cor(data,use="pairwise.complete.obs")
write.table(cor,file="homo.n147.pearsonCor.tsv",sep="\t",append=F,quote=F)
####End R#####
##2.3 Heat map for human and monkey
cd /rd1/user/liym/nucleosome/mutiSpecies/position/heatmap
comm -12 <(grep -vE "random|Un|gl|hap|M|summit" ../../../human_tissues/brain/danpos_rst/pooled/brain_final.sorted.smooth.positions.xls|cut -f1,4|sort) <(grep -vE "random|Un|gl|hap|M|summit" ../../../human_tissues/muscle/danpos_rst/pooled/muscle_final.sorted.smooth.positions.xls|cut -f1,4|sort) |awk -v OFS="\t" '{print $1,$2-1,$2,$1"_"$2}' >human.consensus.nucleosome.bed3 
liftOver human.consensus.nucleosome.bed3 /data/liftover/hg19/hg19ToRheMac2.over.chain hg19TorheMac2.consensus.nucleosome.bed3 hg19TorheMac2.unmapped
CrossMap.py bigwig /data/liftover/hg19/hg19ToRheMac2.over.chain /rd1/user/liym/nucleosome/human_tissues/brain/danpos_rst/pooled/brain_final.sorted.smooth.norm.bw human.brainTorheMac2
CrossMap.py bigwig /data/liftover/hg19/hg19ToRheMac2.over.chain /rd1/user/liym/nucleosome/human_tissues/muscle/danpos_rst/pooled/muscle_final.sorted.smooth.norm.bw human.muscleTorheMac2
awk -v OFS="\t" '{if($1!="chr"){print $1"_"$4,$6}}' ../../../human_tissues/brain/danpos_rst/pooled/brain_final.sorted.smooth.positions.xls |join.pl -i1 human.consensus.nucleosome.bed3 -f1 4 |cut -f1-4,6 >tmp.bed 
awk -v OFS="\t" '{if($1!="chr"){print $1"_"$4,$6}}' ../../../human_tissues/muscle/danpos_rst/pooled/muscle_final.sorted.smooth.positions.xls|join.pl -i1 tmp.bed -f1 4|awk -v OFS="\t" '{print $1,$2,$3,$4,($5+$7)/2}'|sort -k5,5n >human.consensus.nucleosome.bed4+
join.pl -i1 hg19TorheMac2.consensus.nucleosome.bed3 -f1 4 -i2 human.consensus.nucleosome.bed4+ -f2 4 -o1 >hg19TorheMac2.consensus.nucleosome.fuzzinessSort.bed4 
computeMatrix reference-point -R hg19TorheMac2.consensus.nucleosome.fuzzinessSort.bed4 --missingDataAsZero -S human.brainTorheMac2.bw human.muscleTorheMac2.bw /rd1/user/liym/nucleosome/monkey_tissues/brain/danpos_rst/pooled/monkey_tissues_brain_final.sorted.smooth.norm.bw /rd1/user/liym/nucleosome/monkey_tissues/heart/danpos_rst/pooled/monkey_tissues_heart_final.sorted.smooth.norm.bw /rd1/user/liym/nucleosome/monkey_tissues/kidney/danpos_rst/pooled/monkey_tissues_kidney_final.sorted.smooth.norm.bw /rd1/user/liym/nucleosome/monkey_tissues/liver/danpos_rst/pooled/monkey_tissues_liver_final.sorted.smooth.norm.bw /rd1/user/liym/nucleosome/monkey_tissues/muscle/danpos_rst/pooled/monkey_tissues_muscle_final.sorted.smooth.norm.bw --referencePoint center -b 200 -a 200 -p 10 -out nucleosome.dyad.fl200.NC.matrix.gz
plotHeatmap -m nucleosome.dyad.fl200.NC.matrix.gz -out human.monkey.nucleosome.dyad.fl200.NC.heatmap.v5.pdf --refPointLabel "0" --colorMap Blues --sortRegions no --whatToShow "heatmap and colorbar" --samplesLabel H_brain H_muscle M_brain M_heart M_kidney M_liver M_muscle --missingDataColor blue 
#3 Calculate exon-intron NC log ratio and GC log ratio
## 3.1 GC log ratio(GC difference) for 5 species
cd /mnt/share/liym/data/structure
perl $scriptPath/NC_difference.pl -b hg19/hg19_refGene_internal_uniqExon.bed6+ -w /mnt/share/liym/data/gcContent/hg19.gc5Base.bw -s 6 -u 150 -d 150 >hg19/hg19_refGene_internal_uniqExon_GCdiff.bed6+
perl $scriptPath/NC_difference.pl -b rheMac2/rheMac2.ensGene.internal.uniqExon.bed6+ -w /mnt/share/liym/data/gcContent/rheMac2.gc5Base.bw -s 6 -u 150 -d 150 >rheMac2/rheMac2.ensGene.internal.uniqExon.GC.diff.bed6+
perl $scriptPath/NC_difference.pl -b tupBel1/tupBel1.ensGene.internal.uniqExon.bed6+ -w /mnt/share/liym/data/gcContent/tupBel1.gc5Base.bw -s 6 -u 150 -d 150 >tupBel1/tupBel1.ensGene.internal.uniqExon.GC.diff.bed6+
perl $scriptPath/NC_difference.pl -b mm9/mm9_refGene_internal_uniqExon.bed6 -w /mnt/share/liym/data/gcContent/mm9.gc5Base.bw -s 6 -u 150 -d 150 > mm9/mm9_refGene_internal_uniqExon.GCdiff.bed6+
perl $scriptPath/NC_difference.pl -b susScr3/susScr3.ensGene.internal.uniqExon.bed6+ -w /mnt/share/liym/data/gcContent/susScr3.gc5Base.bw -s 6 -u 150 -d 150 >susScr3/susScr3.ensGene.internal.exon.GC.Diff.tsv
## 3.2 NC log ratio
cd /rd1/user/liym/nucleosome/mutiSpecies/exon_intron
mkdir hg19 rheMac2 tupBel1 mm9 susScr3
cd hg19
perl /mnt/share/liym/nucleosome/scripts/NC_log_ratio.pl -b /mnt/share/liym/data/structure/hg19/hg19_refGene_internal_uniqExon.bed6+ -w /rd1/user/liym/nucleosome/human_tissues/muscle/danpos_rst/pooled/muscle_final.sorted.smooth.norm.bw -p 0.00001 -s 6 -u 150 -d 150 > muscle_logNC.tsv 2> logNC.log &
perl /mnt/share/liym/nucleosome/scripts/NC_log_ratio.pl -b /mnt/share/liym/data/structure/hg19/hg19_refGene_internal_uniqExon.bed6+ -w /rd1/user/liym/nucleosome/human_tissues/brain/danpos_rst/pooled/brain_final.sorted.smooth.bw -p 0.00001 -s 6 -u 150 -d 150 > brain_logNC.tsv 2> logNC.log
paste <(sort -k4,4 brain_logNC.tsv) <(sort -k4,4 muscle_logNC.tsv) |cut -f1-6,11,23 |sort|uniq >2tissues.uplogNC.tsv
paste <(sort -k4,4 brain_logNC.tsv) <(sort -k4,4 muscle_logNC.tsv) |cut -f1-6,12,24 |sort|uniq >2tissues.downlogNC.tsv
cd ../rheMac2
ls /rd1/user/liym/nucleosome/monkey_tissues/*/danpos_rst/pooled/*final*norm.bw|while read file;do
        tissue=$(echo $file|cut -f7 -d '/')
perl /mnt/share/liym/nucleosome/scripts/NC_log_ratio.pl -b /mnt/share/liym/data/structure/rheMac2/rheMac2.ensGene.internal.uniqExon.bed6+ -w $file -p 0.00001 -s 6 -u 150 -d 150 >${tissue}_logNC.tsv 2>logNC.log
done;
paste *_logNC.tsv |cut -f1-6,11,23,35,47,59 >5tissue.uplogNC.tsv
paste *_logNC.tsv |cut -f1-6,12,24,36,48,60 >5tissue.downlogNC.tsv
cd ../tupBel1
ls /rd1/user/liym/nucleosome/treeShrew_tissues/*/danpos_rst/pooled/*final*norm.bw |while read file;do tissue=$(echo $file|cut -f7 -d '/');perl /mnt/share/liym/nucleosome/scripts/NC_log_ratio.pl -b /mnt/share/liym/data/structure/tupBel1/tupBel1.ensGene.internal.uniqExon.bed6+ -w $file -p 0.00001 -s 6 -u 150 -d 150 >${tissue}_logNC.tsv 2> logNC.log;done;
paste *_logNC.tsv |cut -f1-6,11,23,35,47 >4tissues.uplogNC.tsv
paste *_logNC.tsv |cut -f1-6,12,24,36,48 >4tissues.downlogNC.tsv
cd ../mm9
ls /rd1/user/liym/nucleosome/mouse_tissues/*/*bw |while read file;do tissue=$(echo $file|cut -f7 -d '/');perl /mnt/share/liym/nucleosome/scripts/NC_log_ratio.pl -b /mnt/share/liym/data/structure/mm9/mm9_refGene_internal_uniqExon.bed6 -w $file -p 0.00001 -s 6 -u 150 -d 150 > ${tissue}_logNC.tsv 2> logNC.log;done;
ls *.tsv|while read file;do awk -v OFS="\t" '{print $1,$2,$3,$1":"$4":"$2"-"$3,"0",$4,$6,$7}' $file >tmp;mv tmp $file;done;
paste *_logNC.tsv|cut -f1-6,7,15,23,31,39 >5tissue.uplogNC.tsv
paste *_logNC.tsv|cut -f1-6,8,16,24,32,40 >5tissue.downlogNC.tsv
cd ../susScr3
ls /rd1/user/liym/nucleosome/pig_tissue_data/*/danpos_rst/pooled/*final*norm.bw|while read file;do tissue=$(echo $file|cut -f7 -d '/'); perl /mnt/share/liym/nucleosome/scripts/NC_log_ratio.pl -b /mnt/share/liym/data/structure/susScr3/susScr3.ensGene.internal.uniqExon.bed6+ -w $file -p 0.00001 -s 6 -u 150 -d 150 > ${tissue}_logNC.tsv 2> logNC.log;done;
paste *_logNC.tsv|cut -f1-6,7,15,23,31,39 >5tissues.uplogNC.tsv 
paste *_logNC.tsv|cut -f1-6,8,16,24,32,40 >5tissues.downlogNC.tsv 
## 3.3 logNC & logGC
mkdir logNC_logGC & cd logNC_logGC
ls ../*/*_logNC.tsv |grep -v "susScr3"|while read file;do 
        species=$(echo $file|cut -f2 -d '/');
        tissue=$(echo $file|cut -f3 -d '/'|sed 's/_logNC.tsv//');
        join.pl -i1 /rd1/brick/liym/data/structure/${species}/*GC*diff.bed6+ -f1 4 -i2 $file -f2 4 |cut -f4,14,26 >${species}-${tissue}.GCdiff.logNC.tsv;
done;
ls ../susScr3/*_logNC.tsv|while read file;do
            tissue=$(echo $file|cut -f3 -d '/'|sed 's/_ss_logNC.tsv//');
                join.pl -i1 /rd1/brick/liym/data/structure/susScr3/susScr3.ensGene.internal.exon.GC.Diff.SS.tsv -f1 4 -i2 $file -f2 4 |cut -f4,14,26 >susScr3-${tissue}.GCdiff.logNC.tsv;
done;
ls *.tsv|while read file;do sort $file|uniq|awk -v OFS="\t" '{if($2<=0 && $2>(-20)){print $1,$2,$3,(int($2/10)-1)*10"-"int($2/10)*10}else if($2<=(-20)){print $1,$2,$3,"<=-20"}else if($2>0 && $2<30){print $1,$2,$3,int($2/10)*10"-"(int($2/10)+1)*10}else{print $1,$2,$3,">=30"}}' >tmp;mv tmp $file;done;

#4 Identification of human new exons @Jupiter
##4.1 Data preparation
cd ~/nucleosome/assignExonAge
###4.1.1 Gene annotation
mkdir species_exon && cd species_exon
mkdir hg19 rheMac2 tupBel1 mm9 susScr3
ls -d */|sed 's;/;;'|while read dir;do 
        cd $dir;
        ls /mnt/share/liym/data/structure/${dir}/*.gpe|grep -v "Filtered"|while read file;do
                prefix=$(basename $file|cut -f2 -d '.');
                perl /mnt/share/liym/bin/gpeGetExons.pl -i $file >${prefix}.exon.bed6
        done;
        cd ../;
done;
ls -d */|while read dir;do
        cd $dir;
        cat *.exon.bed6|awk -v OFS="\t" '{print $1,$2,$3,$1":"$6":"$2"-"$3,$5,$6}' |sort|uniq >merged.uniq.exon.bed 
        cd ../
        done;
##### hg19 internal exons
perl ~/bin/gpeGetInternalExons.pl -i /mnt/share/liym/data/structure/hg19/gff/all.unique.gpe |awk -v OFS="\t" '{print $1,$2,$3,$1":"$6":"$2"-"$3,$5,$6}' |sort|uniq >hg19.all.internal.exon.bed6
###4.1.2 RNA-seq data
mkdir GTEx_hg19 rheMac2 tupBel1 mm9 susScr3
####Generate junction and coverage file & copy from pluto 
##4.2 Get human initiatial exons list
    #Get internal unique exon list
cd ~/nucleosome/assignExonAge
mkdir list_v3 && cd list_v3
ln -s ../species_exon/hg19.all.internal.exon.bed6 hg19.initial.exon.bed6 
##4.3 Lift human exons to other four species and assign exon status
mkdir pairwiseAln && cd pairwiseAln
###4.3.1 lift hg19 to rheMac2, tupBel1, mm9, susScr3
mafLiftOver2.pl -r hg19.initial.exon.bed6 -q rheMac2 -f bed6 -o hg19_rheMac2.bed6+ -q rheMac2 -f bed6 -o hg19_rheMac2.bed6+ /data/conservation/pairwise_alignment/hg19/rheMac2/net.maf/all.maf > hg19_rheMac2.tsv
cat /data/conservation/pairwise_alignment/hg19/TupBel1/*.maf | mafLiftOver2.pl -r hg19.initial.exon.bed6 -q tupBel1 -f bed6 -o hg19_tupBel1.bed6+ > hg19_tupBel1.tsv
cat /data/conservation/pairwise_alignment/hg19/mm9/*.maf | mafLiftOver2.pl -r hg19.initial.exon.bed6 -q mm9 -f bed6 -o hg19_mm9.bed6+ > hg19_mm9.tsv
cat ~/data/conservation/pairwise_alignment/hg19/susScr3/*.maf | mafLiftOver2.pl -r hg19.initial.exon.bed6 -q susScr3 -f bed6 -o hg19_susScr3.bed6+ > hg19_susScr3.tsv
###4.3.2 Select those exons that are lifted as a whole and stastics their depth and junction reads
mkdir rheMac2 tupBel1 mm9 susScr3
cd rheMac2
mkdir genomeDepth
awk '$4==$14 && $13=="1/1"' ../hg19_rheMac2.bed6+ >tmp.bed6+
ls ~/nucleosome/assignExonAge/rheMac2/rhesus*/cov/*.bw |while read file;do
        indiv=$(echo $file|cut -f8 -d '/');
        tissue=$(basename $file|sed 's/.bg.bw//'|sed 's/.uniq//')
        perl ~/bin/bigWigSummaryForBed.pl -bed tmp.bed6+ -bw $file >${indiv}_${tissue}.depth 2>err
        perl ~/bin/bigWigSumForGenome.pl -c /data/chr.size/rheMac2.size -bw $file >genomeDepth/${indiv}_${tissue}.total
done;
ls *.depth|while read file;do
        indiv=$(echo $file|cut -f1 -d '_')
        tissue=$(echo $file|cut -f2 -d '_'|cut -f1 -d '.')
        perl ~/nucleosome/scripts/sum_junctionReadsForExon.pl -e $file -j ~/nucleosome/assignExonAge/rheMac2/${indiv}/junc/${tissue}.bed >${indiv}_${tissue}.depth.junc
        done;
myAssign(){
    ls *.depth.junc|while read file;do
        prefix=$(echo $file|sed 's/.depth.junc//')
        total=$(cat genomeDepth/${prefix}.total|tr -d '\n');
        awk -v var=$total '{for(i=1;i<=14;i++){printf $i"\t"}print ($15/var)*1e9"\t"$16"\t"$17}' ${prefix}.depth.junc >${prefix}.dpkm.junc
        done;
    bedtools intersect -loj -s -a tmp.bed6+ -b ~/nucleosome/assignExonAge/species_exon/${species}/merged.uniq.exon.bed >intersect.tsv 
    awk -v tag=$tag '{if($3==$17 && $2==$16){for(i=1;i<=14;i++){printf $i"\t"}print tag"\tE"}}' intersect.tsv |sort|uniq >${species}.tag.annotated.bed12+ 
    comm -23 <(cut -f1-14 intersect.tsv|sort|uniq) <(awk '$16>0' intersect.tsv|cut -f1-14|sort|uniq)|awk '{print $0"\t-\tE"}' >>${species}.tag.annotated.bed12+ 
    awk '$15<=0.2 && $16==0 && $17==0' *.dpkm.junc|cut -f1-14|sort|uniq -c|awk -v value=$sample '$1==value'|awk '{for(i=2;i<=NF;i++){printf $i"\t"}print "-\t"$1}' >${species}.tag.bed12+
    awk '$16>=1 && $17>=1' *.dpkm.junc|cut -f1-14|sort|uniq -c|awk '$1>=2'|awk -v tag=$tag '{for(i=2;i<=NF;i++){printf $i"\t"}print tag"\t"$1}' >>${species}.tag.bed12+
    perl ~/nucleosome/scripts/annotatedAsExon.pl ${species}.tag.bed12+ ${species}.tag.annotated.bed12+ >${species}.tag.final.bed12+
}
species="rheMac2"
sample=17
tag="R"
myAssign
rm *depth *depth.junc

cd ../tupBel1
mkdir genomeDepth
awk '$4==$14 && $13=="1/1"' ../hg19_tupBel1.bed6+ >tmp.bed6+
ls ~/nucleosome/assignExonAge/tupBel1/*/*/cov/*.bw|while read file;do
        indiv=$(echo $file|cut -f8 -d '/');
        tissue=$(echo $file|cut -f9 -d '/');
        perl ~/bin/bigWigSummaryForBed.pl -bed tmp.bed6+ -bw $file >${indiv}.${tissue}.depth 2>err
        perl ~/bin/bigWigSumForGenome.pl -c /data/chr.size/tupBel1.size -bw $file >genomeDepth/${indiv}.${tissue}.total
done;
ls *.depth|while read file;do
        indiv=$(echo $file|cut -f1 -d '.')
        tissue=$(echo $file|cut -f2 -d '.')
        perl ~/nucleosome/scripts/sum_junctionReadsForExon.pl -e $file -j ~/nucleosome/assignExonAge/tupBel1/${indiv}/${tissue}/junc/junc.bed12 >${indiv}.${tissue}.depth.junc
done;
species="tupBel1"
sample=12
tag="T"

cd ../mm9
mkdir genomeDepth
awk '$4==$14 && $13=="1/1"' ../hg19_mm9.bed6+ >tmp.bed6+
ls ~/nucleosome/assignExonAge/mm9/*/*/*.bw|while read file;do
        indiv=$(echo $file|cut -f8 -d '/');
        tissue=$(echo $file|cut -f9 -d '/');
        perl ~/bin/bigWigSummaryForBed.pl -bed tmp.bed6+ -bw $file >${indiv}.${tissue}.depth 2>err
        perl ~/bin/bigWigSumForGenome.pl -c /data/chr.size/mm9.size -bw $file >genomeDepth/${indiv}.${tissue}.total
done;
ls *.depth|while read file;do
        indiv=$(echo $file|cut -f1 -d '.')
        tissue=$(echo $file|cut -f2 -d '.')
        perl ~/nucleosome/scripts/sum_junctionReadsForExon.pl -e $file -j ~/nucleosome/assignExonAge/mm9/${indiv}/${tissue}/*junc* >${indiv}.${tissue}.depth.junc
 done;
species="mm9"
sample=13
tag="M"


cd ../susScr3
mkdir genomeDepth
awk '$4==$14 && $13=="1/1"' ../hg19_susScr3.bed6+ >tmp.bed6+
ls ~/nucleosome/assignExonAge/susScr3/*/*/cov/*.bw|while read file;do
        indiv=$(echo $file|cut -f8 -d '/');
        tissue=$(echo $file|cut -f9 -d '/');
        perl ~/bin/bigWigSummaryForBed.pl -bed tmp.bed6+ -bw $file >${indiv}.${tissue}.depth 2>err
        perl ~/bin/bigWigSumForGenome.pl -c /data/chr.size/susScr3.size -bw $file >genomeDepth/${indiv}.${tissue}.total
done;
ls *.depth|while read file;do
        indiv=$(echo $file|cut -f1 -d '.')
        tissue=$(echo $file|cut -f2 -d '.')
        perl ~/nucleosome/scripts/sum_junctionReadsForExon.pl -e $file -j ~/nucleosome/assignExonAge/susScr3/${indiv}/${tissue}/junc/junc.bed12 >${indiv}.${tissue}.depth.junc
done;
species="susScr3"
sample=17
tag="S"
###4.3.3 Assign ages to huaman exons
cd ~/nucleosome/assignExonAge/list_v3/pairwiseAln
perl ~/nucleosome/scripts/assignExonAge.pl -hg hg19.initial.exon.bed6 -o rheMac2/rheMac2.tag.final.bed12+,tupBel1/tupBel1.tag.final.bed12+,mm9/mm9.tag.final.bed12+,susScr3/susScr3.tag.final.bed12+ >exon_age.tsv
awk '$7=="H----"' exon_age.tsv >hg19.exonAge.bed6+ 
ls ../*exonAge.bed6+|grep -v "hg19"|while read file;do awk '$7=="H----"' $file|cut -f1-3,6|sort|uniq -c|awk -v OFS="\t" '$1>1{print $2,$3,$4,"","0",$5}'|bedtools intersect -s -f 1.0 -r -a $file -b stdin |cut -f4 >>potential.Dup.list;done;
sort potential.Dup.list |uniq tmp;mv tmp potential.Dup.list;
awk '$7=="H----"' ../hg19.exonAge.bed6+|sort -k4,4|join -v 1 -t $'\t' -1 4 -2 1 - potential.Dup.list |grep -v "random"|awk -v OFS="\t" '{print $2,$3,$4,$1,$5,$6}' >hg19.H.bed6+ 
sed 's/chr//' hg19.H.bed6+ >hg19.H.noChr.bed6
cat hg19.H.noChr.bed6 >hg19.H.startJunc.tsv;
cat hg19.H.noChr.bed6 >hg19.H.endJunc.tsv;
ls /rd1/user/liym/nucleosome/assignExonAge/GTEx_hg19/junc_all/*.bed12 >hg19.juncFileName.txt  
ls /rd1/user/liym/nucleosome/assignExonAge/GTEx_hg19/junc/inhouse_* >>hg19.juncFileName.txt 
cat hg19.juncFileName.txt|while read file;do
    perl /mnt/share/liym/bin/sum_junctionReadsForExon.pl -e hg19.H.noChr.bed6 -j $file >tmp;
    paste hg19.H.startJunc.tsv <(cut -f7 tmp) >tmp3.tsv; 
    mv tmp3.tsv hg19.H.startJunc.tsv;
    paste hg19.H.endJunc.tsv <(cut -f8 tmp) >tmp5.tsv;
    mv tmp5.tsv hg19.H.endJunc.tsv;
done;
cat <(awk '{sum=0;for(i=7;i<=NF;i++){sum+=$i}if(sum>0){print $4}}' hg19.H.startJunc.tsv) <(awk '{sum=0;for(i=7;i<=NF;i++){sum+=$i}if(sum>0){print $4}}' hg19.H.endJunc.tsv)|sort|uniq|awk '{print "chr"$0}' >GTEx.supported.list
join.pl -i1 ../../species_exon/hg19/hg19.all.internal.exon.bed6 -f1 4 -i2 GTEx.supported.list -o1 |cut -f4 >final.list
join.pl -i1 ../../species_exon/hg19/hg19.all.internal.exon.bed6 -f1 4 -i2 final.list -o1 >hg19.H.bed6
join.pl -i1 hg19_rheMac2.bed6+ -f1 4 -i2 final.list -o1 |cut -f1-6 >rheMac2.H.bed6
join.pl -i1 hg19_tupBel1.bed6+ -f1 4 -i2 final.list -o1 |cut -f1-6 >tupBel1.H.bed6
join.pl -i1 hg19_mm9.bed6+ -f1 4 -i2 final.list -o1 |cut -f1-6 >mm9.H.bed6
join.pl -i1 hg19_susScr3.bed6+ -f1 4 -i2 final.list -o1 |cut -f1-6 >susScr3.H.bed6
#5 Nucleosome mediated exon evolution @pluto
cd /rd1/user/liym/nucleosome/mutiSpecies/exonAge
## 5.1 NC log ratio & GC log ratio & splice score
scp liym@jupiter:~/nucleosome/assignExonAge/list_v2/pairwiseAln/*H.bed6 .
mkdir NCratio && cd NCratio
ls ../*.bed6|while read file;do 
        species=$(basename $file|sed 's/.H.bed6//')
        ls /rd1/user/liym/nucleosome/mutiSpecies/data/${species}/*.norm.bw|while read ncFile;do
                tissue=$(basename $ncFile|sed 's/.norm.bw//')
                perl /mnt/share/liym/nucleosome/scripts/NC_log_ratio.pl -b $file -w $ncFile -u 150 -d 150 -p 0.00001 >${species}/${tissue}.logNC.tsv 2>logNC.log;
        done;
done;
paste *brain.logNC.tsv|cut -f4,7,15,23,31,39 |awk -v OFS="\t" '{print $1,$2,$4,$6,$3,$5}' >brain.5tissues.final.tsv 
paste *brain.logNC.tsv|cut -f4,8,16,24,32,40 |awk -v OFS="\t" '{print $1,$2,$4,$6,$3,$5}' >brain.5tissues.final.downLogNC.tsv 
mkdir GCratio && cd GCratio
ls ../*H.bed6|while read file;do
        species=$(basename $file|sed 's/.H.bed6//');
        perl /mnt/share/liym/nucleosome/scripts/NC_log_ratio.pl -b $file -w /mnt/share/liym/data/gcContent/${species}.gc5Base.bw -u 150 -d 150 >GCratio/${species}.logGC.tsv 2>GCratio/${species}.logGC.log
        done;
mkdir ss && cd ss
ls ../*H.bed6|while read file;do
        species=$(basename $file|sed 's/.H.bed6//');
        perl /mnt/share/liym/nucleosome/scripts/spliceSite_score.pl -b $file -f /mnt/share/liym/data/genome/${species}/${species}.fa >ss/${species}.ss.bed6+ 2>ss/${species}.ss.log
        done;
paste <(cut -f4,7,8 ss/hg19.ss.bed6+) <(cut -f8 ss/rheMac2.ss.bed6+) <(cut -f8 ss/tupBel1.ss.bed6+) <(cut -f8 ss/mm9.ss.bed6+) <(cut -f8 ss/susScr3.ss.bed6+) >H----.summary.5ss.tsv 
paste <(cut -f4,7,9 ss/hg19.ss.bed6+) <(cut -f9 ss/rheMac2.ss.bed6+) <(cut -f9 ss/tupBel1.ss.bed6+) <(cut -f9 ss/mm9.ss.bed6+) <(cut -f9 ss/susScr3.ss.bed6+) >H----.summary.3ss.tsv 
paste <(cut -f4,7,8 GCratio/hg19.logGC.tsv) <(cut -f8 GCratio/rheMac2.logGC.tsv) <(cut -f8 GCratio/tupBel1.logGC.tsv) <(cut -f8 GCratio/mm9.logGC.tsv) <(cut -f8 GCratio/susScr3.logGC.tsv) >H----.summary.uplogGC.tsv
paste <(cut -f4,7,9 GCratio/hg19.logGC.tsv) <(cut -f9 GCratio/rheMac2.logGC.tsv) <(cut -f9 GCratio/tupBel1.logGC.tsv) <(cut -f9 GCratio/mm9.logGC.tsv) <(cut -f9 GCratio/susScr3.logGC.tsv) |join.pl -i2 final.list -o1 >H----.summary.downlogGC.final.tsv

## 5.2 splice motif
### 5.2.1 Assign common ancestor sequence
cd /rd1/user/liym/data
cat /rd1/brick/share/data/conservation/multiz/hg19/*.maf|mafSpeciesSubset stdin species.list hg19.rheMac2.calJac1.maf
perl /mnt/share/liym/nucleosome/scripts/commonAncestorFromMaf_3species.pl -m hg19.rheMac2.calJac1.maf -t hg19 -q rheMac2 -og calJac1 --o1 hg19.ancestor.3species.fa --o2 rheMac2.ancestor.3species.fa
cat /rd1/brick/share/data/conservation/multiz/hg19/*.maf|mafSpeciesSubset stdin species.list.5species 5species.maf
perl /mnt/share/liym/bin/commonAncestorFromMaf.pl -m 5species.maf -t hg19 -q rheMac2,mm9,monDom5,bosTau4 -o hg19.ancestor.5species.fa
### 5.2.2 WebLogo
cd /rd1/user/liym/nucleosome/mutiSpecies/exonAge/final_list/ss_motif
ls ../*H.bed6|while read file;do
        prefix=$(basename $file|sed 's/.H.bed6//');
        awk -v OFS="\t" '{if($6=="+"){print $1,$3-3,$3+20,$4,$5,$6}else{print $1,$2-20,$2+3,$4,$5,$6}}' $file >${prefix}.3ss.bed6 ##20 bases in intron,3 bases in exon
done;
fastaFromBed -name -s -fi /mnt/share/liym/data/genome/mm9/mm9.fa -bed mm9.H----.3ss.bed6 -fo H----.mm9.3ss.fa
fastaFromBed -name -s -fi /mnt/share/liym/data/genome/hg19/hg19.fa -bed hg19.H----.3ss.bed6 -fo H----.hg19.3ss.fa
fastaFromBed -name -s -fi /rd1/user/liym/data/hg19.ancestor.5species.fa -bed hg19.H----.3ss.bed6 -fo H----.hg19.ancestor.3ss.fa 
fastaFromBed -name -s -fi /mnt/share/liym/data/genome/rheMac2/rheMac2.fa -bed rheMac2.H----.3ss.bed6 -fo H----.rheMac2.3ss.fa
perl -ne '
        open IN,"$ARGV[0]";
        my $name;
        while(chomp(my $line=<IN>)){
                if($line=~/^>/){
                        $line=~s/>//;   
                        $name=$line;
                }else{
                        if($line=~/[ATGC]/){
                                print "$name\n";
                        }
                }
        }
' H----.hg19.ancestor.3ss.fa >exonWithAncestor.txt
comm -12 <(sort exonWithAncestor.txt) <(sort ../final.list)|join.pl -i1 mm9.H----.3ss.bed6 -f1 4 -o1|fastaFromBed -name -s -fi /mnt/share/liym/data/genome/mm9/mm9.fa -bed stdin -fo H----.mm9.3ss.filtered.fa;
comm -12 <(sort exonWithAncestor.txt) <(sort ../final.list)|join.pl -i1 hg19.H----.3ss.bed6 -f1 4 -o1|fastaFromBed -name -s -fi /mnt/share/liym/data/genome/hg19/hg19.fa -bed stdin -fo H----.hg19.3ss.filtered.fa;
comm -12 <(sort exonWithAncestor.txt) <(sort ../final.list)|join.pl -i1 rheMac2.H----.3ss.bed6 -f1 4 -o1|fastaFromBed -name -s -fi /mnt/share/liym/data/genome/rheMac2/rheMac2.fa -bed stdin -fo H----.rheMac2.3ss.filtered.fa
prefix=H----.mm9.3ss.filtered;
python /mnt/share/liym/tools/weblogo-master/weblogo -f ${prefix}.fa -D fasta -o ${prefix}.filtered.pdf -F pdf -A dna -i -20 -l -20 -u 2 -t ancestor --number-interval 1 --fineprint "" --errorbars False -C green A 'Adenine' -C orange G 'guanine' -C red T 'Thymine' -C blue C 'Cytosine'
prefix=H----.hg19.3ss.filtered;
prefix=H----.rheMac2.3ss.filtered;

## 5.3 Nucleosome's influence to GC content
cd /rd1/user/liym/nucleosome/mutiSpecies/GC_evolution/hg19_regions
perl /mnt/share/liym/bin/bedToIntervals.pl -b /mnt/share/liym/data/structure/hg19/hg19.refGene.intergenic.bed -l 100 -d 100| awk -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' >hg19.intergenic.100bp.bed
liftOver -bedPlus=3 ../hg19.intergenic.100bp.bed /data/liftover/hg19/hg19ToRheMac2.over.chain rheMac2.intergenic.bed rheMac2.intergenic.unmapped &
cd liftOver/hg19_rheMac2
join.pl -i1 ../../hg19.intergenic.100bp.bed -f1 4 -i2 ../rheMac2.intergenic.bed -f2 4 | awk '$7-$6==100' > HR.intergenic.bed
ls --color=auto /rd1/user/liym/nucleosome/mutiSpecies/data/hg19/*.bw | while read file; do
    tissue=$(basename $file|sed 's/.norm.bw//'); cut -f1-4 HR.intergenic.bed | bwtool summary -skip-median -keep-bed stdin $file intergenic.hg19.${tissue}.NC.tsv & cut -f1-4 HR.intron.bed | bwtool summary -skip-median -keep-bed stdin $file intron.hg19.${tissue}.NC.tsv;
        done &
ls --color=auto /rd1/user/liym/nucleosome/mutiSpecies/data/rheMac2/*.bw | while read file; do
    tissue=$(basename $file|sed 's/.norm.bw//'); cut -f5-8 HR.intergenic.bed | bwtool summary -skip-median -keep-bed stdin $file intergenic.rheMac2.${tissue}.NC.tsv & cut -f5-8 HR.intron.bed | bwtool summary -skip-median -keep-bed stdin $file intron.rheMac2.${tissue}.NC.tsv;
        done &
paste <(cut -f1-4,9 intergenic.rheMac2.brain.NC.tsv) <(cut -f9 intergenic.rheMac2.heart.NC.tsv) <(cut -f9 intergenic.rheMac2.kidney.NC.tsv) <(cut -f9 intergenic.rheMac2.liver.NC.tsv) <(cut -f9 intergenic.rheMac2.muscle.NC.tsv) >intergenic.rheMac2.bhklm.NC.tsv
paste <(cut -f1-4,9 intergenic.hg19.brain.NC.tsv) <(cut -f9 intergenic.hg19.muscle.NC.tsv) >intergenic.hg19.bm.NC.tsv 
ls *hg19.bm.NC.tsv|sed 's/.tsv//'|while read file;do awk -v OFS="\t" '{print $1,$2,$3,$4,($5+$6)/2}' ${file}.tsv >${file}.mean.tsv;done;
ls *rheMac2.bhklm.NC.tsv|sed 's/.tsv//'|while read file;do awk -v OFS="\t" '{print $1,$2,$3,$4,($5+$6+$7+$8+$9)/5}' ${file}.tsv >${file}.mean.tsv;done;
paste <(awk -v OFS="\t" '{print $1,$2,$3,$4,($5+$6)/2}' intergenic.hg19.bm.NC.tsv) <(awk -v OFS="\t" '{print $1,$2,$3,$4,($5+$6+$7+$8+$9)/5}' intergenic.rheMac2.bhklm.NC.tsv) >intergenic.hg19.rheMac2.tissueMean.NC.tsv
cut -f1-4 HR.intergenic.bed | perl /mnt/share/liym/bin/GC_content.pl -f /rd1/user/liym/data/hg19.ancestor.3species.fa > intergenic.ancestor.GC.tsv
awk '$5>2 && $10>2' intergenic.hg19.rheMac2.tissueMean.NC.tsv |cut -f1-4|join.pl -i1 intergenic.ancestor.hg19.GC.tsv -f1 4 -f2 4 -o1 >intergenic.class1.anestor.hg19.GC.tsv &
awk '$5>2 && $10<0.5' intergenic.hg19.rheMac2.tissueMean.NC.tsv |cut -f1-4|join.pl -i1 intergenic.ancestor.hg19.GC.tsv -f1 4 -f2 4 -o1 >intergenic.class2.anestor.hg19.GC.tsv &
awk '$5<0.5 && $10>2' intergenic.hg19.rheMac2.tissueMean.NC.tsv|cut -f1-4 | join.pl -i1 intergenic.ancestor.hg19.GC.tsv -f1 4 -f2 4 -o1 > intergenic.class3.anestor.hg19.GC.tsv &
awk '$5<0.5 && $10<0.5' intergenic.hg19.rheMac2.tissueMean.NC.tsv|cut -f1-4 | join.pl -i1 intergenic.ancestor.hg19.GC.tsv -f1 4 -f2 4 -o1 > intergenic.class4.anestor.hg19.GC.tsv &
###Control ancestor GC content
grep -v "NA" ../intergenic.class4.anestor.hg19.GC.tsv >intergenic.class4.ancestor.hg19.GC.tsv
grep -v "NA" ../intergenic.class1.anestor.hg19.GC.tsv >intergenic.class1.ancestor.hg19.GC.tsv
perl /mnt/share/liym/bin/getComparableDataList.pl -i1 intergenic.class4.ancestor.hg19.GC.tsv -f1 5 -i2 intergenic.class1.ancestor.hg19.GC.tsv -f2 5 -d 0.003 >intergenic.class4.class1.comparible.ancestor.newList
cut -f1-4 intergenic.class4.class1.comparible.ancestor.newList | perl /mnt/share/liym/bin/divRate_pooled.pl -a /rd1/user/liym/data/hg19.ancestor.3species.fa -f /mnt/share/liym/data/genome/hg19/hg19.fa > intergenic.class4.divRate.new.tsv
cut -f7-10 intergenic.class4.class1.comparible.ancestor.newList | perl /mnt/share/liym/bin/divRate_pooled.pl -a /rd1/user/liym/data/hg19.ancestor.3species.fa -f /mnt/share/liym/data/genome/hg19/hg19.fa > intergenic.class1.divRate.new.tsv
mkdir bootstrap
cut -f1-4 intergenic.class4.class1.comparible.ancestor.newList|perl /mnt/share/liym/bin/region_divRate.pl -a /rd1/user/liym/data/hg19.ancestor.3species.fa -f /mnt/share/liym/data/genome/hg19/hg19.fa >bootstrap/intergenic.class4.regionDiv.tsv
cut -f7-10 intergenic.class4.class1.comparible.ancestor.newList|perl /mnt/share/liym/bin/region_divRate.pl -a /rd1/user/liym/data/hg19.ancestor.3species.fa -f /mnt/share/liym/data/genome/hg19/hg19.fa >bootstrap/intergenic.class1.regionDiv.tsv
sh /mnt/share/liym/nucleosome/scripts/bootstrap_pooledDivRate.v2.sh -i bootstrap/intergenic.class1.regionDiv.tsv -n 1000 -o bootstrap/intergenic.class1.divRate.boot.tsv -t bootstrap/intergenic.class1
sh /mnt/share/liym/nucleosome/scripts/bootstrap_pooledDivRate.v2.sh -i bootstrap/intergenic.class4.regionDiv.tsv -n 1000 -o bootstrap/intergenic.class4.divRate.boot.tsv -t bootstrap/intergenic.class4
mkdir rheMac2
join.pl -i1 ../../HR.intergenic.bed -f1 4 -i2 ../intergenic.class4.class1.comparible.ancestor.newList -f2 4 -o1|cut -f5-8 >intergenic.class4.rheMac2.bed4; 
cut -f7-10 ../intergenic.class4.class1.comparible.ancestor.newList|join.pl -i1 ../../HR.intergenic.bed -f1 4 -f2 4 -o1 |cut -f5-8 >intergenic.class1.rheMac2.bed4;
perl /mnt/share/liym/bin/region_divRate.pl -b intergenic.class1.rheMac2.bed4 -a /rd1/user/liym/data/rheMac2.ancestor.3species.fa -f /mnt/share/liym/data/genome/rheMac2/rheMac2.fa >intergenic.class1.rheMac2.regionDiv2.tsv;
perl /mnt/share/liym/bin/region_divRate.pl -b intergenic.class4.rheMac2.bed4 -a /rd1/user/liym/data/rheMac2.ancestor.3species.fa -f /mnt/share/liym/data/genome/rheMac2/rheMac2.fa >intergenic.class4.rheMac2.regionDiv2.tsv;
sh /mnt/share/liym/nucleosome/scripts/bootstrap_pooledDivRate.v2.sh -i intergenic.class1.rheMac2.regionDiv2.tsv -n 1000 -o intergenic.class1.rheMac2.divRate.boot2.tsv -t intergenic.class1;
sh /mnt/share/liym/nucleosome/scripts/bootstrap_pooledDivRate.v2.sh -i intergenic.class4.rheMac2.regionDiv2.tsv -n 1000 -o intergenic.class4.rheMac2.divRate.boot2.tsv -t intergenic.class4;

