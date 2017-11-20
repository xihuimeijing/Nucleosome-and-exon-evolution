#!/usr/bin/env bash

usage(){
  cat << EOF
Description:
    The tool wraps steps of sra decompressing, fastq quality control, fastq pre-processing, mapping, a set of post-mapping evaluations and expression level calculating:
    
    Step		        Program
    sra decompressing:		fastq-dump
    fastq quality control:	fastqc
    fastq pre-processing:       fqTrimer.pl, fqPeFilter.pl
    mapping:			For tophat1.2: tophat1.2, bowtie0.12.3 or later (not 1.x), bowtie-build
                                For tophat 2: tophat2, bowtie2, bowtie2-align and bowtie2-build
    post-mapping evaluations:	samtools, bamtools, Picard (only needed by tophat2), mutationRate_tophat.pl, bedtools,
                                samSsEva.pl (only needed by strand-specific sequencing),
                                samFlag.pl (only needed by strand-specific sequencing), coverage.R, mutRate_tophat.R,
                                gpeFeature.pl, gpeMerge.pl, gpe2bed.pl, fragmentation.py, fragmentation.R and samContam.pl
    calculate RB Score:         calRnaSeqRbScore.pl
    calculate expression level: geneRPKM.pl, transRPKM.pl, geneFPKM.pl
    find junction:              sam2junc.pl
    calculate coverage          wiggles
    
    Version and Modules:
        for perl: >=5.12, Sky's pm, Bio::DB::Fasta
        for python: =2.7, pysam, argparse
    
    So, please make sure above programs are in your PATH.
    Before running this tool, you must export the minimal anchor length (minAnchor) and the arguments needed by tophat (tophatArgs).
    This tool decompresses (unless -d specified) SRA files under -i directory to fastq files under -o directory,
    Pool all fastq/fq files under -i and -o directory to raw1.fq and raw2.fq under -o directory.
    If fqTrimerArgs1/2 specified, trim raw1/2.fq and store into trim1/2.fq.
    If fqFilterArgs specified, filter above fqs and store into flt1/2.fq.
    Pre-processed fqs are then fed to tophat for mapping.

Usage:
    export fqDumpArgs='All fastq-dump arguments' (Optional. Ignored when -d option is specified. Note: use single quotes here, not double quotes!)
    export fqTrimerArgs1='All fqTrimer.pl arguments used to trim read1' (Optional. No triming when not specified)
    export fqTrimerArgs2='All fqTrimer.pl arguments used to trim read2' (Optional. No triming when not specified)
    export fqFilterArgs='All fqPeFilter.pl arguments' (Optional. No filtering when not specified)
    export minAnchor='The minimal anchor length (read length/10 is recommended)'
    export tophatArgs="All Tophat2 arguments (with index prefix). --out-dir argument will be set as tophat_out under -o directory. --library-type will be set by -l argument of this shell script. --min-anchor will be set as the exported minAnchor in the previous command"
    Note: above arguments must be exported into $(basename $0) using export command
    $(basename $0) Options >$(basename $0 .sh).log 2>$(basename $0 .sh).err
	
Example:
    export fqDumpArgs='--split-3 --defline-seq @\$ac.\$si/\$ri --defline-qual +' (Note: use single quotes here, not double quotes!)
    export fqTrimerArgs1='-l G+'
    export fqTrimerArgs2='-r A+'
    export fqFilterArgs='-f 0.1 -b N -q 0.5 -c 5 -p J'
    export minAnchor=4
    # for tophat1.2
    export tophatArgs="--allow-indels --mate-inner 270 --mate-std-dev 75 --microexon-search --coverage-search --num-threads 5 \
--segment-length 20 --rg-id SEX196335 --rg-sample SRS369772_brain --rg-center Burge --rg-platform HiSeq2000 \
/data/genome/bowtie1/rheMac2/rheMac2"
    # for tophat2
    export tophatArgs="--mate-inner 270 --mate-std-dev 75 --microexon-search --coverage-search --num-threads 5 \
--segment-length 20 --rg-id SEX196335 --rg-sample SRS369772_brain --rg-center Burge --rg-platform HiSeq2000 \
/data/genome/bowtie1/rheMac2/rheMac2"
    $(basename $0) -f /data/genome/fna/rheMac2/all.fa -g /data/genome/structure/gpe/rheMac2.rb2Gene.gpe >$(basename $0 .sh).log 2>$(basename $0 .sh).err

Options:
    -t  1.2/2   The tophat version (1.2 or 2)
    -i  DIR     Input directory containing all fq/fastq files (with _1.fq/fastq and _2.fq/fastq suffix)
                or SRA files with sra suffix[Current Directory]
    -o  DIR     Output directory for all results[Current Directory]
    -d          Do not run fastq-dump for sra files under directory specified by -i
    -e          Only run evaluation and its following steps
    -m          Don't run evaluation
    -f  FILE    Reference sequence in fasta format
    -g  FILE    A gene structure file in gpe format
    -s  INT     The slopping length from the exon-intron joint to intron[Half of minAnchor]
    -l  STR     Sequencing library type, it can be
                    fr-unstranded: for Standard Illumina (default)
                    fr-firststrand: for dUTP, NSR, NNSR
                    fr-secondstrand: for Ligation, Standard SOLiD and Illumina Directional Protocol
    -r  DOU     The RNA Integrity Number[0]
    -p  DIR     The path where Picard jars (CollectInsertSizeMetrics.jar) located[/data/tool/picard]
    -h          Show this help information
EOF
    exit 0
}

[ $1 ] || usage

echo "--------------- Begin at $(date) ---------------"

inputDir=$PWD
outputDir=$PWD
picardPath=/data/tool/picard
libType=fr-unstranded
RIN=0
ssRate=0
perlPath=/mnt/share/liym/bin/rnaSeq_zhangsj

while getopts "ht:i:o:demf:g:r:p:s:l:" OPTION
do
    case $OPTION in
        h) usage;;
        t) tophatVer=$OPTARG;;
        i) inputDir=$OPTARG;;
        o) outputDir=$OPTARG;;
        d) noDecompress=1;;
        e) evaluationOnly=1;;
        m) noEvaluation=1;;
        f) refFa=$OPTARG;;
        g) gpeStructure=$OPTARG;;
        s) slop=$OPTARG;;
        l) libType=$OPTARG;;
        r) RIN=$OPTARG;;
        p) picardPath=$OPTARG;;
        ?) usage;;
    esac
done
shift $((OPTIND - 1))

if [ -z "$slop" ];then
  if [ -z "$minAnchor" ];then
    echo 'Please specify -l or export minAnchor' >&2
    exit 1
  else
    slop=$(($minAnchor/2))
  fi
fi

if [ "X$libType" != "Xfr-unstranded" ] && [ "X$libType" != "Xfr-firststrand" ] && [ "X$libType" != "Xfr-secondstrand" ];then
  echo "Please specify correct libType by -l" >&2
  exit 1
fi

inputDir=$(readlink -f $inputDir) # get full path
outputDir=$(readlink -f $outputDir)
inputDir=$(readlink -f $inputDir) # get full path
outDir=$(readlink -f $outputDir)
if [ ! -e $refFa ];then
  echo "The Reference sequences file $refFa dosen't exist" >&2
  exit 1
fi
refFa=$(readlink -f $refFa)
if [ ! -e $gpeStructure ];then
  echo "The gene structure file file $gpeStructure dosen't exist" >&2
  exit 1
fi
gpeStructure=$(readlink -f $gpeStructure)
picardPath=$(readlink -f $picardPath)

mkdir -p $outputDir
if [ $? -ne 0 ];then
  echo "Can't create directory $outputDir" >&2
  exit 1
fi
cd $outputDir

if [ -d log ];then
  [ -d log_old ] && rm -rf log_old
  mv log log_old
fi
mkdir log

if [ -z $evaluationOnly ];then
  if [ $(ls $inputDir | grep '\.sra$' | wc -l) -gt 0 ] && [ -z $noDecompress ];then
    echo "Uncompressing SRA File... at $(date)"
    ls $inputDir/*.sra | while read file;do
      fastq-dump $fqDumpArgs $file	# -Q 2 can output as sanger format
    done
    if [ $? -ne 0 ];then
      echo 'Something wrong with running fastq-dump.' >&2
      exit 1
    fi
  fi
  if [ $(ls $inputDir | grep -P '_1.*\.f(ast)?q$' | wc -l) -eq 0 ] && [ $(ls $inputDir | grep -P '_2.*\.f(ast)?q$' | wc -l) -eq 0 ] && \
  [ $(ls | grep -P '_1.*\.f(ast)?q$' | wc -l) -eq 0 ] && [ $(ls | grep -P '_2.*\.f(ast)?q$' | wc -l) -eq 0 ];then
    echo "No fa/fastq file(s) in '$inputDir' and '$outputDir'. Exit" >&2
    exit 1
  fi
  
  [ -e raw1.fq ] && rm raw1.fq
  [ -e raw2.fq ] && rm raw2.fq
  [ -e trim1.fq ] && rm trim1.fq
  [ -e trim2.fq ] && rm trim2.fq
  [ -e flt1.fq ] && rm flt1.fq
  [ -e flt2.fq ] && rm flt2.fq
  [ $(ls | grep '^final1\.fq$' | wc -l) -gt 0 ] && rm final1.fq
  [ $(ls | grep '^final2\.fq$' | wc -l) -gt 0 ] && rm final2.fq
  
  if [ $(ls $inputDir | grep '_[12].*\.fastq$' | wc -l) -gt 0 ];then
    if [ $(ls $inputDir | grep '_1.*\.fastq$' | wc -l) -eq $(ls $inputDir | grep '_2.*\.fastq$' | wc -l) ];then
      cat $inputDir/*_1*.fastq >raw1.fq
      cat $inputDir/*_2*.fastq >raw2.fq
    else
      echo "Unpaired read 1 and 2 fastq files in '$inputDir'. Exit" >&2
      exit 1
    fi
  fi
  if [ $(ls $inputDir | grep '_[12].*\.fq$' | wc -l) -gt 0 ];then
    if [ $(ls $inputDir | grep '_1.*\.fq$' | wc -l) -eq $(ls $inputDir | grep '_2.*\.fq$' | wc -l) ];then
      cat $inputDir/*_1*.fq >>raw1.fq
      cat $inputDir/*_2*.fq >>raw2.fq
    else
      echo "Unpaired read 1 and 2 fq files in '$inputDir'. Exit" >&2
      exit 1
    fi
  fi
  if [ "X$inputDir" != "X$outputDir" ];then
    if [ $(ls | grep '_[12].*\.fastq$' | wc -l) -gt 0 ];then
      if [ $(ls | grep '_1.*\.fastq$' | wc -l) -eq $(ls | grep '_2.*\.fastq$' | wc -l) ];then
        cat *_1*.fastq >>raw1.fq
        cat *_2*.fastq >>raw2.fq
      else
        echo "Unpaired read 1 and 2 fastq files in $PWD. Exit" >&2
        exit 1
      fi
    fi
    if [ $(ls | grep '_[12].*\.fq$' | wc -l) -gt 0 ];then
      if [ $(ls | grep '_1.*\.fq$' | wc -l) -eq $(ls | grep '_2.*\.fq$' | wc -l) ];then
        cat *_1*.fq >>raw1.fq
        cat *_2*.fq >>raw2.fq
      else
        echo "Unpaired read 1 and 2 fq files in $(pwd). Exit" >&2
        exit 1
      fi
    fi
  fi
  
  currentFq1='raw1.fq'
  currentFq2='raw2.fq'
  echo 'Run FastQC for Raw Fastq (Background Runing)'
  fastqc -q -t 2 --extract $currentFq1 $currentFq2 &
  [ $? -eq 0 ] || exit 1
  
  if [ "$fqTrimerArgs1" ];then
    fqTrimer.pl $fqTrimerArgs1 $currentFq1 >trim1.fq
    currentFq1='trim1.fq'
    echo 'Run FastQC for Trimmed Fastq1 (Background Runing)'
    fastqc -q --extract $currentFq1 &
  fi
  [ $? -eq 0 ] || exit 1
  
  if [ "$fqTrimerArgs2" ];then
    fqTrimer.pl $fqTrimerArgs2 $currentFq2 >trim2.fq
    currentFq2='trim2.fq'
    echo 'Run FastQC for Trimmed fastq2 (Background Runing)'
    fastqc -q --extract $currentFq2 &
  fi
  [ $? -eq 0 ] || exit 1
  
  if [ "$fqFilterArgs" ];then
    echo "Running fqPeFilter.pl $fqFilterArgs --i1 $currentFq1 --i2 $currentFq2 --o1 flt1.fq --o2 flt2.fq 2>log/fqFilter.log..."
    fqPeFilter.pl $fqFilterArgs --i1 $currentFq1 --i2 $currentFq2 --o1 flt1.fq --o2 flt2.fq 2>log/fqFilter.log
    currentFq1='flt1.fq'
    currentFq2='flt2.fq'
    echo 'Run FastQC for Filtered Fastq (Background Runing)'
    fastqc -q -t 2 --extract $currentFq1 $currentFq2 &
  fi
  [ $? -eq 0 ] || exit 1
  
  ln -s $currentFq1 final1.fq
  ln -s $currentFq2 final2.fq
  [ $(ls | grep '^final1_fastqc$' | wc -l) -eq 0 ] && ln -s $(basename $currentFq1 .fq)_fastqc final1_fastqc
  [ $(ls | grep '^final2_fastqc$' | wc -l) -eq 0 ] && ln -s $(basename $currentFq2 .fq)_fastqc final2_fastqc
  
  if [ "X$tophatVer" == "X1.2" ];then
    echo "Running tophat1.2 -a $minAnchor --library-type $libType $tophatArgs final1.fq final2.fq 2>log/tophat_run.log... at $(date)"
    tophat1.2 -a $minAnchor --library-type $libType $tophatArgs final1.fq final2.fq 2>log/tophat_run.log
  else
    if [ "X$tophatVer" == "X2" ];then
      echo "Running tophat2 -a $minAnchor --library-type $libType $tophatArgs final1.fq final2.fq 2>log/tophat_run.log... at $(date)"
      tophat2 -a $minAnchor --library-type $libType $tophatArgs final1.fq final2.fq 2>log/tophat_run.log
    else
      echo 'Please specify the correct tophat version' >&2
      exit 1
    fi
  fi

  if [ $? -ne 0 ];then
    echo 'Something wrong with running tophat. Please view log/tophat_run.log for detail' 1>&2
    exit 1
  fi
  
  echo "Filtering Tophat Result... at $(date)"
  if [ "X$tophatVer" == "X1.2" ];then
    samtools view -bu -q 255 -@ 5 tophat_out/accepted_hits.bam | samtools sort -@ 5 -m 10G - uniq.sorted 2>log/samtoolsSort.log
  else
    # MAPQ (50: NH=1; 3: 1<=NH<=2; 1: 3<=NH<=4; 0: NH>=5)
    samtools view -bu -q 50 -@ 5 tophat_out/accepted_hits.bam | samtools sort -@ 5 -m 10G - uniq.sorted 2>log/samtoolsSort.log
  fi
  samtools index uniq.sorted.bam
  
  echo 'Run FastQC for uniq.sorted.bam (Background Runing)'
  fastqc -q --extract -f bam uniq.sorted.bam &
fi

if [ -z "$noEvaluation" ];then
    [ $? -eq 0 ] || exit 1
    if [ -d evaluation ];then
      [ -d evaluation_old ] && rm -r evaluation_old
      mv evaluation evaluation_old
    fi
    mkdir -p evaluation && cd evaluation
    
    echo 'Mapping Statistics (Background Running)'
    bamtools stats -in $outputDir/uniq.sorted.bam -insert >bamStats.txt &
    
    if [ "X$tophatVer" == "X2" ];then
      echo 'Insert Size Distribution Plot (Background Running)'
      java -jar $picardPath/CollectInsertSizeMetrics.jar H=insertSize.pdf I=$outputDir/uniq.sorted.bam O=insertSize.log 2>$outputDir/log/insertSize.log &
    fi
    
    if [ $(head -n 1 $gpeStructure | awk '{print NF}') -eq 15 ];then
      ln -s $gpeStructure myGpe.gpe
    else
      if [ $(head -n 1 $gpeStructure | awk '{print NF}') -eq 16 ];then
        cut -f2- $gpeStructure >myGpe.gpe
      else
        echo "Please offer the gene structure file in correct gpe format" >&2
        exit 1
      fi
    fi
    gpeStructure=$outputDir/evaluation/myGpe.gpe
    
    echo 'Fragmentation Evaluation (Background Running)'
    $perlPath/fragmentation.py --gpe $gpeStructure --bam $outputDir/uniq.sorted.bam >fragmentation.tsv && $perlPath/fragmentation.R <fragmentation.tsv &

    echo 'Mutation Rate Evaluation (Background Running)'
    $perlPath/mutationRate_tophat.pl --sam $outputDir/uniq.sorted.bam --ref $refFa >mutRate.tsv 2>$outputDir/log/mutRate.log && $perlPath/mutRate_tophat.R <mutRate.tsv &

    echo 'Coverage Evaluation (Background Running)'
    if [ "X$libType" == "Xfr-unstranded" ];then
      #-a and -b are distinct for different bedtools version
	  /mnt/share/liym/tools/bedtools2/bin/bedtools coverage -b $outputDir/uniq.sorted.bam -a <(gpeFeature.pl -e $gpeStructure | cut -f1-6) -hist -split >cov.tsv
    else
      if [ "X$libType" == "Xfr-secondstrand" ];then
        $perlPath/samFlag.pl -2 --rev $outputDir/uniq.sorted.bam | samtools view -buS - 2>/dev/null | /mnt/share/liym/tools/bedtools2/bin/bedtools coverage -s -b /dev/stdin -a <($perlPath/gpeFeature.pl -e $gpeStructure | cut -f1-6) -hist -split >cov.tsv
      else
        $perlPath/samFlag.pl -1 --rev $outputDir/uniq.sorted.bam | samtools view -buS - 2>/dev/null | /mnt/share/liym/tools/bedtools2/bin/bedtools coverage -s -b /dev/stdin -a <($perlPath/gpeFeature.pl -e $gpeStructure | cut -f1-6) -hist -split >cov.tsv
      fi
    fi && grep ^all cov.tsv | cut -f2,3,5 >cov_all.tsv && $perlPath/coverage.R <cov_all.tsv >cov_all_cum.tsv 
    
    echo 'Duplicate Evaluation'
    #java -jar $picardPath/picard.jar MarkDuplicates I=../uniq.sorted.bam O=uniq.sorted.rmDup.bam M=markDup.metrix.txt 2> markDup.err &
	#samtools rmdup $outputDir/uniq.sorted.bam uniq.rmdup.bam >/dev/null 2>rmdup.txt &
    #dupFraction=$(cut -f6 -d ' ' rmdup.txt)
        
    echo 'DNA Contamination Evaluation (Background Running)'
    if [ "X$libType" == "Xfr-unstranded" ];then
      $perlPath/samContam.pl -l $libType -g $gpeStructure -s $slop $outputDir/uniq.sorted.bam >contam.tsv 2>$outputDir/log/contam.log
    else
      $perlPath/samContam.pl -l $libType -g $gpeStructure -s $slop $outputDir/uniq.sorted.bam >contam.tsv 2>$outputDir/log/contam.log 
      echo 'Strand-specific Evaluation (Background Running)'
      $perlPath/gpeMerge.pl -n $gpeStructure | $perlPath/gpe2bed.pl >structure.bed
      bedtools subtract -a <(awk '$6=="+"' structure.bed) -b <(awk '$6=="-"' structure.bed) >noAntisenseRegion.bed
      bedtools subtract -a <(awk '$6=="-"' structure.bed) -b <(awk '$6=="+"' structure.bed) >>noAntisenseRegion.bed
      $perlPath/samSsEva.pl -b noAntisenseRegion.bed -l $libType $outputDir/uniq.sorted.bam >ssEva.tsv 2>$outputDir/log/ssEva.log && ssRate=$(grep "^Correct Strand-specific Reads:" ssEva.tsv | grep -oP "[0-9.]+%" | sed 's/%$//')
    fi && echo "Run calRnaSeqRbScore.pl -r $RIN -s $ssRate -u $(grep "^Total reads:" bamStats.txt | grep -oP "\d+") --rpkmExon2Intron $(sed -n '5p;6q' contam.tsv | cut -f 1) --rpkmExon2Intergenic $(sed -n '6p;7q' contam.tsv | cut -f 1) -d $(awk '$1==10{print $2}' cov_all_cum.tsv) --fragmentation fragmentation.tsv -m $(awk '$5<=0.01' mutRate.tsv | wc -l),$(awk '$9<=0.01' mutRate.tsv | wc -l) -f $outputDir/final2*_fastqc/fastqc_data.txt $outputDir/final1*_fastqc/fastqc_data.txt >rbScore.tsv 2>$outputDir/log/rbScore.err" && $perlPath/calRnaSeqRbScore.pl -r $RIN -s $ssRate -u $(grep "^Total reads:" bamStats.txt | grep -oP "\d+") --rpkmExon2Intron $(sed -n '5p;6q' contam.tsv | cut -f 1) --rpkmExon2Intergenic $(sed -n '6p;7q' contam.tsv | cut -f 1) -d $(awk '$1==10{print $2}' cov_all_cum.tsv) --fragmentation fragmentation.tsv -m $(awk '$5<=0.01' mutRate.tsv | wc -l),$(awk '$9<=0.01' mutRate.tsv | wc -l) -f $outputDir/final2*_fastqc/fastqc_data.txt $outputDir/final1*_fastqc/fastqc_data.txt >rbScore.tsv 2>$outputDir/log/rbScore.err &
    
    cd $outputDir
fi

#echo 'Find Junctions (Background Running)'
#mkdir junc
#sam2junc.pl -l $libType uniq.sorted.bam >junc/junctions.bed12 2>log/sam2junc.log &
##
#echo 'Calculate Coverage (Background Running)'
#mkdir cov
#samtools view uniq.sorted.bam | wiggles /dev/stdin cov/coverage.bg 2>/dev/null &

echo 'Calculate Expression Level...'
#$perlPath/geneRPKM.pl -l $libType -g $gpeStructure -s $slop uniq.sorted.bam >RPKM.gene.bed6+ 2>log/RPKM.gene.log & 
##transRPKM.pl -l $libType -g $gpeStructure -s $slop uniq.sorted.bam >RPKM.trans.bed6+ 2>log/RPKM.trans.log &
$perlPath/geneFPKM.pl -l $libType -g $gpeStructure -s $slop uniq.sorted.bam >FPKM.gene.bed6+ 2>log/FPKM.gene.log
#transFPKM.pl -l $libType -g $gpeStructure -s $slop uniq.sorted.bam >FPKM.trans.bed6+ 2>log/FPKM.trans.log &

echo "--------------- End at $(date) ---------------"
