#!/bin/sh
usage(){
  cat << EOF
Description: Statistics the difference between nucleosome occupancy of the target region and its flanking regions.
Usage: sh $0 -i <IN.bed> -w <NC.bw> -u <INT> -d <INT> -t <STRING> -o<OUT.tsv>
Author: Yumei Li,2017-7-11
Output: input_filelds targetRegion_NC up_NC down_NC
Options:
        i     FILE       The input file in bed6 format.
        w     FILE       The nucleosome occupancy data in bigWig format.
        u     INT        The upstream length of the target region.
	d     INT        The downstream length of the target region.
        t     STRING     The prefix of tmp file.
        o     FILE       The output file name.
        h                print this help message    
EOF
    exit 0
}
[ $1 ] || usage

while getopts "hi:w:u:d:t:o:" OPTION
do
    case $OPTION in
        h) usage;;
        i) inFile=$OPTARG;;
        w) bwFile=$OPTARG;;
        u) up=$OPTARG;;
        d) down=$OPTARG;;
        t) prefix=$OPTARG;;
        o) outFile=$OPTARG;;
        ?) usage;;
    esac
done
awk -v OFS="\t" '{print $1,$2,$3}' $inFile|bwtool summary -skip-median stdin $bwFile stdout|cut -f8 >${prefix}.tsv
awk -v OFS="\t" -v up=$up -v down=$down '{if($6=="+"){print $1,$2-up,$2}else{print $1,$3,$3+down}}' $inFile |awk -v OFS="\t" '{if($2<0){print $1,"0",$3}else{print $0}}'|bwtool summary -skip-median stdin $bwFile stdout|cut -f8 >${prefix}.up.tsv
awk -v OFS="\t" -v up=$up -v down=$down '{if($6=="+"){print $1,$3,$3+down}else{print $1,$2-up,$2}}' $inFile |awk -v OFS="\t" '{if($2<0){print $1,"0",$3}else{print $0}}'|bwtool summary -skip-median stdin $bwFile stdout|cut -f8 >${prefix}.down.tsv
paste $inFile ${prefix}.tsv ${prefix}.up.tsv ${prefix}.down.tsv |awk -v OFS="\t" '{print $0,$(NF-2)-$(NF-1),$(NF-2)-$NF}' >$outFile

