#!/bin/sh

usage(){
  cat << EOF
  Description: This script will calculate the profile surrounding TSS using bwtool aggregate.
  Usage: sh $0 -u 1500 -d 1500 -g <*.gpe> -b <*.bw> -o <TSS_profile.tsv>
  Note: Make sure that bwtool is in your current working PATH.
  Options:
    -u  INT     The upstream length of TSS to analysize.[default:1500]
    -d  INT     The downstream length of TSS to analysize.[default:1500]
    -g  FILE    The gene structure file in GPE format.
    -b  FILE    The value stored in bigWig file.
    -o  FILE    The output file name.
    -h          Show this help information
EOF
    exit 0
}

[ $1 ] || usage

while getopts "hu:d:g:b:o:" OPTION
do
    case $OPTION in
        h) usage;;
        u) upstream=$OPTARG;;
        d) downstream=$OPTARG;;
        g) gpeFile=$OPTARG;;
        b) bwFile=$OPTARG;;
        o) output=$OPTARG;;
        ?) usage;;
    esac
done

if [ $(head -n 1 $gpeFile | awk '{print NF}') -eq 15 ];then
    awk -v OFS="\t" '{if($3=="+"){print $2,$4,$4+1,$12,"0",$3}else{print $2,$5-1,$5,$12,"0",$3}}' $gpeFile >TSS.bed
else
    if [ $(head -n 1 $gpeFile | awk '{print NF}') -eq 16 ];then  
      cut -f2- $gpeFile | awk -v OFS="\t" '{if($3=="+"){print $2,$4,$4+1,$12,"0",$3}else{print $2,$5-1,$5,$12,"0",$3}}' $gpeFile >TSS.bed
    else
      echo "Please offer the gene structure file in correct gpe format" >&2
      exit 1
    fi
fi
bwtool agg ${downstream}:$upstream TSS.bed $bwFile $output
rm TSS.bed

