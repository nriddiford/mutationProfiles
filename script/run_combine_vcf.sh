#!/usr/bin/bash

usage() {
    echo "
  usage:   run_trinucs.sh [options]
  options:
    -f    Freebayes directory
    -s    SomaticSeq directory
    -o    output directory
    -h    show this message
  "
}

freebayes=/Volumes/perso/Analysis/Analysis/Freebayes/vcf
somaticSeq=/Volumes/perso/Analysis/Analysis/somaticSeq/vcf
output_dir=$(pwd)

while getopts 'f:s:o:' flag; do
  case "${flag}" in
    s)  somaticSeq=${OPTARG};;
    f)  freebayes=${OPTARG};;
    o)  output_dir=${OPTARG};;
    h)  usage
    exit 0 ;;
  esac
done

# if [[ $# -eq 0 ]]
# then
#   usage
#   exit 0
# fi

if [ -d $somaticSeq ]
then
  echo $somaticSeq
  for f in ${somaticSeq}/*_somatic_snps.vcf
  do
    output_base=$(basename "$f" | cut -d '_' -f1)
    echo $output_base
    echo $f
    echo "$freebayes/${output_base}_freebayes.vcf.gz"
    python script/combinevcf.py --somseq $f --freebayes $freebayes/${output_base}_freebayes.vcf.gz
  done
fi
