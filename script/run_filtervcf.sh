#!/usr/bin/bash

usage() {
    echo "
  usage:   run_filtervcf.sh [options]
  options:
    -d    VCF directory
    -p    panel of normals
    -o    output directory
    -h    show this message
  "
}

dir=$(pwd)
panel_of_normals=/Volumes/perso/Analysis/Analysis/Mutect2/panel_of_normals.vcf.gz
output_dir=$(pwd)


while getopts 'd:p:o:' flag; do
  case "${flag}" in
    s)  dir=${OPTARG};;
    f)  panel_of_normals=${OPTARG};;
    o)  output_dir=${OPTARG};;
    h)  usage
    exit 0 ;;
  esac
done


if [ -d $dir ]
then
  echo $dir
  for f in ${dir}/*_merged.vcf
  do
    output_base=$(basename "$f" | cut -d '_' -f 1)
    echo "python script/filtervcf.py -v $f --panel-of-normals $panel_of_normals"
    python script/filtervcf.py -v $f --panel-of-normals $panel_of_normals
  done
fi
