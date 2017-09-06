#!/bin/sh

usage() {
    echo "
usage:   run_trinucs.sh [options]
options:
  -v    read from varscan native output
  -h    show this message
"
}

varscan=0

while getopts 'vh' flag; do
  case "${flag}" in
    v)  varscan=1 ;;
    h)  usage
        exit 0 ;;
  esac
done

rm *.txt

if [[ $varscan -eq 1 ]]
then
  for varscan_file in data/*.snp
  do
    perl trinucs.pl -g /Users/Nick_curie/Documents/Curie/Data/Genomes/dmel_6.12.fa -i $varscan_file
  done
fi

if [[ $varscan -eq 0 ]]
then
  for vcf in data/*.vcf
  do
    perl trinucs.pl -g /Users/Nick_curie/Documents/Curie/Data/Genomes/dmel_6.12.fa -v $vcf
  done
fi
