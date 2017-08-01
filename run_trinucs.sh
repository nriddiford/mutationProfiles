#!/bin/sh

usage() {
    echo "
  usage:   run_trinucs.sh [options]
  options:
    -v    read from varscan native output
    -g    path to genome.fasta
    -h    show this message
  "
}

varscan=0

genome=/Users/Nick_curie/Documents/Curie/Data/Genomes/dmel_6.12.fa

while getopts 'vhg:' flag; do
  case "${flag}" in
    v)  varscan=1 ;;
    g)  genome=${OPTARG};;
    h)  usage
        exit 0 ;;
  esac
done


if [[ -f "data/GW.snv.dist.txt" ]]
then
  echo "Cleaning up old files"
  rm data/GW.snv.dist.txt
  #rm data/GW.trinucs.txt
  #rm data/chroms.trinucs.txt
fi

if [[ $varscan -eq 1 ]]
then
  for varscan_file in data/*.snp
  do
    perl script/trinucs.pl -g $genome -i $varscan_file
  done
fi

if [[ $varscan -eq 0 ]]
then
  for vcf in data/*.vcf
  do
    perl script/trinucs.pl -g $genome -v $vcf
  done
fi

mv *.txt data/
