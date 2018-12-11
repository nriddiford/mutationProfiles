#!/bin/sh

usage() {
    echo "
  usage:   run_trinucs.sh [options]
  options:
    -v    varscan2
    -m    mutect2 file
    -g    path to genome.fasta
    -h    show this message
  "
}

varscan=0
mutect=0

genome=/Users/Nick_curie/Documents/Curie/Data/Genomes/dmel_6.12.fa

while getopts 'vmhg:' flag; do
  case "${flag}" in
    v)  varscan=1 ;;
    m)  mutect=1 ;;
    g)  genome=${OPTARG};;
    h)  usage
        exit 0 ;;
  esac
done


if [[ -f "data/combined_snvs.txt" ]]
then
  echo "Cleaning up old files"
  rm data/combined_snvs.txt
fi

if [[ $varscan -eq 1 ]]
then
  for vcf in data/raw/*varscan*.vcf
  do
    echo "perl script/vcffilter.pl -v $vcf -s varscan -o data"
    perl script/vcffilter.pl -v $vcf -s varscan -o data
  done

  for filt_vcf in data/*varscan_filt.vcf
  do
    echo "perl script/trinucs.pl -g $genome -v $filt_vcf -d data"
    perl script/trinucs.pl -g $genome -v $filt_vcf -d data
  done
fi


if [[ $mutect -eq 1 ]]
then
  for vcf in data/raw/*mutect*.vcf
  do
    echo "perl script/vcffilter.pl -v $vcf -s mutect -o data"
    perl script/vcffilter.pl -v $vcf -s mutect -o data
  done

  for filt_vcf in data/*mutect_filt.vcf
  do
    echo "perl script/trinucs.pl -g $genome -v $filt_vcf -d data"
    perl script/trinucs.pl -g $genome -v $filt_vcf -d data
  done
fi


# mv combined_snvs.txt data
