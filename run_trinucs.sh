#!/bin/sh

usage() {
    echo "
  usage:   run_trinucs.sh [options]
  options:
    -v    varscan2
    -m    mutect2 file
    -a    annotate variants
    -g    path to genome.fasta
    -f    path to feature.gtf
    -h    show this message
  "
}

varscan=0
mutect=0
indel=0
annotate=0

genome=/Users/Nick_curie/Documents/Curie/Data/Genomes/dmel_6.12.fa
#genome=/Users/Nick/Documents/Curie/Data/Genomes/Dmel_v6.12/Dmel_6.12.fasta # home

features=/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf
#features=/Users/Nick/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf # home

while getopts 'vmaihg:' flag; do
  case "${flag}" in
    v)  varscan=1 ;;
    m)  mutect=1 ;;
    i)  indel=1 ;;
    g)  genome=${OPTARG};;
    a)  annotate=1 ;;
    h)  usage
        exit 0 ;;
  esac
done

if [[ $# -eq 0 ]]
then
  usage
  exit 0
fi

# if [[ -f "data/combined_snvs.txt" ]]
# then
#   echo "Cleaning up old files"
#   rm data/combined_snvs.txt
# fi

if [[ $varscan -eq 1 ]]
then
  for vcf in data/raw/snpEff/*_varscan_ann.vcf
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

if [[ $indel -eq 1 ]]
then
  for vcf in data/raw/indel/snpEff/*.vcf
  do
    echo "perl script/vcffilter.pl -v $vcf -s indel -o data"
    perl script/vcffilter.pl -v $vcf -s indel -o data
  done

  for filt_vcf in data/*indel_filt.vcf
  do
    echo "perl script/trinucs.pl -g $genome -v $filt_vcf -t indel -d data -o combined_indels.txt"
    perl script/trinucs.pl -g $genome -v $filt_vcf -t indel -d data -o combined_indels.txt
  done
fi

if [[ $mutect -eq 1 ]]
then
  for vcf in data/raw/snpEff/*_mutect_ann.vcf
  do
    # bname=echo ${vcf##*/}
    # name=${bname%_mutect_ann.vcf}
    name=$(basename "$vcf" | cut -d '_' -f1)
    echo "bcftools norm -Ov -m-any $vcf > data/${name}_mutect_norm.vcf"
    bcftools norm -Ov -m-any $vcf > data/${name}_mutect_norm.vcf
    echo "perl script/vcffilter.pl -v data/${name}_mutect_norm.vcf -s mutect -o data"
    perl script/vcffilter.pl -v data/${name}_mutect_norm.vcf -s mutect -o data
  done

  for filt_vcf in data/*mutect_filt.vcf
  do
    echo "perl script/trinucs.pl -g $genome -v $filt_vcf -c mutect -d data"
    perl script/trinucs.pl -g $genome -v $filt_vcf -c mutect -d data
  done
fi

if [[ $annotate -eq 1 ]]
then
  echo "Annotating vars"
  echo "perl script/snv2gene.pl -i data/combined_snvs.txt -f $features"
  perl script/snv2gene.pl -i data/combined_snvs.txt -f $features
fi

# mv combined_snvs.txt data
