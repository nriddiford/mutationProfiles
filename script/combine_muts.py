from __future__ import division

import fnmatch
import json
import ntpath
import os
import sys
from collections import defaultdict
from optparse import OptionParser

import pandas as pd


def extract_svs(options):
    print("Extracting structural variants from %s" % (options.svs))

    df = pd.read_csv(options.svs, delimiter="\t", index_col=False, na_filter=False)
    print(df.head())
    # for index, row in df.iterrows():
    #     print(row)


    return True

def combine_vars(options):
    svs = extract_svs(options)

    return True


def write_vcf(vars, options):
    tumour, normal = find_normal(options)
    df = pd.DataFrame(vars, columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', tumour, normal])

    chroms = ['2L', '2R', '3L', '3R', '4', 'X', 'Y']

    df = df.sort_values(by=['#CHROM', 'POS'])
    df = df[df['#CHROM'].isin(chroms)]

    output_VCF = '_'.join([tumour, "merged.vcf"])
    with open(output_VCF, 'w') as vcf:
        vcf.write("##fileformat=VCFv4.1\n")
        vcf.write("##source=combinevcf.py\n")
        vcf.write("##somaticSeq=" + options.somseq + "\n")
        vcf.write("##Freebayes=" + options.freebayes + "\n")

    df.to_csv(output_VCF, sep="\t", mode='a', index=False)

    return True


def get_args():
    parser = OptionParser()

    parser.add_option("--structural_variants", dest="svs", action="store", help="Structural variants file")
    parser.add_option("--snvs", dest="snvs", action="store", help="SNV file")
    parser.add_option("--indels", dest="indels", action="store", help="INDEL file")
    parser.add_option("--transposable-elements", dest="tes", action="store", help="TE file")
    parser.add_option("--drivers", dest="drivers", action="store", help="Drivers")
    parser.add_option("--annotate", dest="annotate", action="store", help="Genes to look for mutations in")

    parser.add_option("-o", "--out_file", dest="out_file", action="store", help="File to write annotated vars to")
    parser.set_defaults(svs='/Users/Nick_curie/Desktop/script_test/alleleFreqs/all_samples_snvs_23719.txt',
                        snvs='/Users/Nick_curie/Desktop/script_test/mutationProfiles/data/annotated_snvs.txt',
                        indels='/Users/Nick_curie/Desktop/script_test/mutationProfiles/data/annotated_indels.txt',
                        tes='/Users/Nick_curie/Desktop/Analysis_pipelines/TEs/all_bps.txt', drivers=['N', 'kuz', "Dl"],
                        annotate='/Users/Nick_curie/Desktop/gene2bed/bed/dna_damage_repair.bed')

    options, args = parser.parse_args()

    if not options.annotate:
        parser.print_help()
        print ()
        sys.exit("[!] Must provide both a somaticSeq and Freebayes .vcf files. Exiting.")
    else:
        return options, args


def main():
    options, args = get_args()

    try:
        combine_vars(options)
    except IOError as err:
        sys.stderr.write("IOError " + str(err) + "\n")
        return


if __name__ == "__main__":
    sys.exit(main())
