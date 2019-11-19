from __future__ import division

import fnmatch
import json
import ntpath
import os
import sys
from collections import defaultdict
from optparse import OptionParser

import pandas as pd
import vcf


def find_normal(options):
    sample = ntpath.basename(options.somseq).split("_")[0]
    config = pd.read_csv(options.config, delimiter="\t", index_col=False, na_filter=False, names=['sample', 'assay'])

    samples = config['sample'].tolist()
    it = iter(samples)

    s = {}
    for x in it:
        s[x] = next(it)

    if s[sample]:
        print("Tumour: %s" % sample)
        print("Normal: %s" % s[sample])
        return sample, s[sample]
    else:
        print("Cant find corresponding normal sample for %s" % sample)


def extract_freebayes(options, tumour, normal):
    vcf_reader = vcf.Reader(open(options.freebayes, 'r'))

    freebayes_vars = defaultdict(lambda: defaultdict(dict))
    vars = []
    count = 0

    for record in vcf_reader:
        count += 1
        vaf = round((record.genotype(tumour)['AO'] / (record.genotype(tumour)['AO'] + record.genotype(tumour)['RO'])), 2)
        freebayes_vars[record.CHROM][record.POS] = [record.REF, str(record.ALT[0]), record.INFO, record.genotype(tumour)['GT'], vaf]

    print("%s freebayes vars" % (count))
    return freebayes_vars


def extract_somseq(options, tumour, normal):
    vcf_reader = vcf.Reader(open(options.somseq, 'r'))

    somseq_vars = defaultdict(lambda: defaultdict(dict))

    count = 0

    for record in vcf_reader:
        count += 1
        filter = 'PASS'
        if record.FILTER:
            filter = record.FILTER[0]

        somseq_vars[record.CHROM][record.POS] = [record.REF, str(record.ALT[0]), record.INFO['MVSK'], record.genotype(tumour)['GT'], record.genotype(tumour)['VAF'], filter]

    # print(json.dumps(somseq_vars, indent=4, sort_keys=True))
    print("%s somatic seq vars" % (count))

    return somseq_vars, vcf_reader


def combine_vars(options):

    tumour, normal = find_normal(options)
    freebayes_vars = extract_freebayes(options, tumour, normal)

    somseq_vars, reader = extract_somseq(options, tumour, normal)

    merged_vars = []
    intersect = 0
    for c in somseq_vars:
        for p in somseq_vars[c]:
            ref, alt, callers, gt, vaf, filter = somseq_vars[c][p]
            if freebayes_vars[c][p]:
                intersect += 1
                callers.append(1)
                del freebayes_vars[c][p]
            else:
                callers.append(0)

            numtools = callers.count(1)
            callers = str(callers).strip('[]')

            info = 'MVSKF=%s;TOOLS=%s' % (callers.replace(' ', ''), numtools)
            merged_vars.append([c, p, '.', ref, alt, 0, filter, info, 'GT:VAF', ':'.join([gt, str(vaf)])])

    print("%s vars in both vcfs" % (intersect))
    print("Merged overlaps %s " % (len(merged_vars)))

    additional_vars = 0
    for c in freebayes_vars:
        for p in freebayes_vars[c]:
            if freebayes_vars[c][p]:
                additional_vars += 1
                ref, alt, info, gt, vaf = freebayes_vars[c][p]
                merged_vars.append([c, p, '.', ref, alt, 0, 'PASS', 'MVSKF=0,0,0,0,1;TOOLS=1', 'GT:VAF', ':'.join([gt, str(vaf)])])

    print("Adding %s vars to somseq" % (additional_vars))
    # print(json.dumps(merged_vars, indent=4, sort_keys=True))

    write_vcf(merged_vars, options, tumour)

    return True


def write_vcf(vars, options, tumour):
    dfObj = pd.DataFrame(vars, columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', tumour])
    dfObj = dfObj.sort_values(by=['#CHROM', 'POS'])

    output_VCF = '_'.join([tumour, "merged.vcf"])
    with open(output_VCF, 'w') as vcf:
        vcf.write("##fileformat=VCFv4.1\n")
        vcf.write("##source=combinevcf.py\n")
        vcf.write("##somaticSeq=" + options.somseq + "\n")
        vcf.write("##Freebayes=" + options.freebayes + "\n")

    dfObj.to_csv(output_VCF, sep="\t", mode='a', index=False)


def get_args():
    parser = OptionParser()

    parser.add_option("--somseq", dest="somseq", action="store", help="somaticSeq vcf file")
    parser.add_option("--freebayes", dest="freebayes", action="store", help="Freebayes vcf file")
    parser.add_option("--config", dest="config", action="store", help="mapping for tumour/normal samples")
    parser.add_option("-o", "--out_file", dest="out_file", action="store", help="File to write annotated vars to")
    parser.set_defaults(config='/Users/Nick_curie/Desktop/script_test/alleleFreqs/data/samples.tsv')

    options, args = parser.parse_args()

    if not options.somseq or not options.freebayes:
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
