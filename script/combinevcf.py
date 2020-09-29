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
    mode = 'r'
    if options.freebayes.endswith('.gz'):
        mode = 'rb'

    vcf_reader = vcf.Reader(open(options.freebayes, mode))

    freebayes_vars = defaultdict(lambda: defaultdict(dict))
    vars = []
    count = 0

    for record in vcf_reader:
        count += 1
        tumour_vaf = round((record.genotype(tumour)['AO'] / (record.genotype(tumour)['AO'] + record.genotype(tumour)['RO'])), 2)
        normal_vaf = round((record.genotype(normal)['AO'] / (record.genotype(normal)['AO'] + record.genotype(normal)['RO'])), 2)
        tumour_gt, normal_gt = (record.genotype(tumour)['GT'], record.genotype(normal)['GT'])
        tumour_dp, normal_dp = (record.genotype(tumour)['DP'], record.genotype(normal)['DP'])

        t_dict = {'GT': tumour_gt, 'VAF': tumour_vaf, 'DP': tumour_dp}
        n_dict = {'GT': normal_gt, 'VAF': normal_vaf, 'DP': normal_dp}

        freebayes_vars[record.CHROM][record.POS] = [record.REF, str(record.ALT[0]), record.INFO, t_dict, n_dict]

    print("%s freebayes vars" % (count))
    return freebayes_vars


def extract_somseq(options, tumour, normal):
    mode = 'r'
    if options.somseq.endswith('.gz'):
        mode = 'rb'

    vcf_reader = vcf.Reader(open(options.somseq, mode))

    somseq_vars = defaultdict(lambda: defaultdict(dict))

    count = 0

    for record in vcf_reader:
        count += 1
        filter = 'PASS'
        if record.FILTER:
            filter = record.FILTER[0]

        tumour_gt, normal_gt = (record.genotype(tumour)['GT'], record.genotype(normal)['GT'])
        tumour_vaf, normal_vaf = (record.genotype(tumour)['VAF'], record.genotype(normal)['VAF'])
        tumour_dp, normal_dp = (record.genotype(tumour)['DP4'], record.genotype(normal)['DP4'])



        t_dict = {'GT': tumour_gt, 'VAF': tumour_vaf, 'DP': sum(tumour_dp)}
        n_dict = {'GT': normal_gt, 'VAF': normal_vaf, 'DP': sum(normal_dp)}

        somseq_vars[record.CHROM][record.POS] = [record.REF, str(record.ALT[0]), record.INFO['MVSK'], t_dict, n_dict, filter]

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
            ref, alt, callers, tumour, normal, filter = somseq_vars[c][p]
            if freebayes_vars[c][p]:
                intersect += 1
                callers.append(1)
                del freebayes_vars[c][p]
            else:
                callers.append(0)

            numtools = callers.count(1)
            callers = str(callers).strip('[]')

            info = 'MVSKF=%s;TOOLS=%s' % (callers.replace(' ', ''), numtools)
            t_parts = ':'.join(map(str, [tumour['GT'], tumour['VAF'], tumour['DP'] ]))
            n_parts = ':'.join(map(str, [normal['GT'], normal['VAF'], normal['DP'] ]))

            merged_vars.append([c, p, '.', ref, alt, 0, filter, info, 'GT:VAF:DP', t_parts, n_parts])

    print("%s vars in both vcfs" % (intersect))
    print("Merged overlaps %s " % (len(merged_vars)))

    additional_vars = 0
    for c in freebayes_vars:
        for p in freebayes_vars[c]:
            if freebayes_vars[c][p]:
                additional_vars += 1
                ref, alt, info, t_dict, n_dict = freebayes_vars[c][p]
                t_parts = ':'.join(map(str, [tumour['GT'], tumour['VAF'], tumour['DP']]))
                n_parts = ':'.join(map(str, [normal['GT'], normal['VAF'], normal['DP'] ]))

                merged_vars.append([c, p, '.', ref, alt, 0, 'PASS', 'MVSKF=0,0,0,0,1;TOOLS=1', 'GT:VAF:DP', t_parts, n_parts])

    print("Adding %s vars to somseq" % (additional_vars))

    # print(json.dumps(merged_vars, indent=4, sort_keys=True))

    write_vcf(merged_vars, options)

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
