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
    sample = ntpath.basename(options.vcf_in).split("_")[0]
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


def filter_all(pon, options):
    dir = os.path.abspath(options.dir)
    print("Looking for .vcf files in directory: %s" % (dir))

    for file in os.listdir(dir):
        if file.endswith("_merged.vcf"):
            options.vcf_in = os.path.join(dir, file)
            print(options.vcf_in)
            filter_vars(pon, options)

    return True


def filter_vars(pon, options):
    tumour, normal = find_normal(options)

    mode = 'r'
    if options.vcf_in.endswith('.gz'):
        mode = 'rb'
    vcf_reader = vcf.Reader(open(options.vcf_in, mode))

    count = 0
    pon_filtered = 0
    depth_filtered = 0

    somatic_vars = []

    for record in vcf_reader:
        count += 1
        key = '_'.join(map(str, [record.CHROM, record.POS, str(record.ALT[0])]))

        if key in pon:
            pon_filtered += 1
            continue

        if record.genotype(normal)['DP'] < 20:
            depth_filtered += 1
            continue

        filter = 'PASS'
        if record.FILTER:
            filter = record.FILTER[0]

        tumour_gt, normal_gt = (record.genotype(tumour)['GT'], record.genotype(normal)['GT'])
        tumour_vaf, normal_vaf = (record.genotype(tumour)['VAF'], record.genotype(normal)['VAF'])
        tumour_dp, normal_dp = (record.genotype(tumour)['DP'], record.genotype(normal)['DP'])

        t_dict = {'GT': tumour_gt, 'VAF': tumour_vaf, 'DP': tumour_dp}
        n_dict = {'GT': normal_gt, 'VAF': normal_vaf, 'DP': normal_dp}

        t_parts = ':'.join(map(str, [t_dict['GT'], t_dict['VAF'], t_dict['DP'] ]))
        n_parts = ':'.join(map(str, [n_dict['GT'], n_dict['VAF'], n_dict['DP'] ]))

        callers = str(record.INFO['MVSKF']).strip('[]')
        callers = callers.replace(' ', '').replace("'", "")
        ntools = record.INFO['TOOLS'][0]

        info = 'MVSKF=%s;TOOLS=%s' % (callers, ntools)

        somatic_vars.append([record.CHROM, record.POS, '.', record.REF, str(record.ALT[0]), 0, filter, info, 'GT:VAF:DP', t_parts, n_parts])

    print("o %s/%s PON vars filtered" % (pon_filtered, count))
    count = count - pon_filtered
    print("o %s/%s low depth (< 20 in normal) vars filtered" % (depth_filtered, count))
    count = count - depth_filtered
    print("o %s remaining vars in '%s'" % (count, options.vcf_in))

    # print(json.dumps(somatic_vars, indent=4, sort_keys=True))

    write_vcf(somatic_vars, options, tumour, normal)

    return True


def write_vcf(vars, options, tumour, normal):
    header = get_header(options)

    df = pd.DataFrame(vars, columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', tumour, normal])

    chroms = ['2L', '2R', '3L', '3R', '4', 'X', 'Y']

    df = df.sort_values(by=['#CHROM', 'POS'])
    df = df[df['#CHROM'].isin(chroms)]

    vcf_out_file = '_'.join([tumour, "consensus_filt.vcf"])

    with open(vcf_out_file, 'w') as vcf:
        vcf.write("\n".join(header))

    df.to_csv(vcf_out_file, sep="\t", mode='a', index=False)

    return True


def get_header(options):
    header = []

    with open(options.vcf_in, 'r') as vcf:
        for line in vcf:
            line = line.strip()
            if line.startswith("##"):
                header.append(line)
    header.append("##PON_filtered=TRUE\n")

    return header


def extract_germline(options):
    print("Reading panel of normals file '%s'" % (options.pon))
    mode = 'r'
    if options.pon.endswith('.gz'):
        mode = 'rb'
    vcf_reader = vcf.Reader(open(options.pon, mode))

    germline_vars = defaultdict(int)

    count = 0

    for record in vcf_reader:
        if record.CHROM not in options.chroms: continue
        count += 1
        if count % 1e5 == 0:
            print("...Read %s variants from '%s'" % (count, options.pon))
        key = '_'.join(map(str, [record.CHROM, record.POS, record.ALT[0]]))
        germline_vars[key] += 1

    # print(json.dumps(germline_vars, indent=4, sort_keys=True))
    print("o Finished reading %s germline variants from panel of normals" % (count))
    return germline_vars


def get_args():
    parser = OptionParser()

    parser.add_option("-v", "--vcf", dest="vcf_in", action="store", help="vcf file")
    parser.add_option("-d", "--dir", dest="dir", action="store", help="Directory containing vcf files to filter")
    parser.add_option("--panel-of-normals", dest="pon", action="store", help="PON file [Produced by Mutect2")
    parser.add_option("--chroms", dest="chroms", action="store_true", help="Filter on chroms")
    parser.add_option("--config", dest="config", action="store", help="mapping for tumour/normal samples")
    parser.add_option("-o", "--out_file", dest="out_file", action="store", help="File to write annotated vars to")
    parser.set_defaults(pon='/Volumes/perso/Analysis/Analysis/Mutect2/panel_of_normals.vcf.gz',
                        config='/Users/Nick_curie/Desktop/script_test/alleleFreqs/data/samples.txt',
                        chroms=['2L', '2R', '3L', '3R', '4', 'X', 'Y'])

    options, args = parser.parse_args()

    if not options.vcf_in and not options.dir:
        parser.print_help()
        print()
        sys.exit("[!] Must provide a vcf file [-v, --vcf] or directory of vcf files [-d, --dir] as input. Exiting.")
    else:
        return options, args


def main():
    options, args = get_args()

    try:
        pon = extract_germline(options)
        if options.dir:
            filter_all(pon, options)
        else:
            filter_vars(pon, options)
    except IOError as err:
        sys.stderr.write("IOError " + str(err) + "\n")
        return


if __name__ == "__main__":
    sys.exit(main())
