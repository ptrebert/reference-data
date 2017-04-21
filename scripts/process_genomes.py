#!/usr/bin/env python3
# coding=utf-8

"""
Read 2bit genome files and create various output files
"""

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import re as re

import twobitreader as tbr


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--task', '-tk', type=str, choices=['sizes'], dest='task', required=True)
    parser.add_argument('--2bit', '-2b', type=str, dest='genome', required=True)
    parser.add_argument('--assembly', '-assm', type=str, default='auto', dest='assembly')

    comgroup = parser.add_argument_group('Output chromosome sizes')
    comgroup.add_argument('--tsv-format', '-tsv', type=str, default='', dest='tsv')
    comgroup.add_argument('--bed-format', '-bed', type=str, default='', dest='bed')
    comgroup.add_argument('--select', '-s', type=str, default='"[a-zA-Z0-9\._]+"', dest='select')
    comgroup.add_argument('--sort-order', '-so', type=str, choices=['karyo', 'size'], default='karyo',
                          dest='sortorder', help='Sort by karyotype: natural name sort of chromosomes')

    args = parser.parse_args()
    return args


def check_special_chroms(name):
    """
    :param name:
    :return:
    """
    # check for chimp exceptions
    if name in ['chr2A', '2A', '2a', 'chr2a']:
        return 21
    elif name in ['chr2B', '2B', '2b', 'chr2b']:
        return 22
    elif name in ['chrX', 'X', 'chrZ', 'Z']:
        return 1000
    elif name in ['chrY', 'Y', 'chrW', 'W']:
        return 2000
    elif name in ['chrM', 'M', 'chrMT', 'MT']:
        return 3000
    elif name in ['hs37d5']:
        return 4000
    else:
        return -1


def chrom_karyo_sort(chroms):
    """
    :param chroms:
    :return:
    """
    ordered = []
    unordered = []
    for cname, size in chroms:
        try:
            ord = int(cname.lower().strip('chr'))
            ordered.append((cname, size, ord * 10))
        except ValueError:
            ord = check_special_chroms(cname)
            if ord > 0:
                ordered.append((cname, size, ord))
            else:
                unordered.append((cname, size, -1))
    unordered = sorted(unordered, key=lambda x: x[1], reverse=True)
    ordered = sorted(ordered, key=lambda x: x[2])
    ordered.extend(unordered)
    return [(t[0], t[1]) for t in ordered]


def write_chrom_sizes(data, path, assm):
    """
    :param data:
    :param path:
    :param assm:
    :return:
    """
    if not path:
        return 0
    if path == 'auto':
        wd = os.getcwd()
        fn = '_'.join([assm, 'chrom', 'sizes'])
        if len(data[0]) == 2:
            fn += '.tsv'
        else:
            fn += '.bed'
        path = os.path.join(wd, fn)
    with open(path, 'w') as out:
        _ = out.write('\n'.join(['\t'.join(t) for t in data]) + '\n')
    return 0


def output_chromosome_sizes(args):
    """
    :param args:
    :return:
    """
    if args.assembly == 'auto':
        assm = os.path.basename(args.genome).split('.')[0].strip()
    else:
        assm = args.assembly
    out = []
    selector = re.compile(args.select.strip('"'))
    genome = tbr.TwoBitFile(args.genome)
    for chrom, size in genome.sequence_sizes().items():
        if selector.match(chrom) is not None:
            out.append((chrom, size))
    if args.sortorder == 'size':
        tsv_out = sorted(out, key=lambda x: x[1], reverse=True)
    else:
        tsv_out = chrom_karyo_sort(out)
    tsv_out = [(x[0], str(x[1])) for x in tsv_out]
    bed_out = [(x[0], '0', str(x[1])) for x in tsv_out]
    _ = write_chrom_sizes(tsv_out, args.tsv, assm)
    _ = write_chrom_sizes(bed_out, args.bed, assm)
    return


if __name__ == '__main__':
    try:
        args = parse_command_line()
        choices = {'sizes': output_chromosome_sizes}
        run = choices[args.task]
        run(args)
    except Exception as err:
        trb.print_exc()
        sys.stderr.write('\nError: {}\n'.format(err))
        sys.exit(1)
    else:
        sys.exit(0)