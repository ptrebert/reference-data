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

    args = parser.parse_args()
    return args


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
    tsv_out = sorted(out, key=lambda x: x[1], reverse=True)
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