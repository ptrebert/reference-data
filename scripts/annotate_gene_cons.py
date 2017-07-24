#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import collections as col

import pandas as pd
import numpy as np


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--gene-body', '-gb', type=str, dest='genebody')
    parser.add_argument('--cons-scores', '-cs', type=str, dest='consscores')
    parser.add_argument('--window-size', '-ws', type=int, nargs='+', default=[50, 100, 200], dest='windowsize')
    parser.add_argument('--output', '-o', type=str, dest='output')

    args = parser.parse_args()
    return args


def load_genes(fpath):
    """
    :param fpath:
    :return:
    """
    genes = []
    with pd.HDFStore(fpath, 'r') as hdf:
        for k in hdf.keys():
            if k == '/metadata':
                continue
            chrom = os.path.split(k)[-1]
            data = hdf[k]
            data['chrom'] = chrom
            genes.append(data)
    genes = pd.concat(genes, axis=0, ignore_index=True)
    genes = genes.replace(['-', '+'], [-1, 1], inplace=False)
    genes['tss'] = 0
    genes.loc[genes.strand > 0, 'tss'] = genes.loc[genes.strand > 0, 'start']
    genes.loc[genes.strand < 0, 'tss'] = genes.loc[genes.strand < 0, 'end']
    genes.drop(['score'], axis=1, inplace=True)
    return genes


def load_cons_scores(fpath, chrom):
    """
    :param fpath:
    :param chrom:
    :return:
    """
    with pd.HDFStore(fpath, 'r') as hdf:
        load_key = [k for k in hdf.keys() if k.endswith(chrom)]
        load_key = load_key[0]
        data = hdf[load_key]
    return data


def annotate_conservation(genes, scorefile, windows):
    """
    :param genes:
    :param scorefile:
    :param windows:
    :return:
    """
    for chrom in genes['chrom'].unique():
        scores = load_cons_scores(scorefile, chrom)
        collector = col.defaultdict(list)
        for gene in genes.loc[genes['chrom'] == chrom, :].itertuples():
            for w in windows:
                avg = scores[gene.tss - w:gene.tss + w].mean()
                collector['ext{}'.format(w)].append(avg)
        for w in windows:
            genes.loc[genes['chrom'] == chrom, 'ext{}'.format(w)] = collector['ext{}'.format(w)]
    for w in windows:
        genes['rnk{}'.format(w)] = genes['ext{}'.format(w)].rank(method='min', pct=True, ascending=True)
    for w in windows:
        genes['grp{}'.format(w)] = 0
    steps = list(range(5, 100, 5))
    steps.extend([99, 99.9])
    steps = np.array(steps, dtype=np.float16)
    steps /= 100.
    for s in steps:
        for w in windows:
            genes.loc[genes['rnk{}'.format(w)] >= s, 'grp{}'.format(w)] += 1
    return genes


def main():
    """
    :return:
    """
    args = parse_command_line()
    genes = load_genes(args.genebody)
    for w in args.windowsize:
        genes['ext{}'.format(w)] = 0
    genes = annotate_conservation(genes, args.consscores, args.windowsize)
    with pd.HDFStore(args.output, 'w', complib='blosc', complevel=9) as hdf:
        hdf.put('genes', genes, format='table')
        hdf.flush()
    return


if __name__ == '__main__':
    try:
        main()
    except Exception as err:
        trb.print_exc(file=sys.stderr)
        sys.stderr.write('\nError: {}\n'.format(str(err)))
        sys.exit(1)
    else:
        sys.exit(0)
