#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import gzip as gz
import argparse as argp

import pandas as pd
import numpy as np


SPECIES_MAP = {'human': 'hsa',
               'mouse': 'mmu',
               'cow': 'bta',
               'chicken': 'gga',
               'dog': 'cfa',
               'pig': 'ssc'}


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, nargs='+', dest='inputfiles')
    parser.add_argument('--annotations', '-ann', type=str, dest='annotations')
    parser.add_argument('--output', '-o', type=str, dest='output')
    args = parser.parse_args()
    return args


def read_gene_model(basedir, species):
    """
    :param basedir:
    :param species:
    :return:
    """
    abbr = SPECIES_MAP[species]
    modelfile = os.listdir(basedir, '{}_*.genes.bed.gz'.format(abbr))
    assert len(modelfile) == 1, 'No gene set identified: {} / {}'.format(species, abbr)
    fpath = os.path.join(basedir, modelfile[0])
    genes = pd.read_csv(fpath, skip_blank_lines=True, header=0,
                        na_values='n/a', sep='\t',
                        usecols=['#chrom', 'name', 'symbol'])
    new_cols = []
    for c in genes.columns:
        if c == '#chrom':
            new_cols.append('chrom')
        else:
            new_cols.append('{}_{}'.format(species, c))
    genes.columns = new_cols
    return genes


def read_hcop_table(fpath):
    """
    :param fpath:
    :return:
    """
    df = pd.read_csv(fpath, delimiter='\t', header=0, skip_blank_lines=True,
                     na_values='-', dtype=str)
    species = set([c.split('_')[0] for c in df.columns])
    assert len(species) == 2, 'Could not identify species: {}'.format(df.columns.tolist())
    assert 'human' in species, 'Human is not part of species set: {}'.format(species)
    other = [s for s in species if s != 'human'][0]
    df.dropna(axis=1, how='any', inplace=True,
              subset=['human_ensembl_gene', '{}_ensembl_gene'.format(other)])
    new_cols = []
    for c in df.columns:
        if c.endswith('_gene_name'):
            spec, _, _ = c.split('_')
            new_cols.append('{}_name'.format(spec))
        else:
            new_cols.append(c)
    df.columns = new_cols
    return df, other


def dump_raw_hcop(dataset, other, outpath, mode):
    """
    :param dataset:
    :param other:
    :param other:
    :param mode:
    :return:
    """
    group = os.path.join('raw', 'human', other)
    with pd.HDFStore(outpath, mode, complib='blosc', complevel=9) as hdf:
        hdf.put(group, dataset, format='fixed')
        hdf.flush()
    return


def main():
    """
    :return:
    """
    args = parse_command_line()
    outmode = 'w'
    for fp in args.inputfiles:
        data, other = read_hcop_table(fp)
        dump_raw_hcop(data, other, args.output, outmode)
        outmode = outmode.replace('w', 'a')
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