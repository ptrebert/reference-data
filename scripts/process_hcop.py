#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import gzip as gz
import argparse as argp
import fnmatch as fnm

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
    parser.add_argument('--annotation', '-ann', type=str, dest='annotation')
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
    modelfiles = os.listdir(basedir)
    modelfile = fnm.filter(modelfiles, '{}_*.genes.bed.gz'.format(abbr))
    assert len(modelfile) == 1, 'No gene set identified: {} / {}'.format(species, abbr)
    fpath = os.path.join(basedir, modelfile[0])
    genes = pd.read_csv(fpath, skip_blank_lines=True, header=0,
                        na_values='n/a', sep='\t', dtype=str,
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
    df.dropna(axis=0, how='any', inplace=True,
              subset=['human_ensembl_gene', '{}_ensembl_gene'.format(other)])
    new_cols = []
    for c in df.columns:
        if c.endswith('_ensembl_gene'):
            spec, _, _ = c.split('_')
            new_cols.append('{}_name'.format(spec))
        else:
            new_cols.append(c)
    df.columns = new_cols
    return df, other


def dump_hcop_data(dataset, other, group_root, outpath, mode):
    """
    :param dataset:
    :param other:
    :param other:
    :param mode:
    :return:
    """
    if other:
        group = os.path.join(group_root, 'human', other)
    else:
        group = group_root
    with pd.HDFStore(outpath, mode, complib='blosc', complevel=9) as hdf:
        hdf.put(group, dataset, format='table')
        hdf.flush()
    return


def select_ortholog_pairs(hcop, human, other, species, locfilter):
    """
    :param hcop:
    :param human:
    :param other:
    :param locfilter:
    :return:
    """
    other_name = '{}_name'.format(species)

    human_select = human['chrom'].str.match(locfilter, as_indexer=True)
    other_select = other['chrom'].str.match(locfilter, as_indexer=True)
    human_subset = human.loc[human_select, 'human_name'].unique()
    other_subset = other.loc[other_select, other_name].unique()

    human_idx = hcop['human_name'].isin(human_subset)
    other_idx = hcop[other_name].isin(other_subset)
    joined_idx = np.logical_and(human_idx, other_idx)

    hcop_subset = hcop.loc[joined_idx, ['human_name', other_name]]
    uniq_matches = hcop.groupby('human_name').count()
    uniq_human_names = uniq_matches.loc[uniq_matches[other_name] == 1, :].index
    # sanity check
    raw_size = uniq_human_names.size
    uniq_size = uniq_human_names.unique().size
    assert raw_size == uniq_size, 'Duplicate gene names still in list: {} vs {}'.format(raw_size, uniq_size)

    uniq_pairs = hcop_subset.loc[hcop_subset['human_name'].isin(uniq_human_names), :]
    uniq_pairs = uniq_pairs.merge(human, how='left', on='human_name')
    uniq_pairs.drop('chrom', axis=1, inplace=True)

    uniq_pairs = uniq_pairs.merge(other, how='left', on=other_name)
    uniq_pairs.drop('chrom', axis=1, inplace=True)
    uniq_pairs.reset_index(drop=True, inplace=True)
    return uniq_pairs


def select_ortholog_groups(fpath, group_root):
    """
    :param fpath:
    :param group_root:
    :return:
    """
    if not group_root.startswith('/'):
        group_root = '/' + group_root
    groups = None
    with pd.HDFStore(fpath, 'r') as hdf:
        load_keys = [k for k in hdf.keys() if k.startswith(group_root)]
        for k in load_keys:
            if groups is None:
                groups = hdf[k]
            else:
                groups = groups.merge(hdf[k], on=['human_name', 'human_symbol'],
                                      how='outer', suffixes=('', ''))
    groups.dropna(axis=0, how='any', inplace=True)
    groups.reset_index(drop=True, inplace=True)
    return groups


def main():
    """
    :return:
    """
    args = parse_command_line()
    outmode = 'w'
    human_genes = read_gene_model(args.annotation, 'human')
    for fp in args.inputfiles:
        data, other = read_hcop_table(fp)
        dump_hcop_data(data, other, 'raw', args.output, outmode)
        outmode = outmode.replace('w', 'a')
        other_genes = read_gene_model(args.annotation, other)
        pairs = select_ortholog_pairs(data.copy(), human_genes.copy(), other_genes.copy(),
                                      other, 'chr[0-9XYZW]+$')
        dump_hcop_data(pairs, other, 'augo/pairs', args.output, outmode)
        pairs = select_ortholog_pairs(data.copy(), human_genes.copy(), other_genes.copy(),
                                      other, 'chr[0-9]+$')
        dump_hcop_data(pairs, other, 'auto/pairs', args.output, outmode)

    groups = select_ortholog_groups(args.output, 'auto/pairs')
    dump_hcop_data(groups, '', 'auto/groups', args.output, outmode)

    groups = select_ortholog_groups(args.output, 'augo/pairs')
    dump_hcop_data(groups, '', 'augo/groups', args.output, outmode)
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