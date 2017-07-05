#!/usr/bin/env python
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import re as re
import argparse as argp
import itertools as itt
import multiprocessing as mp
import collections as col
import pickle as pck

import twobitreader as tbr
import pandas as pd
import numpy as np
import numpy.random as rng


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--enhancer', '-e', type=str, dest='enhancer')
    parser.add_argument('--genes', '-g', type=str, dest='genes')
    parser.add_argument('--sequence', '-s', type=str, dest='sequence')
    parser.add_argument('--chrom', '-c', type=str, dest='chromosomes', default='"chr[0-9XY]+$"')
    parser.add_argument('--workers', '-w', type=int, dest='workers', default=4)
    mgrp = parser.add_mutually_exclusive_group()
    mgrp.add_argument('--train', action='store_true', default=False, dest='train')
    mgrp.add_argument('--test', action='store_true', default=False, dest='test')
    parser.add_argument('--output', '-o', type=str, dest='output')
    parser.add_argument('--cv-index', '-i', type=str, dest='cvindex')
    args = parser.parse_args()
    return args


def load_genes(fpath, chrom, subset):
    """
    :param fpath:
    :param chrom:
    :return:
    """
    data = None
    with pd.HDFStore(fpath, 'r') as hdf:
        for k in hdf.keys():
            if k.endswith(chrom):
                data = hdf[k]
                data['tss'] = data['start'] + ((data['end'] - data['start']) / 2).round(decimals=0)
                data = data.loc[:, ['tss', 'start', 'end', 'name', 'symbol']]
                if subset:
                    name_idx = data['name'].isin(subset)
                    symbol_idx = data['symbol'].isin(subset)
                    joined_idx = np.logical_or(name_idx, symbol_idx)
                    data = data.loc[joined_idx, :]
                    if data.empty:
                        data = None

                break
    if data is not None:
        data.sort_values(by='tss', axis=0, ascending=True, inplace=True)
        assert 'tss' in data.columns, 'Missing required tss column: {}'.format(data.columns)
    return data


def load_enhancer(fpath, chrom):
    """
    :param fpath:
    :param chrom:
    :return:
    """
    assert os.path.isfile(fpath), 'Invalid path to enhancer file: {}'.format(fpath)
    with open(fpath, 'r') as bedfile:
        header = bedfile.readline().split()
        header[0] = 'chrom'
        bedfile.seek(0)
        enh = pd.read_csv(bedfile, header=0, delimiter='\t', names=header)
    enh = enh.loc[enh['chrom'] == chrom, :]
    if enh.empty:
        enh = None
        assoc_genes = {}
    elif 'name' in enh.columns or 'symbol' in enh.columns:
        enh.sort_values(by=['start', 'end'], axis=0, ascending=True, inplace=True)
        enh['mid'] = enh['start'] + ((enh['end'] - enh['start']) / 2).round(decimals=0)
        assoc_genes = set(itt.chain.from_iterable(enh['name'].str.split(',').tolist()))
        assoc_genes = assoc_genes.union(set(itt.chain.from_iterable(enh['symbol'].str.split(',').tolist())))
        assoc_genes = assoc_genes - {'n/a'}
    else:
        enh['mid'] = enh['start'] + ((enh['end'] - enh['start']) / 2).round(decimals=0)
        assoc_genes = {}
    if enh is not None:
        assert 'mid' in enh.columns, 'Missing required mid column: {}'.format(enh.columns)
    return enh, assoc_genes


def load_sequence(fpath, chrom, re_chrom=''):
    """
    :param fpath:
    :param chrom:
    :param re_chrom:
    :return:
    """
    genome = tbr.TwoBitFile(fpath)
    if chrom == 'all':
        chrom_match = re.compile(re_chrom.strip('"'))
        ret = [c for c in genome.keys() if chrom_match.match(c) is not None]
    else:
        obj = genome[chrom]
        ret = str(obj).lower()
    return ret


def comp_seq_features(seq, suffix):
    """
    :param seq:
    :param suffix:
    :return:
    """
    slen = float(len(seq))
    total_gc = (seq.count('c') + seq.count('g')) / slen
    total_cpg = seq.count('cg') / (slen / 2)
    # https://doi.org/10.1101/111955
    total_tpa = seq.count('ta') / (slen / 2)
    total_cpa = seq.count('ca') / (slen / 2)
    total_gpa = seq.count('ga') / (slen / 2)
    rec = {'ftgc_pct_GC_' + suffix: total_gc * 100.,
           'ftdim_pct_CpG_' + suffix: total_cpg * 100.,
           'ftdim_pct_TpA_' + suffix: total_tpa * 100.,
           'ftdim_pct_CpA_' + suffix: total_cpa * 100.,
           'ftdim_pct_GpA_' + suffix: total_gpa * 100.}
    return rec


def prepare_training_data(params):
    """
    :param params:
    :return:
    """
    chrom, cargs, lock = params
    assert cargs.train, 'Wrong mode'
    enh_set, assoc_genes = load_enhancer(cargs.enhancer, chrom)
    gene_set = load_genes(cargs.genes, chrom, assoc_genes)
    if enh_set is None or gene_set is None:
        return chrom, None
    seq = load_sequence(cargs.sequence, chrom)
    records = []
    for gene in gene_set.itertuples():
        base_feat = {'name': gene.name, 'symbol': gene.symbol}
        base_feat.update(comp_seq_features(seq[gene.start:gene.end], 'reg5p'))
        gene_nbh = gene_set.query('tss >= @gene.tss - 1000000 & tss <= @gene.tss + 1000000')
        base_feat['ftnbh_abs_tss'] = gene_nbh.shape[0] - 1
        enh_nbh = enh_set.query('mid >= @gene.tss - 1000000 & mid <= @gene.tss + 1000000')
        if enh_nbh.empty:
            continue
        base_feat['ftnbh_abs_enh'] = enh_nbh.shape[0]
        for enhancer in enh_nbh.itertuples():
            pair_feat = dict(base_feat)
            enh_names = enhancer.name.split(',')
            enh_symbols = enhancer.symbol.split(',')
            if gene.name in enh_names:
                idx = enh_names.index(gene.name)
                pair_feat['output'] = float(enhancer.assoc_score.split(',')[idx])
            elif gene.symbol in enh_symbols:
                idx = enh_symbols.index(gene.symbol)
                pair_feat['output'] = float(enhancer.assoc_score.split(',')[idx])
            else:
                pair_feat['output'] = 0.
            pair_feat.update(comp_seq_features(seq[enhancer.start:enhancer.end], 'enh'))
            pair_feat['ftenh_abs_score'] = enhancer.enhancer_score
            pair_feat['ftenh_abs_elite'] = enhancer.is_elite
            pair_feat['enh_name'] = enhancer.GHid
            pair_feat['enh_start'] = enhancer.start
            pair_feat['enh_end'] = enhancer.end
            pair_feat['ftenh_abs_dist'] = enhancer.mid - gene.tss
            pair_feat['ftcons_abs_enh_mean'] = enhancer.ftcons_enh_abs_mean
            pair_feat['ftcons_abs_enh_median'] = enhancer.ftcons_enh_abs_median
            pair_feat['ftcons_abs_enh_max'] = enhancer.ftcons_enh_abs_max
            pair_feat['ftcons_abs_enh_min'] = enhancer.ftcons_enh_abs_min
            records.append(pair_feat)
    dataset = pd.DataFrame.from_records(records, index=range(len(records)))
    # shuffle once
    rand_index = dataset.index.values
    rng.shuffle(rand_index)
    dataset = dataset.loc[rand_index, :]
    dataset.reset_index(drop=True, inplace=True)
    with lock:
        with pd.HDFStore(cargs.output, 'a', complib='blosc', complevel=9) as hdf:
            hdf.put(chrom, dataset, format='table')
            hdf.flush()
    return chrom, dataset.shape[0]


def prepare_testing_data(params):
    """
    :param params:
    :return:
    """
    chrom, cargs, lock = params
    assert cargs.test, 'Wrong mode'
    enh_set, assoc_genes = load_enhancer(cargs.enhancer, chrom)
    gene_set = load_genes(cargs.genes, chrom, set())
    if enh_set is None or gene_set is None:
        return chrom, None
    seq = load_sequence(cargs.sequence, chrom)
    records = []
    for gene in gene_set.itertuples():
        base_feat = {'name': gene.name, 'symbol': gene.symbol}
        base_feat.update(comp_seq_features(seq[gene.start:gene.end], 'reg5p'))
        gene_nbh = gene_set.query('tss >= @gene.tss - 1000000 & tss <= @gene.tss + 1000000')
        base_feat['ftnbh_abs_tss'] = gene_nbh.shape[0] - 1
        enh_nbh = enh_set.query('mid >= @gene.tss - 1000000 & mid <= @gene.tss + 1000000')
        if enh_nbh.empty:
            continue
        base_feat['ftnbh_abs_enh'] = enh_nbh.shape[0]
        for enhancer in enh_nbh.itertuples():
            pair_feat = dict(base_feat)
            pair_feat.update(comp_seq_features(seq[enhancer.start:enhancer.end], 'enh'))
            pair_feat['ftenh_abs_score'] = enhancer.enhancer_score
            pair_feat['ftenh_abs_elite'] = enhancer.is_elite
            pair_feat['enh_name'] = enhancer.GHid
            pair_feat['enh_start'] = enhancer.start
            pair_feat['enh_end'] = enhancer.end
            pair_feat['ftenh_abs_dist'] = enhancer.mid - gene.tss
            pair_feat['ftcons_abs_enh_mean'] = enhancer.ftcons_enh_abs_mean
            pair_feat['ftcons_abs_enh_median'] = enhancer.ftcons_enh_abs_median
            pair_feat['ftcons_abs_enh_max'] = enhancer.ftcons_enh_abs_max
            pair_feat['ftcons_abs_enh_min'] = enhancer.ftcons_enh_abs_min
            records.append(pair_feat)
    dataset = pd.DataFrame.from_records(records, index=range(len(records)))
    with lock:
        with pd.HDFStore(cargs.output, 'a', complib='blosc', complevel=9) as hdf:
            hdf.put(chrom, dataset, format='table')
            hdf.flush()
    return chrom, dataset.shape[0]


def main(args):
    """
    :param args:
    :return:
    """
    chromosomes = load_sequence(args.sequence, 'all', args.chromosomes)
    data_splits = 16  # 15 models, 1 validation set
    valid_idx = data_splits - 1
    model_indices = dict((k, col.defaultdict(list)) for k in range(0, valid_idx))
    model_indices['validation'] = col.defaultdict(list)
    with mp.Manager() as mng:
        lock = mng.Lock()
        params = [(c, args, lock) for c in chromosomes]
        with mp.Pool(args.workers) as pool:
            if args.train:
                resit = pool.imap_unordered(prepare_training_data, params)
            else:
                resit = pool.imap_unordered(prepare_testing_data, params)
            for chrom, dsize in resit:
                if dsize is None:
                    print('No data for chromosome: {}'.format(chrom))
                    continue
                if args.train:
                    subset_size = dsize // data_splits
                    for midx in range(0, valid_idx):
                        model_indices[midx][chrom].append(midx * subset_size)
                        model_indices[midx][chrom].append(midx * subset_size + subset_size)
                    assert dsize - (valid_idx * subset_size + subset_size) < data_splits,\
                        'Miscalcd index values, end up at: {} - {}'.format(valid_idx * subset_size + subset_size, dsize)
                    model_indices['validation'][chrom].append(valid_idx * subset_size)
                else:
                    pass
    if args.train:
        for midx in model_indices.keys():
            if midx == 'validation':
                continue
            outpath = args.cvindex + '_{}.pck'.format(midx)
            to_dump = model_indices[midx]
            to_dump['validation'] = model_indices['validation']
            with open(outpath, 'wb') as dump:
                pck.dump(to_dump, dump)
    return


if __name__ == '__main__':
    mp.set_start_method('forkserver')
    try:
        args = parse_command_line()
        with pd.HDFStore(args.output, 'w', complib='blosc', complevel=9):
            pass
        main(args)
    except Exception:
        trb.print_exc()
        raise
    else:
        sys.exit(0)
