#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import multiprocessing as mp
import csv as csv

import numpy as np
import pandas as pd


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--map-1', '-m1', type=str, dest='map1')
    parser.add_argument('--map-2', '-m2', type=str, dest='map2')
    parser.add_argument('--map-ref', '-mr', type=str, choices=['target', 'query'], dest='mapref')
    parser.add_argument('--workers', '-wrk', type=int, default=2, dest='workers')
    parser.add_argument('--output', '-o', type=str, dest='output')
    args = parser.parse_args()
    return args


def compare_chrom_mask(c1, c2):
    """
    Masks are given as bool arrays with
    True/1/masked = not conserved
    False/0/unmasked = conserved
    :param c1: np.array
    :param c2: np.array
    :return:
    """
    assert c1.size == c2.size, 'Expected chromosome masks of same size: {} vs {}'.format(c1.size, c2.size)
    raw_cov_c1 = c1.size - c1.sum()
    raw_cov_c2 = c2.size - c2.sum()
    joined_cov = c1.size - (c1 | c2).sum()
    return raw_cov_c1, raw_cov_c2, joined_cov


def load_mask_groups(fpath):
    """
    :param fpath:
    :return:
    """
    with pd.HDFStore(fpath, 'r') as hdf:
        md = hdf['metadata']
        target = md.loc[md.key == 'target', 'value'][0]
        query = md.loc[md.key == 'query', 'value'][1]
        trg_groups = [g for g in hdf.keys() if g.startswith('/target/cons/mask')]
        qry_groups = [g for g in hdf.keys() if g.startswith('/query/cons/mask')]
        trg_groups = [(g.rsplit('/', 1)[-1], g) for g in trg_groups]
        qry_groups = [(g.rsplit('/', 1)[-1], g) for g in qry_groups]
    return target, trg_groups, query, qry_groups


def build_argument_list(groups1, groups2, map1, map2):
    """
    :param groups1:
    :param groups2:
    :param map1:
    :param map2:
    :return:
    """
    params = []
    for g1, g2 in zip(sorted(groups1), sorted(groups2)):
        assert g1[0] == g2[0], 'Different chromosome names detected: {} vs {}'.format(g1, g2)
        tmp = {'map1': map1, 'map2': map2, 'load_group': g1[1], 'chrom': g1[0]}
        params.append(tmp)
    return params


def process_chrom(params):
    """
    :param params:
    :return:
    """
    with pd.HDFStore(params['map1'], 'r') as hdf:
        mask1 = hdf[params['load_group']]
    with pd.HDFStore(params['map2'], 'r') as hdf:
        mask2 = hdf[params['load_group']]
    cov1, cov2, joined = compare_chrom_mask(mask1, mask2)
    return params['chrom'], mask1.size, cov1, cov2, joined


def main():
    """
    :return:
    """
    args = parse_command_line()
    trg1, trg1_groups, qry1, qry1_groups = load_mask_groups(args.map1)
    trg2, trg2_groups, qry2, qry2_groups = load_mask_groups(args.map2)
    if args.mapref == 'target':
        assm1, assm2 = trg1, trg2
        groups1, groups2 = trg1_groups, trg2_groups
    else:
        assm1, assm2 = qry1, qry2
        groups1, groups2 = qry1_groups, qry2_groups
    assert assm1 == assm2, 'Mismatch for map reference {}: {} vs {}'.format(args.mapref, qry1, qry2)
    assert len(groups1) == len(groups2), \
        'Different number of chromosomes: {} vs {}'.format(groups1, groups2)
    paramlist = build_argument_list(groups1, groups2, args.map1, args.map2)
    total_size = 0
    total_ovl = 0
    total_cov1 = 0
    total_cov2 = 0
    outbuffer = []
    with mp.Pool(args.workers) as pool:
        m1_short = os.path.basename(args.map1).split('.')[0]
        m2_short = os.path.basename(args.map2).split('.')[0]
        resit = pool.imap_unordered(process_chrom, paramlist)
        for res in resit:
            chrom, size, cov1, cov2, ovl = res
            total_size += size
            total_ovl += ovl
            total_cov1 += cov1
            total_cov2 += cov2
            row = {'reference': assm1, 'map1': m1_short, 'map2': m2_short, 'chrom': chrom,
                   'size': size, 'cons_pos_1': cov1, 'cons_pct_1': np.round(cov1 / size, 3),
                   'cons_pos_2': cov2, 'cons_pct_2': np.round(cov2 / size, 3),
                   'cons_pos_ovl': ovl, 'cons_pct_ovl': np.round(ovl / size, 3)}
            outbuffer.append(row)
    outbuffer = sorted(outbuffer, key=lambda d: d['size'], reverse=True)
    row = {'reference': assm1, 'map1': m1_short, 'map2': m2_short, 'chrom': 'genome',
           'size': total_size, 'cons_pos_1': total_cov1, 'cons_pct_1': np.round(total_cov1 / total_size, 3),
           'cons_pos_2': total_cov2, 'cons_pct_2': np.round(total_cov2 / total_size, 3),
           'cons_pos_ovl': total_ovl, 'cons_pct_ovl': np.round(total_ovl / total_size, 3)}
    outbuffer.append(row)
    with open(args.output, 'w') as out:
        writer = csv.DictWriter(out, fieldnames=['reference', 'map1', 'map2', 'chrom',
                                                 'size', 'cons_pos_1', 'cons_pct_1',
                                                 'cons_pos_2', 'cons_pct_2', 'cons_pos_ovl',
                                                 'cons_pct_ovl'], delimiter='\t')
        writer.writeheader()
        writer.writerows(outbuffer)
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
