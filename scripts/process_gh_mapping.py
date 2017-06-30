#!/usr/bin/env python
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import csv as csv
import functools as fnt
import collections as col
import multiprocessing as mp

import numpy as np
import pandas as pd
import intervaltree as ivt


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, dest='inputfile')
    parser.add_argument('--cons', '-c', type=str, dest='conservation')
    parser.add_argument('--workers', '-w', type=int, default=4, dest='workers')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')
    args = parser.parse_args()
    return args


def compute_weights(cons, enh):
    """
    :param enh:
    :param cons:
    :return:
    """
    s, e = enh['start'], enh['end']
    enh['ftcons_enh_abs_min'] = np.round(cons[s:e].min(), 2)
    enh['ftcons_enh_abs_max'] = np.round(cons[s:e].max(), 2)
    enh['ftcons_enh_abs_mean'] = np.round(cons[s:e].mean(), 2)
    enh['ftcons_enh_abs_median'] = np.round(cons[s:e].median(), 2)
    total_score = sum([float(s) for s in enh['assoc_score'].split(',')])
    enh['weight'] = 1000 - np.round(float(enh['enhancer_score']) * total_score, 2)
    return enh


def process_mapped_enhancer(params):
    """
    :param params:
    :return:
    """
    inputfile, consfile, chrom = params

    with pd.HDFStore(consfile, 'r') as hdf:
        cons_scores = hdf[chrom]
    comp_wt = fnt.partial(compute_weights, cons_scores)
    header = ['chrom', 'start', 'end', 'GHid', 'enhancer_score', 'is_elite', 'cluster_id',
              'name', 'symbol', 'assoc_score', 'enh_gene_dist']
    regions = []
    with open(inputfile, 'r', newline='') as infile:
        rows = csv.DictReader(infile, delimiter='\t', fieldnames=header)
        for r in rows:
            if r['chrom'] == chrom:
                r['start'] = int(r['start'])
                r['end'] = int(r['end'])
                l = r['end'] - r['start']
                if l < 2:
                    continue
                r = comp_wt(r)
                regions.append(comp_wt(r))
    ivtree = ivt.IntervalTree()
    for r in regions:
        ivtree[r['start']:r['end']] = r['GHid'], r['ftcons_enh_abs_mean'], r['weight']
    regions = sorted(regions, key=lambda d: d['ftcons_enh_abs_mean'])
    blacklist = set()
    whitelist = set()
    for item in regions:
        if item['GHid'] in blacklist:
            continue
        overlaps = ivtree[item['start']:item['end']]
        if len(overlaps) == 1:
            whitelist.add(overlaps.pop().data[0])
        elif len(overlaps) > 1:
            overlaps = [o for o in sorted(overlaps, key=lambda i: (i.data[1], i.data[2])) if o.data[0] not in blacklist]
            whitelist.add(overlaps[0].data[0])
            [blacklist.add(o.data[0]) for o in overlaps[1:]]
        else:
            raise AssertionError('No self-overlap in tree: {}'.format(item))
    regions = sorted([r for r in regions if r['GHid'] in whitelist], key=lambda x: (x['start'], x['end']))
    return regions


def main(args):
    """
    :param args:
    :return:
    """
    with pd.HDFStore(args.conservation, 'r') as hdf:
        chroms = [k.strip('/') for k in hdf.keys()]
    params = [(args.inputfile, args.conservation, c) for c in chroms]
    header = ['chrom', 'start', 'end', 'GHid', 'enhancer_score', 'is_elite',
              'name', 'symbol', 'assoc_score', 'ftcons_enh_abs_mean',
              'ftcons_enh_abs_median', 'ftcons_enh_abs_min', 'ftcons_enh_abs_max']
    with mp.Pool(args.workers) as pool:
        res = pool.imap_unordered(process_mapped_enhancer, params)
        outbuffer = []
        for regions in res:
            outbuffer.extend(regions)
        outbuffer = sorted(outbuffer, key=lambda d: (d['chrom'], d['start'], d['end']))
        with open(args.outputfile, 'w') as out:
            _ = out.write('#')
            writer = csv.DictWriter(out, fieldnames=header, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            writer.writerows(outbuffer)

    return


if __name__ == '__main__':
    try:
        args = parse_command_line()
        main(args)
    except Exception as err:
        trb.print_exc()
        raise err
    else:
        sys.exit(0)
