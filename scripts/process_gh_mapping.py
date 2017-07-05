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
    parser.add_argument('--direction', '-d', type=str, choices=['target', 'query'], dest='direction')
    parser.add_argument('--chromosomes', '-chr', type=str, dest='chromosomes')
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


def process_target_enhancer(args):
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


def process_annotated_enhancer(params):
    """
    :param params:
    :return:
    """
    enh_file, chrom = params
    header = ['chrom', 'start', 'end', 'GHid', 'enhancer_score', 'is_elite',
              'ftcons_enh_abs_mean', 'ftcons_enh_abs_median',
              'ftcons_enh_abs_min', 'ftcons_enh_abs_max']
    enh_collect = col.defaultdict(list)
    with open(enh_file, 'r') as infile:
        rows = csv.DictReader(infile, delimiter='\t', fieldnames=header)
        for row in rows:
            if row['chrom'] == chrom:
                row['start'] = int(row['start'])
                row['end'] = int(row['end'])
                row['enhancer_score'] = float(row['enhancer_score'])
                row['ftcons_enh_abs_mean'] = float(row['ftcons_enh_abs_mean'])
                enh_collect[row['GHid']].append(row)
    enh_collect = merge_split_enhancers(enh_collect)

    ivtree = ivt.IntervalTree()
    for r in enh_collect:
        ivtree[r['start']:r['end']] = r['GHid'], r['ftcons_enh_abs_mean'], r['ftcons_enh_abs_min']
    enh_collect = sorted(enh_collect, key=lambda d: d['ftcons_enh_abs_mean'])
    blacklist = set()
    whitelist = set()
    for item in enh_collect:
        ghid = item['GHid']
        if ghid in blacklist or ghid in whitelist:
            continue
        overlaps = ivtree[item['start']:item['end']]
        if len(overlaps) == 1:
            # that is: only self overlap
            whitelist.add(ghid)
            continue
        elif len(overlaps) > 1:
            if any([o.data[0] in whitelist for o in overlaps if o.data[0] != ghid]):
                # region overlaps with a whitelist region -> blacklist
                blacklist.add(ghid)
                continue
            overlaps = [o for o in sorted(overlaps, key=lambda i: (i.data[1], i.data[2])) if o.data[0] not in blacklist]
            if overlaps[0].data[0] == ghid:
                # the query region has highest conservation
                # others can safely be blacklisted
                whitelist.add(ghid)
                [blacklist.add(o.data[0]) for o in overlaps[1:]]
            else:
                # another region is selected; could be that among
                # the remaining regions, others might also be feasible
                blacklist.add(ghid)
                whitelist.add(overlaps[0].data[0])
        else:
            raise AssertionError('No self-overlap in tree: {}'.format(item))
    enh_collect = sorted([r for r in enh_collect if r['GHid'] in whitelist], key=lambda x: (x['start'], x['end']))
    return enh_collect


def merge_split_enhancers(collector):
    """
    :param collector:
    :return:
    """
    mrg_collect = []
    for ghid, splits in collector.items():
        if len(splits) == 1:
            mrg_collect.append(splits[0])
            continue
        c = 1
        splits = sorted(splits, key=lambda d: (d['start'], d['end']))
        s, e = splits[0]['start'], splits[1]['end']
        for idx, entry in enumerate(splits[:-1]):
            if splits[idx+1]['start'] <= e + 100:
                s = min(s, splits[idx+1]['start'])
                e = max(e, splits[idx+1]['end'])
            else:
                new_enh = dict(splits[0])
                new_enh['GHid'] = new_enh['GHid'] + '-{}-{}'.format(new_enh['chrom'].strip('chr'), c)
                new_enh['start'] = s
                new_enh['end'] = e
                mrg_collect.append(new_enh)
                c += 1
                s, e = splits[idx+1]['start'], splits[idx+1]['end']
        new_enh = dict(splits[0])
        new_enh['GHid'] = new_enh['GHid'] + '-{}-{}'.format(new_enh['chrom'].strip('chr'), c)
        new_enh['start'] = s
        new_enh['end'] = e
        mrg_collect.append(new_enh)
    mrg_collect = [m for m in mrg_collect if m['end'] - m['start'] > 49]
    return mrg_collect


def process_query_enhancer(args):
    """
    :param args:
    :return:
    """
    with open(args.chromosomes, 'r') as infile:
        chroms = [l.split()[0].strip() for l in infile.readlines()]
    header = ['chrom', 'start', 'end', 'GHid', 'enhancer_score', 'is_elite',
              'ftcons_enh_abs_mean', 'ftcons_enh_abs_median',
              'ftcons_enh_abs_min', 'ftcons_enh_abs_max']
    params = [(args.inputfile, c) for c in chroms]
    with mp.Pool(args.workers) as pool:
        res = pool.imap_unordered(process_annotated_enhancer, params)
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
        if args.direction == 'target':
            process_target_enhancer(args)
        else:
            process_query_enhancer(args)
    except Exception as err:
        trb.print_exc()
        raise err
    else:
        sys.exit(0)
