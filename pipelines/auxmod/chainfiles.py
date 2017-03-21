# coding=utf-8

"""
Some helper functions to process UCSC chain files
"""

import os as os
import json as js
import collections as col

from pipelines.auxmod.auxiliary import read_chromsizes, collect_full_paths


def make_chromosome_string(chromdir):
    """
    :param chromdir:
    :return:
    """
    chromfiles = collect_full_paths(chromdir, '*.tsv')
    select = dict()
    for chrf in chromfiles:
        fp, fn = os.path.split(chrf)
        assm, _ = fn.split('_', 1)
        sizes = read_chromsizes(chrf)
        assm_select = ','.join(list(sizes.keys()))
        select[assm] = assm_select
    assert select, 'No chromosome information read from {}'.format(chromdir)
    return select


def build_chain_filter_commands(chainfiles, chromref, outpath, cmd, jobcall):
    """
    :param chainfiles:
    :param chromref:
    :param cmd:
    :param jobcall:
    :return:
    """
    chrom_select = make_chromosome_string(chromref)
    params = []
    for chf in chainfiles:
        fp, fn = os.path.split(chf)
        target, query = fn.split('.', 1)[0].split('To')
        query = query[:1].lower() + query[1:]
        tselect = chrom_select[target]
        try:
            qselect = chrom_select[query]
        except KeyError:
            continue
        tmp = cmd.format(**{'targetchroms': tselect, 'querychroms': qselect})
        outputname = '{}_to_{}.filt.chain.gz'.format(target, query)
        outputpath = os.path.join(outpath, outputname)
        params.append([chf, outputpath, tmp, jobcall])
    assert params, 'No system calls created to filter chain files: {}'.format(chainfiles)
    return params


def build_symm_filter_commands(chainfiles, chromref, outpath, cmd, jobcall):
    """
    :return:
    """
    chromfiles = collect_full_paths(chromref, '*.tsv')
    assert chromfiles, 'No chromosome files found at location: {}'.format(chromref)
    assm_chrom = dict()
    for chrf in chromfiles:
        assm = os.path.basename(chrf).split('_')[0]
        sizes = read_chromsizes(chrf)
        assm_chrom[assm] = list(sizes.keys())
    params = []
    for chf in chainfiles:
        fn = os.path.basename(chf)
        target, query = fn.split('.', 1)[0].split('_to_')
        chroms = assm_chrom[query]
        for c in chroms:
            outname = '{}_to_{}.{}.symmap.tsv.gz'.format(target, query, c)
            outfull = os.path.join(outpath, outname)
            tmp = cmd.format(**{'chrom': c})
            params.append([chf, outfull, tmp, jobcall])
    if len(chainfiles) > 0:
        assert params, 'No parameters created for chain symmetry filtering'
    return params


def check_lifted_blocks(inputfiles, outputfile):
    """
    :param inputfiles:
    :param outputfile:
    :return:
    """
    if 'lifted' in inputfiles[0]:
        liftfile = inputfiles[0]
        origfile = inputfiles[1]
    else:
        liftfile = inputfiles[1]
        origfile = inputfiles[0]
    cov = col.Counter()
    last_l = []
    last_b = []
    with open(origfile, 'r') as orig:
        with open(liftfile, 'r') as lift:
            while 1:
                if last_b:
                    bl = last_b
                    last_b = []
                else:
                    bl = _parse_block(orig.readline())
                if last_l:
                    ll = last_l
                    last_l = []
                else:
                    ll = _parse_block(lift.readline())
                if bl is None and ll is None:
                    break
                if bl is None:
                    cov['lift_cov'] += ll[3] - ll[2]
                if ll is None:
                    cov['base_cov'] += bl[3] - bl[2]
                if bl[0] == ll[0]:
                    cov['base_cov'] += bl[3] - bl[2]
                    cov['lift_cov'] += ll[3] - ll[2]
                    cov['shared_cov'] += max(min(bl[3], ll[3]) - max(bl[2], ll[2]), 0)
                elif bl[0] > ll[0]:
                    last_b = bl
                    cov['lift_cov'] += ll[3] - ll[2]
                elif bl[0] < ll[0]:
                    last_l = ll
                    cov['base_cov'] += bl[3] - bl[2]
                else:
                    raise RuntimeError('Unconsidered situation: {} and {}'.format(bl, ll))
    with open(outputfile, 'w') as dump:
        js.dump(dump, cov, indent=1, sort_keys=True)
    return outputfile


def _parse_block(line):
    """
    :param line:
    :return:
    """
    if not line.strip():
        return None
    parts = line.strip().split()
    return int(parts[3]), parts[0], int(parts[1]), int(parts[2])
