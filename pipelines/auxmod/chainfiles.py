# coding=utf-8

"""
Some helper functions to process UCSC chain files
"""

import os as os

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
        qselect = chrom_select[query]
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
    assert params, 'No parameters created for chain symmetry filtering'
    return params
