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

