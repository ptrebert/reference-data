# coding=utf-8

import os as os
import io as io
import operator as op
import gzip as gz
import csv as csv


def merge_promoter_files(inputfiles, outputfiles, expmap):
    """
    :param inputfiles:
    :param outputfiles:
    :param expmap:
    :return:
    """
    # names like: Proximal-Prediction-1
    header = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'core_start', 'core_end', 'rgb']
    get_items = op.itemgetter(*('chrom', 'start', 'end', 'core_start', 'core_end'))
    prox_buffer = io.StringIO()
    dist_buffer = io.StringIO()
    for fp in inputfiles:
        dnase, hist = os.path.basename(fp).split('_')[:2]
        cell1 = expmap[dnase]
        cell2 = expmap[hist]
        assert cell1 == cell2, 'Cell type mismatch for file {} - {} - {}'.format(fp, cell1, cell2)
        with gz.open(fp, 'rt', newline='') as infile:
            reader = csv.DictReader(infile, fieldnames=header, delimiter='\t')
            for row in reader:
                vals = get_items(row)
                if row['name'].startswith('Proximal'):

                else:

