# coding=utf-8

import os as os
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
    prox_buffer = []
    dist_buffer = []
    for fp in inputfiles:
        dnase, hist = os.path.basename(fp).split('_')[:2]
        cell1 = expmap[dnase]
        cell2 = expmap[hist]
        assert cell1 == cell2, 'Cell type mismatch for file {} - {} - {}'.format(fp, cell1, cell2)
        with gz.open(fp, 'rt', newline='') as infile:
            reader = csv.DictReader(infile, fieldnames=header, delimiter='\t')
            for row in reader:
                vals = list(get_items(row))
                out = vals[:3] + [cell1, '100', '.'] + vals[3:]
                if row['name'].startswith('Proximal'):
                    prox_buffer.append(tuple(out))
                else:
                    dist_buffer.append(tuple(out))
    header = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'core_start', 'core_end']
    prox_buffer = sorted(prox_buffer, key=lambda x: (x[0], int(x[1]), int(x[2]), int(x[6])))
    dist_buffer = sorted(dist_buffer, key=lambda x: (x[0], int(x[1]), int(x[2]), int(x[6])))
    with open(outputfiles[0], 'w') as out:
        _ = out.write('#' + '\t'.join(header) + '\n')
        _ = out.write('\n'.join(['\t'.join(p) for p in prox_buffer]) + '\n')

    with open(outputfiles[1], 'w') as out:
        _ = out.write('#' + '\t'.join(header) + '\n')
        _ = out.write('\n'.join(['\t'.join(d) for d in dist_buffer]) + '\n')
    return outputfiles
