# coding=utf-8

import os as os
import operator as op
import gzip as gz
import csv as csv
import json as js


def split_promoter_files(inputfiles, outputfile, expmapfile, regselect):
    """
    :param inputfiles:
    :param outputfile:
    :param expmapfile:
    :param regselect:
    :return:
    """
    # names like: Proximal-Prediction-1
    header = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'core_start', 'core_end', 'rgb']
    get_items = op.itemgetter(*('chrom', 'start', 'end', 'core_start', 'core_end'))
    out_buffer = []
    expmap = js.load(open(expmapfile, 'r'))
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
                if row['name'].startswith(regselect):
                    out_buffer.append(tuple(out))

    header = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'core_start', 'core_end']
    out_buffer = sorted(out_buffer, key=lambda x: (x[0], int(x[1]), int(x[2]), int(x[6])))
    with open(outputfile, 'w') as out:
        _ = out.write('#' + '\t'.join(header) + '\n')
        _ = out.write('\n'.join(['\t'.join(d) for d in out_buffer]) + '\n')
    return outputfile


def score_promoter_regions(inputfile, outputfile):
    """
    :param inputfile:
    :param outputfile:
    :return:
    """
    out_buffer = []
    with open(inputfile, 'r') as infile:
        reader = csv.DictReader(infile, delimiter='\t',
                                fieldnames=['#chrom', 'start', 'end', 'name', 'core_start', 'core_end'])
        for row in reader:
            score = str(min(1000, len(row['name'].split('@')) * 25))
            row['itemRgb'] = '162,44,41'
            row['score'] = score
            row['name'] = row['name'].replace('--', '-').replace('@', ',')
            row['strand'] = '.'
            out_buffer.append(row)
    out_buffer = sorted(out_buffer, key=lambda d: (d['#chrom'], int(d['start']), int(d['end'])))
    out_header = ['#chrom', 'start', 'end', 'name', 'score', 'strand', 'core_start', 'core_end', 'itemRgb']
    with open(outputfile, 'w') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=out_header, delimiter='\t')
        writer.writeheader()
        writer.writerows(out_buffer)
    return outputfile
