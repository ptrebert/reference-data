#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import csv as csv
import traceback as trb
import argparse as argp
import collections as col

import pandas as pd


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--genes', '-g', type=str, dest='genes')
    parser.add_argument('--mapping', '-m', type=str, dest='mapping')
    parser.add_argument('--output', '-o', type=str, dest='output')

    args = parser.parse_args()
    return args


def read_gene_transcript_map(fpath):
    """
    :param fpath:
    :return:
    """
    df = pd.read_csv(fpath, sep='\t', header=None, names=['transcript', 'gene'])
    return df


def read_gencode_annotation(fpath, tgmap):
    """
    :param fpath:
    :param tgmap:
    :return:
    """
    df = pd.read_csv(fpath, sep='\t', header=0,
                     usecols=['#chrom', 'start', 'end', 'gene_id', 'feature', 'exon_id'])
    df = df[['#chrom', 'start', 'end', 'gene_id', 'feature', 'exon_id']]
    df.columns = ['chrom', 'start', 'end', 'gene_id', 'feature', 'exon_id']
    df = df.loc[df['feature'] == 'exon', :]
    df['gene_short_id'] = df['gene_id'].str.extract('([A-Z0-9]+)', expand=False)
    df = df.loc[df['gene_short_id'].isin(tgmap['gene']), :].copy()
    df['exon_short_id'] = df['exon_id'].str.extract('([A-Z0-9]+)', expand=False)
    collector = col.defaultdict(list)
    for row in df.itertuples():
        collector[(row.gene_short_id, row.chrom)].append((row.start, row.end))
    return collector


def read_ucsc_ensembl_annotation(fpath, tgmap):
    """
    :param fpath:
    :param tgmap:
    :return:
    """
    df = pd.read_csv(fpath, sep='\t', header=0,
                     usecols=['#chrom', 'start', 'end', 'feature', 'name2', 'exonStarts', 'exonEnds'])
    df = df[['#chrom', 'start', 'end', 'feature', 'name2', 'exonStarts', 'exonEnds']]
    df.columns = ['chrom', 'start', 'end', 'feature', 'name2', 'exonStarts', 'exonEnds']
    df = df.loc[df['feature'] == 'transcript', :]
    df['gene_short_id'] = df['name2'].str.extract('([A-Z0-9]+)', expand=False)
    df = df.loc[df['gene_short_id'].isin(tgmap['gene']), :].copy()
    collector = col.defaultdict(list)
    for row in df.itertuples():
        starts = row.exonStarts.strip(',').split(',')
        ends = row.exonEnds.strip(',').split(',')
        for s, e in zip(starts, ends):
            collector[(row.gene_short_id, row.chrom)].append((s, e))
    return collector


def read_ensembl_annotation(fpath, tgmap):
    """
    :param fpath:
    :param tgmap:
    :return:
    """
    df = pd.read_csv(fpath, sep='\t', header=0,
                     usecols=['#seqid', 'start', 'end', 'type', 'Parent'])
    df = df[['#seqid', 'start', 'end', 'type', 'Parent']]
    df.columns = ['chrom', 'start', 'end', 'type', 'transcript']
    df = df.loc[df['type'] == 'exon', :]
    df = df.loc[df['transcript'].isin(tgmap['transcript']), :]
    df = df.merge(tgmap, how='outer', on='transcript', suffixes=('', ''))
    collector = col.defaultdict(list)
    for row in df.itertuples():
        collector[(row.gene, row.chrom)].append((row.start, row.end))
    return collector


def dump_annotation(annotation, output):
    """
    :param annotation:
    :param output:
    :return:
    """
    header = ['name', 'chrom', 'starts', 'ends']
    rows = []
    for (name, chrom), exons in annotation.items():
        exons = sorted(exons, key=lambda x: int(x[0]))
        starts = ','.join([str(x[0]) for x in exons])
        ends = ','.join([str(x[1]) for x in exons])
        rows.append({'name': name, 'chrom': chrom, 'starts': starts, 'ends': ends})
    rows = sorted(rows, key=lambda d: (d['chrom'], d['name']))
    with open(output, 'w') as dump:
        writer = csv.DictWriter(dump, fieldnames=header, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    return


def main():
    """
    :return:
    """
    args = parse_command_line()
    tgmap = read_gene_transcript_map(args.mapping)
    fn = os.path.basename(args.genes)
    if fn.startswith('gencode'):
        annotation = read_gencode_annotation(args.genes, tgmap)
    elif '_ensembl_' in fn:
        annotation = read_ucsc_ensembl_annotation(args.genes, tgmap)
    elif fn.startswith('Ensembl'):
        annotation = read_ensembl_annotation(args.genes, tgmap)
    else:
        raise ValueError('Unknown gene annotation file: {}'.format(fn))
    dump_annotation(annotation, args.output)
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
