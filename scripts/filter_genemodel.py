#!/usr/bin/env python3
# coding=utf-8

import csv as csv
import gzip as gz
import sys as sys
import traceback as trb
import argparse as argp
import operator as op
import io as io


def det_open_mode(fp, read=True):
    """
    :param fp:
    :param read:
    :return:
    """
    if fp.endswith('.gz'):
        if read:
            opn, mode = gz.open, 'rt'
        else:
            opn, mode = gz.open, 'wt'
    else:
        if read:
            opn, mode = open, 'r'
        else:
            opn, mode = open, 'w'
    return opn, mode


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, dest='input')
    parser.add_argument('--subset', '-sub', type=str, choices=['protein_coding'], dest='subset')
    parser.add_argument('--gene-size', '-gs', type=int, default=0, dest='genesize')
    parser.add_argument('--gene-out', '-go', type=str, dest='geneout')
    parser.add_argument('--trans-out', '-to', type=str, dest='transout')
    parser.add_argument('--map-out', '-mo', type=str, dest='mapout')
    parser.add_argument('--input-type', '-ity', type=str, choices=['gencode', 'ensucsc', 'ensbta'], dest='inputtype')
    args = parser.parse_args()
    return args


def subset_selector(intype):
    """
    :param intype:
    :return:
    """
    selector = {'gencode': 'gene_type', 'ensucsc': 'feature_type',
                'ensbta': 'source'}
    return selector[intype]


def feature_selector(intype):
    """
    :param intype:
    :return:
    """
    selector = {'gencode': 'feature', 'ensucsc': 'feature', 'ensbta': 'type'}
    return selector[intype]


def gene_column_selector(intype):
    """
    :param intype:
    :return:
    """
    selector = {'gencode': op.itemgetter(*('#chrom', 'start', 'end', 'gene_id', 'score', 'strand', 'gene_name')),
                'ensucsc': op.itemgetter(*('#chrom', 'start', 'end', 'name', 'score', 'strand', 'gene_name')),
                'ensbta': op.itemgetter(*('#seqid', 'start', 'end', 'ID', 'score', 'strand', 'Gene_name'))}
    return selector[intype]


def trans_column_selector(intype):
    """
    :param intype:
    :return:
    """
    selector = {'gencode': op.itemgetter(*('#chrom', 'start', 'end', 'transcript_id', 'score', 'strand', 'transcript_name', 'gene_id')),
                'ensucsc': op.itemgetter(*('#chrom', 'start', 'end', 'name', 'score', 'strand', 'gene_name', 'name2')),
                'ensbta': op.itemgetter(*('#seqid', 'start', 'end', 'ID', 'score', 'strand', 'Transcript_name', 'Parent'))}
    return selector[intype]


def select_protein_coding_subset(args):
    """
    :param args:
    :return:
    """
    opn, mode = det_open_mode(args.input, True)
    subset_entry = subset_selector(args.inputtype)
    subset = args.subset
    feature_entry = feature_selector(args.inputtype)
    gene_columns = gene_column_selector(args.inputtype)
    trans_columns = trans_column_selector(args.inputtype)
    genes = []
    transcripts = []
    use_genes = set()
    with opn(args.input, mode) as infile:
        rows = csv.DictReader(infile, delimiter='\t')
        for r in rows:
            if subset in r[subset_entry]:
                feat = r[feature_entry]
                if feat == 'gene':
                    if int(r['end']) - int(r['start']) < args.genesize:
                        continue
                    this_gene = list(gene_columns(r))
                    use_genes.add(this_gene[3])
                    genes.append(this_gene)
                elif feat == 'transcript':
                    this_trans = list(trans_columns(r))
                    transcripts.append(this_trans)
                else:
                    continue
    transcripts = [t for t in transcripts if t[-1] in use_genes]
    assert genes, 'No genes selected for protein coding subset'
    assert transcripts, 'No transcripts selected for protein coding subset'
    mapping = io.StringIO()
    for t in transcripts:
        mapping.write('{}\t{}\n'.format(t[3], t[-1]))
    opn, mode = det_open_mode(args.mapout, False)
    with opn(args.mapout, mode) as outf:
        _ = outf.write(mapping.getvalue())
    gene_header = ['#chrom', 'start', 'end', 'name', 'score', 'strand', 'symbol']
    opn, mode = det_open_mode(args.geneout, False)
    with opn(args.geneout, mode) as outf:
        writer = csv.writer(outf, delimiter='\t')
        writer.writerow(gene_header)
        writer.writerows(genes)
    trans_header = ['#chrom', 'start', 'end', 'name', 'score', 'strand', 'transcript_name', 'gene_id']
    opn, mode = det_open_mode(args.transout, False)
    with opn(args.transout, mode) as outf:
        writer = csv.writer(outf, delimiter='\t')
        writer.writerow(trans_header)
        writer.writerows(transcripts)
    return


def run_genemodel_filter():
    """
    :return:
    """
    args = parse_command_line()
    if args.subset == 'protein_coding':
        _ = select_protein_coding_subset(args)
    else:
        raise NotImplementedError('Unknown subset: {}'.format(args.subset))
    return 0


if __name__ == '__main__':
    try:
        status = run_genemodel_filter()
        sys.exit(status)
    except Exception as err:
        trb.print_exc()
        sys.stderr.write('\nError: {}\n'.format(err))
        sys.exit(1)
