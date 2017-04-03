#!/usr/bin/env python3
# coding=utf-8

import csv as csv
import gzip as gz
import sys as sys
import traceback as trb
import argparse as argp
import operator as op
import collections as col
import difflib as diffl
import io as io
import re as re
import json as js


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
    parser.add_argument('--ccds-only', action='store_true', default=False, dest='ccdsonly',
                        help='From GENCODE: "member of the consensus CDS gene set, '
                             'confirming coding regions between ENSEMBL, UCSC, NCBI and HAVANA."')
    parser.add_argument('--chrom', '-c', type=str, default="(chr)?[0-9XYZW][0-9AB]?$", dest='chrom')
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


def subset_value_selector(intype, subset):
    """
    The official GENCODE sequence files for the protein coding subset contain:
    - protein_coding, nonsense_mediated_decay, non_stop_decay,
    - IG_*_gene, TR_*_gene, polymorphic_pseudogene
    Using a simplified version here to come closer to the other species that
    do not have such a detailed annotation
    :param intype:
    :param subset:
    :return:
    """
    selector = {'gencode': {'protein_coding': ['protein_coding', 'IG_C_gene', 'IG_D_gene',
                                               'IG_J_gene', 'IG_LV_gene', 'IG_V_gene',
                                               'TR_C_gene', 'TR_J_gene', 'TR_V_gene',
                                               'TR_D_gene', 'polymorphic_pseudogene']},
                'ensucsc': {'protein_coding': ['protein_coding']},
                'ensbta': {'protein_coding': ['protein_coding']}}
    return selector[intype][subset]


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
    selector = {'gencode': op.itemgetter(*('#chrom', 'start', 'end', 'transcript_id',
                                           'score', 'strand', 'transcript_name', 'gene_id')),
                'ensucsc': op.itemgetter(*('#chrom', 'start', 'end', 'name', 'score',
                                           'strand', 'gene_name', 'name2')),
                'ensbta': op.itemgetter(*('#seqid', 'start', 'end', 'ID', 'score',
                                          'strand', 'Transcript_name', 'Parent'))}
    return selector[intype]


def exon_column_selector(intype):
    """
    :param intype:
    :return:
    """
    selector = {'gencode': op.itemgetter(*('#chrom', 'start', 'end', 'strand', 'transcript_id', 'exon_id')),
                'ensbta': op.itemgetter(*('#seqid', 'start', 'end', 'strand', 'Parent', 'ID')),
                'ensucsc': op.itemgetter(*('#chrom', 'exonStarts', 'exonEnds'))}
    return selector[intype]


def non_na_selector(intype):
    """
    :param intype:
    :return:
    """
    selector = {'gencode': op.itemgetter(*('gene_id', 'gene_name')),
                'ensucsc': op.itemgetter(*('name', 'gene_name')),
                'ensbta': op.itemgetter(*('ID', 'Gene_name'))}
    return selector[intype]


def select_protein_coding_subset(args):
    """
    :param args:
    :return:
    """
    opn, mode = det_open_mode(args.input, True)
    subset_entry = subset_selector(args.inputtype)
    feature_entry = feature_selector(args.inputtype)
    gene_columns = gene_column_selector(args.inputtype)
    trans_columns = trans_column_selector(args.inputtype)
    exon_columns = exon_column_selector(args.inputtype)
    nonna_values = non_na_selector(args.inputtype)
    subset = subset_value_selector(args.inputtype, args.subset)
    genes = []
    transcripts = []
    use_genes = set()
    basic_set = args.inputtype == 'gencode' and args.ccdsonly
    chrom_filter = re.compile(args.chrom.strip('"'))
    filter_stats = {'biotype': col.Counter(), 'feattype': col.Counter(),
                    'featsize': col.Counter(), 'featsize_names': [],
                    'featloc': col.Counter(), 'featloc_names': [],
                    'featna': col.Counter(), 'featna_names': [],
                    'non_ccds': col.Counter(), 'non_ccds_names': []}
    exon_collector = col.defaultdict(list)
    split_exons = args.inputtype == 'ensucsc'
    with opn(args.input, mode) as infile:
        rows = csv.DictReader(infile, delimiter='\t')
        for r in rows:
            # this weird construct: in bos taurus annotation, the
            # protein coding subset is specified as Intact_Ensembl_protein_coding
            # not sure who came up with that, but I don't like him...
            if r[subset_entry] in subset or any([s in r[subset_entry] for s in subset]):
                feat = r[feature_entry]
                if feat == 'gene':
                    not_na = nonna_values(r)
                    this_gene = list(gene_columns(r))
                    this_gene[3] = this_gene[3].split('.')[0]
                    gene_length = int(this_gene[2]) - int(this_gene[1])
                    if gene_length < args.genesize:
                        filter_stats['featsize']['gene'] += 1
                        filter_stats['featsize_names'].append([this_gene[3], gene_length])
                        continue
                    if chrom_filter.match(this_gene[0]) is None:
                        filter_stats['featloc']['gene'] += 1
                        filter_stats['featloc_names'].append([this_gene[3], this_gene[0]])
                        continue
                    if any([v == 'n/a' for v in not_na]):
                        filter_stats['featna']['gene'] += 1
                        filter_stats['featna_names'].append([this_gene[3]] + list(not_na))
                        continue
                    use_genes.add(this_gene[3])
                    genes.append(this_gene)
                elif feat == 'transcript':
                    this_trans = list(trans_columns(r))
                    this_trans[3] = this_trans[3].split('.')[0]
                    this_trans[-1] = this_trans[-1].split('.')[0]
                    if basic_set:
                        if r['tag'] != 'CCDS':
                            filter_stats['non_ccds']['transcript'] += 1
                            filter_stats['non_ccds_names'].append([this_trans[3], r['tag']])
                            continue
                    if chrom_filter.match(this_trans[0]) is None:
                        filter_stats['featloc']['transcript'] += 1
                        filter_stats['featloc_names'].append([this_trans[3], this_trans[0]])
                        continue
                    transcripts.append(this_trans)
                    if split_exons:
                        all_exons = list(exon_columns(r))
                        all_starts = all_exons[1].strip(',').split(',')
                        all_ends = all_exons[2].strip(',').split(',')
                        assert len(all_starts) == len(all_ends),\
                            'Exon start - ends differ: {} and {}'.format(all_starts, all_ends)
                        for i, (s, e) in enumerate(zip(all_starts, all_ends)):
                            name = 'Exon' + str(i)
                            exon_collector[this_trans[3]].append((int(s), int(e), this_trans[5], name))
                elif feat == 'exon':
                    this_exon = list(exon_columns(r))
                    if chrom_filter.match(this_exon[0]) is None:
                        filter_stats['featloc']['exon'] += 1
                        filter_stats['featloc_names'].append([this_exon[4], this_exon[0]])
                        continue
                    strand = this_exon[3]
                    exon_id = this_exon[-1].split('.')[0]
                    transcript = this_exon[-2].split('.')[0]
                    exon_collector[transcript].append((int(this_exon[1]), int(this_exon[2]), strand, exon_id))
                else:
                    filter_stats['feattype'][feat] += 1
                    continue
            else:
                filter_stats['biotype'][r[subset_entry]] += 1
    assert len(genes) == len(use_genes), 'Lost gene information: {} vs {}'.format(len(genes), len(use_genes))
    filter_stats['raw_transcripts'] = len(transcripts)
    transcripts = [t for t in transcripts if t[-1] in use_genes]
    filter_stats['use_transcripts'] = len(transcripts)
    filter_stats['raw_genes'] = len(genes)
    genes_with_transcripts = set([t[-1] for t in transcripts])
    genes = [g for g in genes if g[3] in genes_with_transcripts]
    filter_stats['use_genes'] = len(genes)
    filter_stats['chromosomes'] = sorted(list(set([g[0] for g in genes])))
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
    trans_header = ['#chrom', 'start', 'end', 'name', 'score', 'strand', 'transcript_name',
                    'gene_id', 'exon_ids', 'exon_starts', 'exon_ends']
    for trans in transcripts:
        ex_info = exon_collector[trans[3]]
        trans_strand = trans[5]
        ex_info = sorted(ex_info)
        # This shit here: the annotation for bosTau7 shows some
        # interesting features, e.g., annotated exons for a transcript
        # on the opposite strand (Ex.: ENSBTAT00000003589 / ENSBTAE00000028999)
        # or exons outside of the transcript boundaries (Ex.: ENSBTAT00000032823 / ENSBTAE00000473587)
        # Exons that do not match the strand of the transcript will simply be dropped
        # Exons outside of the transcript boundaries will be used to adjust transcript coordinates
        ex_info = [e for e in ex_info if e[2] == trans_strand]
        trans.append(','.join([t[3] for t in ex_info]))
        starts = [t[0] for t in ex_info]
        ends = [t[1] for t in ex_info]
        try:
            assert int(trans[1]) <= min(starts), 'Start coordinates mismatch: {} vs {} ({})'.format(starts, trans[1], trans[3])
        except AssertionError:
            sys.stderr.write('\nWarning: exon start before transcript start - adjusting transcript coordinates\n')
            trans[1] = str(min(starts))
        try:
            assert max(ends) <= int(trans[2]), 'End coordinates mismatch: {} vs {} ({})'.format(ends, trans[2], trans[3])
        except AssertionError:
            sys.stderr.write('\nWarning: exon end after transcript end - adjusting transcript coordinates\n')
            trans[2] = str(max(ends))
        trans.append(','.join([str(s) for s in starts]))
        trans.append(','.join([str(e) for e in ends]))
    opn, mode = det_open_mode(args.transout, False)
    with opn(args.transout, mode) as outf:
        writer = csv.writer(outf, delimiter='\t')
        writer.writerow(trans_header)
        writer.writerows(transcripts)
    seqm = diffl.SequenceMatcher(a=args.geneout, b=args.transout)
    stats_out = seqm.find_longest_match(0, len(args.geneout), 0, len(args.transout))
    stats_out = args.geneout[stats_out.a:stats_out.a + stats_out.size]
    stats_out += 'filter-stats.json'
    reorder = col.defaultdict(dict)
    for k, v in filter_stats.items():
        if k.endswith('_names'):
            reorder['items'][k] = v
        else:
            reorder['counters'][k] = v
    with open(stats_out, 'w') as dump:
        js.dump(reorder, dump, indent=1, sort_keys=True)
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
