#!/usr/bin/env python3
# coding=utf-8

import csv as csv
import gzip as gz
import sys as sys
import argparse as argp
import operator as op
import traceback as trb
import collections as col


SHOW_PARSER_WARNING = False
# The following is (hopefully) correct: GENCODE GTF
# files are not BED-like (i.e. not 0-based half-open),
# and I assume the bovine GFF files are neither (but
# I don't have a gold standard to check against)
# Since this converter assumes BED-like output, this at least
# corrects the GENCODE GTF files
ASSUME_OFFBYONE = False
# For the Ensembl UCSC files, I assume 0-based coordinates
# since this is the usual case for data stored in Genome Browser
# tables


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
    parser.add_argument('--output', '-o', type=str, dest='output')
    parser.add_argument('--filter-size', '-fs', type=int, default=0, dest='filtersize')
    parser.add_argument('--filter-type', '-ft', type=str, default='', dest='filtertype')
    parser.add_argument('--show-warnings', '-w', action='store_true', default=False, dest='warnings')
    parser.add_argument('--off-by-one', '-obo', action='store_true', default=False, dest='obo')
    parser.add_argument('--ucsc-ensembl', action='store_true', default=False, dest='ucscensembl')
    parser.add_argument('--ensembl-source', type=str, default='', dest='ensemblsource')
    parser.add_argument('--ensembl-name', type=str, default='', dest='ensemblname')
    args = parser.parse_args()
    return args


def filter_data_entries(rows, types, size):
    """
    :param rows:
    :param types:
    :param size:
    :return:
    """
    if types:
        types = types.strip('"').split(',')
        keys, values = [], []
        for entry in types:
            k, v = entry.split(':')
            keys.append(k.strip())
            values.append(v.strip())
        values = tuple(values)
        getter = op.itemgetter(*keys)
        tmp = []
        for r in rows:
            if getter(r) == values:
                tmp.append(r)
        rows = tmp
    if size > 0:
        rows = [r for r in rows if int(r['end']) - int(r['start']) >= size]
    assert rows, 'No rows left after filtering'
    return rows


def gtf_column_key(c):
    """
    :param c:
    :return:
    """
    cols = {'chrom': 0, 'start': 1, 'end': 2, 'gene_name': 3,
            'score': 4, 'strand': 5, 'frame': 6, 'feature': 7,
            'source': 8}
    return cols.get(c, 9)


def gff3_column_key(c):
    """
    :param c:
    :return:
    """
    cols = {'seqid': 0, 'start': 1, 'end': 2, 'ID': 3,
            'score': 4, 'strand': 5, 'phase': 6, 'type': 7,
            'source': 8}
    return cols.get(c, 9)


def ensembl_column_key(c):
    """
    :param c:
    :return:
    """
    cols = {'chrom': 0, 'start': 1, 'end': 2, 'name': 3,
            'score': 4, 'strand': 5, 'feature': 6,
            'feature_type': 7, 'gene_name': 8}
    return cols.get(c, 9)


def read_ensembl_genes(fpath, source, name):
    """
    :param fpath:
    :param source: content of ensemblSource file; contains transcript type info, e.g., protein coding
    :param name: content of ensemblToGeneName file; contains gene symbol/name info
    :return:
    """
    fieldnames = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
                  'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat',
                  'cdsEndStat', 'exonFrames']
    ensembl_genes = col.defaultdict(list)
    ensembl_rows = []
    opn, mode = det_open_mode(fpath, True)
    with opn(fpath, mode) as infile:
        rows = csv.DictReader(infile, fieldnames=fieldnames, delimiter='\t')
        for r in rows:
            assert r['name'] in source, 'Unknown transcript type: {}'.format(r['name'])
            try:
                gn = name[r['name']]
            except KeyError:
                gn = 'n/a'
            r['feature_type'] = source[r['name']]
            r['gene_name'] = gn
            r['feature'] = 'transcript'
            r['start'] = r['txStart']
            r['end'] = r['txEnd']
            del r['txStart']
            del r['txEnd']
            del r['bin']
            if r['name2'].startswith('ENS'):
                ensembl_genes[r['name2']].append(r)
            ensembl_rows.append(r)
    ensembl_genes = normalize_ensembl_genes(ensembl_genes)
    ensembl_rows.extend(ensembl_genes)
    fields = set()
    for entry in ensembl_rows:
        fields = fields.union(set(list(entry.keys())))
    return fields, ensembl_rows


def normalize_ensembl_genes(ensgenes):
    """
    :param ensgenes:
    :return:
    """
    ens_genes = []
    for ensid, trans in ensgenes.items():
        strands = set(t['strand'] for t in trans)
        assert len(strands) == 1, 'Switching strands: {}'.format(trans)
        strand = strands.pop()
        feat_types = set([t['feature_type'] for t in trans])
        assert len(feat_types) == 1, 'Mixed feature types: {}'.format(trans)
        feat_type = feat_types.pop()
        gene_names = set([t['gene_name'] for t in trans])
        assert len(gene_names) == 1, 'Different gene names: {}'.format(trans)
        gene_name = gene_names.pop()
        chrom = trans[0]['chrom']
        left = str(min([int(t['start']) for t in trans]))
        right = str(max([int(t['end']) for t in trans]))
        infos = {'chrom': chrom, 'start': left, 'end': right, 'strand': strand,
                 'feature_type': feat_type, 'feature': 'gene', 'score': '0',
                 'name': ensid, 'gene_name': gene_name}
        ens_genes.append(infos)
    return ens_genes


def read_value_mapping(fpath):
    """
    :param fpath:
    :return:
    """
    valmap = dict()
    opn, mode = det_open_mode(fpath, True)
    with opn(fpath, mode) as infile:
        for line in infile:
            if not line.strip():
                continue
            try:
                k, v = line.strip().split(maxsplit=1)
                # another jump through a hoop for those funny people
                # out there who think whitespaces should be part of
                # a gene name...
                v = v.replace(' ', '-')
            except ValueError as ve:
                raise ValueError('Cannot process line: {}\nError: {}'.format(line.strip(), str(ve)))
            assert k not in valmap, 'Duplicate key {} in mapping file {}'.format(k, fpath)
            valmap[k] = v
    return valmap


def read_gtf_line(line):
    """
    :param line:
    :return:
    """
    cols = line.strip().split('\t')
    if ASSUME_OFFBYONE:
        cols[3] = str(int(cols[3]) - 1)
    known_fields = set()
    fields = {'chrom': cols[0], 'source': cols[1], 'feature': cols[2],
              'start': cols[3], 'end': cols[4], 'score': cols[5],
              'strand': cols[6], 'frame': cols[7]}
    known_fields.update(fields.keys())
    attrlist = cols[8]
    attributes = dict()
    for attr in attrlist.split(';'):
        if not attr.strip():
            continue
        k, v = attr.strip().split()
        attributes[k] = v.strip('"')
        known_fields.add(k)
    fields.update(attributes)
    return fields, known_fields


def read_gff3_line(line):
    """
    Try to parse a line in the goddamn GFF3 format
    w/o loosing any attributes...
    :param line:
    :return:
    """
    cols = line.strip().split('\t')
    if ASSUME_OFFBYONE:
        cols[3] = str(int(cols[3]) - 1)
    known_fields = set()
    fields = {'seqid': cols[0], 'source': cols[1], 'type': cols[2],
              'start': cols[3], 'end': cols[4], 'score': cols[5],
              'strand': cols[6], 'phase': cols[7]}
    known_fields.update(fields.keys())
    attrlist = cols[8]
    attributes = dict()
    for attr in attrlist.split(';'):
        if not attr.strip():
            continue
        k, v = attr.strip().split('=')
        if k.lower() == 'dbxref':
            try:
                subkey, subvalue = v.split(':')
            except ValueError:
                if SHOW_PARSER_WARNING:
                    sys.stderr.write('\nWarning: skipping Dbxref value {} - no key! Line: {}'.format(v, line.strip()))
                continue
            assert subkey not in attributes, 'Sub-key already in attributes list: {} in line {}'.format(subkey, line)
            attributes[subkey] = subvalue.strip()
            known_fields.add(subkey)
            continue
        elif ',' in v:
            raise ValueError('List of values for key {}: {} in line {}'.format(k, v, line))
        else:
            # who knows what crazy stuff could be here...
            pass
        attributes[k] = v.strip()
        known_fields.add(k)
    fields.update(attributes)
    return fields, known_fields


def convert_gtf_file(args):
    """
    :param args:
    :return:
    """
    rowbuffer = []
    read_gtf = read_gtf_line
    file_keys = set()
    opn, mode = det_open_mode(args.input, True)
    with opn(args.input, mode) as infile:
        for line in infile:
            if not line.strip() or line.startswith('#'):
                continue
            data, keys = read_gtf(line)
            rowbuffer.append(data)
            file_keys.update(keys)
    header = sorted(file_keys, key=gtf_column_key)
    rowbuffer = filter_data_entries(rowbuffer, args.filtertype, args.filtersize)
    opn, mode = det_open_mode(args.output, False)
    with opn(args.output, mode) as outf:
        _ = outf.write('#')
        writer = csv.DictWriter(outf, delimiter='\t', fieldnames=header,
                                restval='n/a')
        writer.writeheader()
        rowbuffer = sorted(rowbuffer, key=lambda x: (x['chrom'], int(x['start']), int(x['end'])))
        writer.writerows(rowbuffer)
    return


def convert_gff3_file(args):
    """
    :param args:
    :return:
    """
    rowbuffer = []
    read_gff3 = read_gff3_line
    file_keys = set()
    opn, mode = det_open_mode(args.input, True)
    with opn(args.input, mode) as infile:
        for line in infile:
            if not line.strip() or line.startswith('#'):
                continue
            data, keys = read_gff3(line)
            rowbuffer.append(data)
            file_keys.update(keys)
    header = sorted(file_keys, key=gff3_column_key)
    rowbuffer = filter_data_entries(rowbuffer, args.filtertype, args.filtersize)
    opn, mode = det_open_mode(args.output, False)
    with opn(args.output, mode) as outf:
        _ = outf.write('#')
        writer = csv.DictWriter(outf, delimiter='\t', fieldnames=header,
                                restval='n/a')
        writer.writeheader()
        rowbuffer = sorted(rowbuffer, key=lambda x: (x['seqid'], int(x['start']), int(x['end'])))
        writer.writerows(rowbuffer)
    return


def convert_ensembl_annotation(args):
    """
    :param args:
    :return:
    """
    ens_source = read_value_mapping(args.ensemblsource)
    ens_name = read_value_mapping(args.ensemblname)
    fields, rows = read_ensembl_genes(args.input, ens_source, ens_name)
    header = sorted(fields, key=ensembl_column_key)
    rowbuffer = filter_data_entries(rows, args.filtertype, args.filtersize)
    opn, mode = det_open_mode(args.output, False)
    with opn(args.output, mode) as outf:
        _ = outf.write('#')
        writer = csv.DictWriter(outf, delimiter='\t', fieldnames=header,
                                restval='n/a')
        writer.writeheader()
        rowbuffer = sorted(rowbuffer, key=lambda x: (x['chrom'], int(x['start']), int(x['end']), x['feature']))
        writer.writerows(rowbuffer)
    return


def run_genemodel_conversion():
    """
    :return:
    """
    args = parse_command_line()
    global SHOW_PARSER_WARNING
    SHOW_PARSER_WARNING = args.warnings
    global ASSUME_OFFBYONE
    ASSUME_OFFBYONE = args.obo
    if '.gtf' in args.input:
        convert_gtf_file(args)
    elif '.gff' in args.input:
        convert_gff3_file(args)
    elif args.ucscensembl:
        convert_ensembl_annotation(args)
    else:
        raise NotImplementedError
    return 0


if __name__ == '__main__':
    try:
        status = run_genemodel_conversion()
        sys.exit(status)
    except Exception as err:
        trb.print_exc()
        sys.stderr.write('\nError: {}\n'.format(err))
        sys.exit(1)
