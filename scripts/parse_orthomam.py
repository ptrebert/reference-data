#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import csv as csv
import traceback as trb
import argparse as argp
import xml.etree.ElementTree as ET

import pandas as pd


def read_orthomam_genes(xmlfiles):
    """
    :return:
    """
    select_species = ['Homo_sapiens', 'Mus_musculus', 'Bos_taurus', 'Canis_lupus', 'Sus_scrofa']
    species_map = dict((k, v) for k, v in zip(select_species, ['human', 'mouse', 'cow', 'dog', 'pig']))

    orthologs = []
    for xmlf in xmlfiles:
        tree = ET.parse(xmlf)
        root = tree.getroot()
        assert root.tag == 'marker', 'Unexpected XML root: {}'.format(root.tag)
        symbol = root.find('symb').text
        if symbol == 'NONAME':
            continue
        gene_desc = root.find('geneDescription')
        ens_id = gene_desc.find('geneId').text
        assert ens_id.startswith('ENSG'), 'Unexpected gene identifier: {}'.format(ens_id)
        desc = gene_desc.find('desc').text
        full_name, source = desc.split('[')
        _, acc = source.split(';Acc:')
        assert acc.startswith('HGNC:'), 'Unexpected source: {}'.format(source)
        full_name = full_name.strip().replace(' ', '_').replace(',', '_').replace('__', '_')
        record = {'human_gene_symbol': symbol, 'human_gene_name': ens_id,
                  'gene_full_name': full_name.strip(), 'HGNC_acc': acc.strip(']')}
        for cds in root.iter('speciesCDS'):
            long_species = cds.attrib['species_web_name']
            if long_species not in select_species:
                continue
            ens_species = cds.find('infoCDS').find('ensid').text
            short_species = species_map[long_species]
            if short_species == 'human':
                assert ens_species == ens_id, 'ID mismatch: {} vs {}'.format(ens_id, ens_species)
            else:
                assert ens_species.startswith('ENS'),\
                    'Unexpected gene name for species {}: {}'.format(short_species, ens_species)
                record['{}_gene_name'.format(short_species)] = ens_species
        orthologs.append(record)
    assert orthologs, 'No ortholog records read'
    orthologs = sorted(orthologs, key=lambda x: x['HGNC_acc'])
    return orthologs


def read_gene_maps(fpaths):
    """
    :param fpaths:
    :return:
    """
    genesets = dict()
    for fp in fpaths:
        genes = set()
        with open(fp, 'r') as infile:
            for line in infile:
                if not line.strip():
                    continue
                _, g = line.strip().split()
                genes.add(g)
        fn = os.path.basename(fp)
        prefix = fn.split('_')[0]
        genesets[prefix] = genes
    return genesets


def annotate_subsets(orthologs, genesets):
    """
    :param orthologs:
    :param genesets:
    :return:
    """
    for rec in orthologs:
        if rec['human_gene_name'] in genesets['hsa'] and \
            rec['mouse_gene_name'] in genesets['mmu'] and \
            rec['cow_gene_name'] in genesets['bta'] and \
            rec['dog_gene_name'] in genesets['cfa'] and \
                rec['pig_gene_name'] in genesets['ssc']:
            orth_five = 1
            orth_four = 1
        elif rec['human_gene_name'] in genesets['hsa'] and \
            rec['mouse_gene_name'] in genesets['mmu'] and \
            rec['cow_gene_name'] in genesets['bta'] and \
                rec['pig_gene_name'] in genesets['ssc']:
            orth_five = 0
            orth_four = 1
        else:
            orth_five = 0
            orth_four = 0
        rec['orth_five'] = orth_five
        rec['orth_four'] = orth_four
    return orthologs


def make_dataframe(orthologs, subgenes):
    """
    :param orthologs:
    :param subgenes:
    :return:
    """
    genemaps = os.listdir(subgenes)
    genemaps = [os.path.join(subgenes, fn) for fn in genemaps if fn.endswith('.map.tsv')]
    genesets = read_gene_maps(genemaps)
    orthologs = annotate_subsets(orthologs, genesets)
    df = pd.DataFrame.from_records(orthologs)
    md = [['raw_orthologs', str(df.shape[0])], ['five_orthologs', str(df['orth_five'].sum())],
          ['four_orthologs', str(df['orth_four'].sum())], ['subset', 'protein_coding']]
    md = pd.DataFrame(md, columns=['key', 'value'])
    return df, md


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, dest='input')
    parser.add_argument('--tsv-output', '-to', type=str, dest='tsvoutput')
    parser.add_argument('--sub-genes', '-sg', type=str, dest='subgenes')
    parser.add_argument('--hdf-output', '-ho', type=str, dest='hdfoutput')
    args = parser.parse_args()
    return args


def main():
    """
    :return:
    """
    args = parse_command_line()
    assert os.path.isdir(args.input), 'Input is not a valid directory folder: {}'.format(args.input)
    all_files = os.listdir(args.input)
    xmlfiles = [os.path.join(args.input, fn) for fn in all_files if fn.endswith('.xml')]
    ortholog_rec = read_orthomam_genes(xmlfiles)
    header = ['HGNC_acc', 'human_gene_symbol', 'gene_full_name', 'human_gene_name',
              'mouse_gene_name', 'cow_gene_name', 'dog_gene_name', 'pig_gene_name']
    with open(args.tsvoutput, 'w') as dump:
        writer = csv.DictWriter(dump, fieldnames=header, restval='n/a', delimiter='\t')
        writer.writeheader()
        writer.writerows(ortholog_rec)

    df, md = make_dataframe(ortholog_rec, args.subgenes)
    with pd.HDFStore(args.hdfoutput, 'w', complevel=9, complib='blosc') as hdf:
        hdf.put('metadata', md, format='table')
        hdf.put('orthologs', df, format='fixed')
        hdf.flush()
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