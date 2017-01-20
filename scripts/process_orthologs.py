#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import gzip as gz
import argparse as argp

import pandas as pd
import numpy as np


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, dest='input')
    parser.add_argument('--model1', '-m1', type=str, dest='model1')
    parser.add_argument('--model2', '-m2', type=str, dest='model2')
    parser.add_argument('--subset1', '-s1', type=str, dest='subset1', default='')
    parser.add_argument('--subset2', '-s2', type=str, dest='subset2', default='')
    parser.add_argument('--subset-name', '-sn', type=str, dest='subsetname', default='')
    parser.add_argument('--output', '-o', type=str, dest='output')
    args = parser.parse_args()
    return args


def get_genemodel_info():
    """
    :return:
    """
    md = {'Ensembl75_liftOver_Btau_4.6.1_genes.bed.gz': {'species': 'cow', 'name_col': 3, 'symbol_col': 9},
          'cfa_canFam3_ensembl_v81.bed.gz': {'species': 'dog', 'name_col': 3, 'symbol_col': 8},
          'gencode.v19.annotation.bed.gz': {'species': 'human', 'name_col': 15, 'symbol_col': 3},
          'gencode.vM1.annotation.bed.gz': {'species': 'mouse', 'name_col': 9, 'symbol_col': 3},
          'ssc_susScr2_ensembl_v63.bed.gz': {'species': 'pig', 'name_col': 3, 'symbol_col': 8}}
    return md


def select_from_genemodel(ortho, spec1, model1, spec2, model2):
    """
    :param ortho:
    :param spec1:
    :param model1:
    :param spec2:
    :param model2:
    :return:
    """
    s1_name, s1_symbol = '{}_gene_name'.format(spec1), '{}_gene_symbol'.format(spec1)
    s2_name, s2_symbol = '{}_gene_name'.format(spec2), '{}_gene_symbol'.format(spec2)
    spec1_col = '{}_ensembl_gene'.format(spec1)
    spec2_col = '{}_ensembl_gene'.format(spec2)
    known_genes = ortho.loc[(ortho[spec1_col].isin(model1[s1_name]) &
                             ortho[spec2_col].isin(model2[s2_name])),
                            [spec1_col, spec2_col]]
    assert not known_genes.empty, 'No orthologs identifiable in gene models'
    known_genes.columns = [s1_name, s2_name]
    # for some reason, there are a handful of duplicates
    known_genes.drop_duplicates(inplace=True)
    # adding gene symbols
    m1_sym = model1.loc[(model1[s1_name].isin(known_genes[s1_name])), [s1_name, s1_symbol]]
    m1_sym.drop_duplicates(inplace=True)
    known_genes = known_genes.merge(m1_sym, on=s1_name)

    m2_sym = model2.loc[(model2[s2_name].isin(known_genes[s2_name])), [s2_name, s2_symbol]]
    m2_sym.drop_duplicates(inplace=True)
    known_genes = known_genes.merge(m2_sym, on=s2_name)
    return known_genes


def limit_to_subset(ortho, spec1, subset1, spec2, subset2):
    """
    :param ortho:
    :param spec1:
    :param subset1:
    :param spec2:
    :param subset2:
    :return:
    """
    s1_name = '{}_gene_name'.format(spec1)
    s2_name = '{}_gene_name'.format(spec2)
    subset = ortho.loc[(ortho[s1_name].isin(subset1.name) & ortho[s2_name].isin(subset2.name)), :]
    return subset


def shorten_gene_id(gid):
    try:
        s = gid.split('.')[0]
    except AttributeError:
        s = gid
    return s


def read_gene_models(path):
    """
    :param path:
    :return:
    """
    datatypes = dict()
    with gz.open(path, 'rt') as model:
        header = model.readline().split('\t')
        for idx, h in enumerate(header):
            if 1 <= idx <= 2:
                datatypes[h] = np.int32
            else:
                datatypes[h] = str
    fp, fn = os.path.split(path)
    infos = get_genemodel_info()[fn]
    spec = infos['species']
    name_col = infos['name_col']
    symbol_col = infos['symbol_col']
    gn_col = '{}_gene_name'.format(spec)
    gs_col = '{}_gene_symbol'.format(spec)
    df = pd.read_csv(path, sep='\t', header=0, na_values='n/a',
                     skip_blank_lines=True, dtype=datatypes)
    gene_id_col = df.columns[int(name_col)]
    gene_ids = df[gene_id_col].tolist()
    gene_ids = [shorten_gene_id(g) for g in gene_ids]
    df.insert(len(df.columns), gn_col, gene_ids)

    gene_sym_col = df.columns[int(symbol_col)]
    gene_symbols = df[gene_sym_col].tolist()
    df.insert(len(df.columns), gs_col, gene_symbols)
    return spec, df


def read_gene_subsets(path):
    """
    :param path:
    :return:
    """
    header = '#chrom	start	end	name	score	strand	symbol'.split('\t')
    datatypes = dict()
    for k, v in zip(header, [str, np.int32, np.int32, str, str, str, str]):
        datatypes[k] = v
    df = pd.read_csv(path, sep='\t', header=0, na_values='n/a',
                     skip_blank_lines=True, dtype=datatypes)
    gene_ids = df['name'].tolist()
    gene_ids = [shorten_gene_id(g) for g in gene_ids]
    df['name'] = gene_ids
    return df


def compute_naive_support(orthologs):
    """
    :param orthologs:
    :return:
    """
    indices = []
    for idx, cn in enumerate(orthologs.columns, start=1):
        if cn.endswith('_assert_ids'):
            indices.append(idx)
    s1, s2 = indices
    support = []
    for row in orthologs.itertuples():
        s = min(len(row[s1].split(',')), len(row[s2].split(',')))
        support.append(s)
    orthologs = orthologs.assign(support=support)
    return orthologs


def main():
    """
    :return:
    """
    args = parse_command_line()
    raw_orthologs = pd.read_csv(args.input, sep='\t', header=0, na_values='-',
                                skip_blank_lines=True, dtype=str)
    raw_orthologs = compute_naive_support(raw_orthologs)
    spec1, model1 = read_gene_models(args.model1)
    subset1, subset2 = None, None
    subset_orthologs = None
    if args.subset1:
        subset1 = read_gene_subsets(args.subset1)
    spec2, model2 = read_gene_models(args.model2)
    if args.subset2:
        subset2 = read_gene_subsets(args.subset2)
    known_orthologs = select_from_genemodel(raw_orthologs, spec1, model1, spec2, model2)
    if args.subset1 and args.subset2:
        subset_orthologs = limit_to_subset(known_orthologs, spec1, subset1, spec2, subset2)
    md = [['ortholog_file', os.path.basename(args.input)], ['species1', spec1], ['species2', spec2],
          ['model_file1', os.path.basename(args.model1)], ['model_file2', os.path.basename(args.model2)],
          ['raw_num_orthologs', str(raw_orthologs.shape[0])], ['num_entries_model1', str(model1.shape[0])],
          ['num_entries_model2', str(model2.shape[0])],
          ['subset_file1', 'N/A' if subset1 is None else os.path.basename(args.subset1)],
          ['subset_file2', 'N/A' if subset2 is None else os.path.basename(args.subset2)],
          ['num_entries_subset1', '0' if subset1 is None else str(subset1.shape[0])],
          ['num_entries_subset2', '0' if subset2 is None else str(subset2.shape[0])],
          ['known_orthologs', str(known_orthologs.shape[0])],
          ['subset_orthologs', '0' if subset1 is None else str(subset_orthologs.shape[0])],
          ['subset_name', 'N/A' if subset1 is None else args.subsetname]]
    metadata = pd.DataFrame(md, columns=['key', 'value'], dtype='object')
    with pd.HDFStore(args.output, 'w') as hdf:
        hdf.put('/metadata', metadata, format='table')
        hdf.put('/orthologs/raw', raw_orthologs, format='table')
        hdf.put('/orthologs/known', known_orthologs, format='fixed')
        if subset_orthologs is not None:
            hdf.put('/orthologs/subset/{}'.format(args.subsetname), subset_orthologs, format='fixed')
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