#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import collections as col
import traceback as trb
import argparse as argp
import fnmatch as fnm
import csv as csv

import pandas as pd
import numpy as np


SPECIES_MAP = {'human': 'hsa',
               'mouse': 'mmu',
               'cow': 'bta',
               'chicken': 'gga',
               'dog': 'cfa',
               'pig': 'ssc'}


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, nargs='+', dest='inputfiles')
    parser.add_argument('--cache-file', '-cf', type=str, dest='cache', default='no_cache')
    parser.add_argument('--annotation', '-ann', type=str, dest='annotation')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')
    parser.add_argument('--strategy', '-s', type=str, choices=['top-down', 'bottom-up'],
                        dest='strategy', default='bottom-up')
    parser.add_argument('--species-table', '-spt', type=str, default='', dest='speciestable')
    args = parser.parse_args()
    return args


def read_species_table(fpath, gene_annot):
    """
    :param fpath:
    :return:
    """
    gene_models = os.listdir(gene_annot)
    lut = {}
    with open(fpath, 'r', newline='') as table:
        rows = csv.DictReader(table, delimiter='\t')
        for row in rows:
            modelfile = fnm.filter(gene_models, '{}_*.genes.bed.gz'.format(row['kegg_org_code']))
            if len(modelfile) == 0:
                # skip this species
                continue
            assert len(modelfile) == 1, \
                'Multiple gene models for species {} / {}'.format(row['common_name', row['kegg_org_code']])
            lut[row['common_name']] = {'tax_id': int(row['taxon_id']), 'code': row['kegg_org_code'],
                                       'genes': os.path.join(gene_annot, modelfile[0])}
            lut[row['kegg_org_code']] = {'tax_id': int(row['taxon_id']), 'name': row['common_name'],
                                         'genes': os.path.join(gene_annot, modelfile[0])}
            lut[int(row['taxon_id'])] = {'name': row['common_name'], 'code': row['kegg_org_code'],
                                         'genes': os.path.join(gene_annot, modelfile[0])}
    return lut


def odb_header(which):
    """
    :param which:
    :return:
    """
    all_headers = {'genes': ['odb_gene_id', 'tax_id', 'protein_id', 'uniprot_id',
                             'gene_name', 'ncbi_gid', 'desc'],
                   'OGs': ['og_id', 'level', 'og_name'],
                   'OG2genes': ['og_id', 'odb_gene_id']}
    return all_headers[which]


def merge_tables(genes, ogroups, mapping, spec_annot):
    """
    :param genes:
    :param ogroups:
    :param mapping:
    :return:
    """
    taxons = set()
    for k, v in spec_annot.items():
        try:
            taxons.add(v['tax_id'])
        except KeyError:
            continue
    taxons = sorted(taxons)
    indexer = genes['tax_id'].isin(taxons)
    species_genes = (genes.loc[indexer, :]).copy(deep=True)
    assert not species_genes.empty, 'No genes selected based on taxon ids: {}'.format(taxons)
    indexer = mapping['odb_gene_id'].isin(species_genes['odb_gene_id'])
    species_map = (mapping.loc[indexer, :]).copy(deep=True)
    indexer = ogroups['og_id'].isin(species_map['og_id'])
    species_og = (ogroups.loc[indexer, :]).copy(deep=True)

    merged = species_og.merge(species_map, on='og_id', how='outer')
    merged = merged.merge(species_genes, on='odb_gene_id', how='outer')
    merged['clade_id'] = merged['og_id'].str.slice(start=0, stop=9)
    return merged


def raw_process_input(args, spec_table):
    """
    :param args:
    :return:
    """
    genes_tab, og_tab, og2genes_tab = None, None, None
    for inpf in args.inputfiles:
        assert 'odb9' in inpf, 'Unexpected file: {}'.format(inpf)
        ttype = os.path.basename(inpf).split('.')[0].split('_')[1]
        if ttype == 'genes':
            assert genes_tab is None, 'Duplicate gene: {}'.format(inpf)
            genes_tab = pd.read_csv(inpf, delimiter='\t', header=None,
                                    usecols=['odb_gene_id', 'tax_id', 'gene_name'],
                                    names=odb_header(ttype), skip_blank_lines=True,
                                    dtype={'tax_id': np.int32, 'gene_name': str, 'odb_gene_id': str})
        elif ttype == 'OGs':
            assert og_tab is None, 'Duplicate OG table: {}'.format(inpf)
            og_tab = pd.read_csv(inpf, delimiter='\t', header=None,
                                 names=odb_header(ttype), skip_blank_lines=True,
                                 usecols=['og_id', 'og_name'], dtype={'og_id': str, 'og_name': str})
        elif ttype == 'OG2genes':
            assert og2genes_tab is None, 'Duplicate map table: {}'.format(inpf)
            og2genes_tab = pd.read_csv(inpf, delimiter='\t', header=None,
                                       names=odb_header(ttype), skip_blank_lines=True,
                                       dtype={'og_id': str, 'odb_gene_id': str})
        else:
            raise ValueError('Unexpected table type: {}'.format(inpf))
    merged = merge_tables(genes_tab, og_tab, og2genes_tab, spec_table)
    return merged


def read_gene_model(fpath, species):
    """
    :param fpath:
    :param species:
    :return:
    """
    genes = pd.read_csv(fpath, skip_blank_lines=True, header=0,
                        na_values='n/a', sep='\t', dtype=str,
                        usecols=['#chrom', 'name', 'symbol'])
    new_cols = []
    for c in genes.columns:
        if c == '#chrom':
            new_cols.append('chrom')
        else:
            new_cols.append('{}_{}'.format(species, c))
    genes.columns = new_cols
    return genes


def select_ortholog_pairs(dataset, a_genes, b_genes, a_name, b_name, locfilter, strategy):
    """
    :param dataset:
    :param a_genes:
    :param b_genes:
    :param a_name:
    :param b_name:
    :param locfilter:
    :param strategy:
    :return:
    """
    agname = '{}_name'.format(a_name)
    a_select = a_genes['chrom'].str.match(locfilter, as_indexer=True)
    a_subset = a_genes.loc[a_select, agname].unique()
    bgname = '{}_name'.format(b_name)
    b_select = b_genes['chrom'].str.match(locfilter, as_indexer=True)
    b_subset = b_genes.loc[b_select, bgname].unique()

    a_idx = dataset['gene_name'].isin(a_subset)
    b_idx = dataset['gene_name'].isin(b_subset)

    a_data = dataset.loc[a_idx, ['clade_id', 'og_id', 'gene_name']]
    a_data.columns = ['clade_id', 'og_id', agname]
    b_data = dataset.loc[b_idx, ['clade_id', 'og_id', 'gene_name']]
    b_data.columns = ['clade_id', 'og_id', bgname]

    merged = a_data.merge(b_data, on=['clade_id', 'og_id'], how='outer')
    merged = merged.dropna(axis=0, how='any', subset=[agname, bgname], inplace=False)

    clade_counts = col.Counter(merged['clade_id'])
    if strategy == 'bottom-up':
        counts = reversed(clade_counts.most_common())
    else:
        counts = clade_counts.most_common()
    remove_ogids = set()
    select_ogids = set()
    select_genes = set()
    for cid, _ in counts:
        subs = merged.loc[merged['clade_id'] == cid, :]
        if subs.empty:
            continue
        ogid_counts = subs.groupby('og_id').count()
        match_idx = ogid_counts[agname] == ogid_counts[bgname]
        rem_ogids = ogid_counts.loc[~match_idx, :].index.tolist()
        remove_ogids.union(rem_ogids)
        check_ogids = ogid_counts.loc[match_idx, :].index.unique()
        for ogid in check_ogids:
            sub_sub = subs.loc[subs['og_id'] == ogid, :]
            if sub_sub[agname].unique().size != sub_sub[agname].size:
                remove_ogids.add(ogid)
                continue
            if sub_sub[bgname].unique().size != sub_sub[bgname].size:
                remove_ogids.add(ogid)
                continue
            if (sub_sub[agname].isin(select_genes)).any():
                remove_ogids.add(ogid)
                continue
            if (sub_sub[bgname].isin(select_genes)).any():
                remove_ogids.add(ogid)
                continue
            select_ogids.add(ogid)
            select_genes = select_genes.union(sub_sub[agname].tolist())
            select_genes = select_genes.union(sub_sub[bgname].tolist())
        merged = merged.loc[~merged['og_id'].isin(remove_ogids), :]
        if merged.empty:
            break
    a_select_ogids = a_data['og_id'].isin(select_ogids)
    a_select_genes = a_data[agname].isin(select_genes)
    a_select_both = np.logical_and(a_select_ogids, a_select_genes)
    a_data_select = a_data.loc[a_select_both, :]
    a_data_select.columns = ['clade_id', 'og_id', agname]
    a_data_select = a_data_select.merge(a_genes, on=agname, how='inner')
    a_data_select.drop('chrom', axis=1, inplace=True)
    asize = a_data_select[agname].size
    auniq = a_data_select[agname].unique().size
    assert asize == auniq,\
        'A genes selected multiple times: {} vs {}'.format(auniq, asize)
    b_select_ogids = b_data['og_id'].isin(select_ogids)
    b_select_genes = b_data[bgname].isin(select_genes)
    b_select_both = np.logical_and(b_select_ogids, b_select_genes)
    b_data_select = b_data.loc[b_select_both, :]
    b_data_select.columns = ['clade_id', 'og_id', bgname]
    b_data_select = b_data_select.merge(b_genes, on=bgname, how='inner')
    b_data_select.drop('chrom', axis=1, inplace=True)
    bsize = b_data_select[bgname].size
    buniq = b_data_select[bgname].unique().size
    assert bsize == buniq,\
        'B genes selected multiple times: {} vs {}'.format(buniq, bsize)
    assert a_data_select.shape == b_data_select.shape, 'Cannot merge: {} vs {}'.format(a_data_select.shape, b_data_select.shape)
    shared_select = a_data_select.merge(b_data_select, on=['og_id', 'clade_id'], copy=True)
    symbol_cols = [c for c in shared_select.columns if c.endswith('symbol')]
    shared_select = shared_select.sort_values(symbol_cols, axis=0, na_position='first')
    # extend OG ids with species
    # why this?
    # OrthoDB OG IDs are not stable between species pairs,
    # apparently they depend on the relation between the two.
    # Hence, merging on OG IDs leads to a ridiculously small
    # dataset for all 6 species. Augment OG IDs with species
    # information and merge on gene names/symbols only
    # (see: select_ortholog_groups)
    ogids = shared_select['og_id'].tolist()
    ogids = [item + '_' + a_name + '_' + b_name for item in ogids]
    shared_select['og_id_{}_{}'.format(a_name, b_name)] = ogids
    cladeids = shared_select['clade_id'].tolist()
    cladeids = [item + '_' + a_name + '_' + b_name for item in cladeids]
    shared_select['clade_id_{}_{}'.format(a_name, b_name)] = cladeids
    shared_select.drop(['og_id', 'clade_id'], axis=1, inplace=True)
    return shared_select


def select_ortholog_groups(fpath, group_root):
    """
    :param fpath:
    :param group_root:
    :return:
    """
    if not group_root.startswith('/'):
        group_root = '/' + group_root
    groups = None
    with pd.HDFStore(fpath, 'r') as hdf:
        load_keys = [k for k in hdf.keys() if k.startswith(group_root)]
        for k in load_keys:
            if groups is None:
                groups = hdf[k]
            else:
                shared_cols = set(groups.columns).intersection(set(hdf[k].columns))
                groups = groups.merge(hdf[k], on=list(shared_cols),
                                      how='outer', suffixes=('', ''))
    if groups is None or groups.empty:
        return None
    groups.dropna(axis=0, how='any', inplace=True)
    groups.reset_index(drop=True, inplace=True)
    return groups


def dump_orthodb_data(dataset, species_a, species_b, group_root, outpath, mode):
    """
    :param dataset:
    :param species_a:
    :param species_b:
    :param group_root:
    :param outpath:
    :param mode:
    :return:
    """
    if species_a and species_b:
        group = os.path.join(group_root, species_a, species_b)
    else:
        group = group_root
    if dataset is None or dataset.empty:
        return
    with pd.HDFStore(outpath, mode, complib='blosc', complevel=9) as hdf:
        hdf.put(group, dataset, format='table')
        hdf.flush()
    return


def main():
    """
    :return:
    """
    args = parse_command_line()
    spec_lut = read_species_table(args.speciestable, args.annotation)
    # create or load cached data
    if args.cache != 'no_cache':
        cache_dir = os.path.dirname(args.inputfiles[0])
        setattr(args, 'cache', os.path.join(cache_dir, args.cache))
    if args.cache == 'no_cache':
        merged = raw_process_input(args, spec_lut)
    elif os.path.isfile(args.cache):
        with pd.HDFStore(args.cache, 'r') as hdf:
            merged = hdf['cache']
    else:
        merged = raw_process_input(args, spec_lut)
        with pd.HDFStore(args.cache, 'w', complib='blosc', complevel=9) as hdf:
            hdf.put('cache', merged, format='table')
    with pd.HDFStore(args.outputfile, 'w', complib='blosc', complevel=9) as hdf:
        hdf.put('/raw', merged, format='table')
    all_species = {v['name'] for _, v in spec_lut.items() if 'name' in v}
    assert len(all_species) > 0, 'No species selected from annotation table'
    #for primary in ['human', 'mouse']:
    # assert that human and mouse are always first in pair
    primary_species = ['human', 'mouse'] + sorted(set(all_species) - {'human', 'mouse'})
    for primary in primary_species:
        prime_genes = read_gene_model(spec_lut[primary]['genes'], primary)
        for secondary in sorted(all_species):
            if primary == secondary:
                continue
            sec_genes = read_gene_model(spec_lut[secondary]['genes'], secondary)
            if sec_genes is None:
                continue
            pairs = select_ortholog_pairs(merged.copy(), prime_genes.copy(), sec_genes.copy(),
                                          primary, secondary, 'chr[0-9A-Z]+$', args.strategy)
            dump_orthodb_data(pairs, primary, secondary, 'augo/pairs', args.outputfile, 'a')
            pairs = select_ortholog_pairs(merged.copy(), prime_genes.copy(), sec_genes.copy(),
                                          primary, secondary, 'chr[0-9A-F]+$', args.strategy)
            dump_orthodb_data(pairs, primary, secondary, 'auto/pairs', args.outputfile, 'a')

    groups = select_ortholog_groups(args.outputfile, 'auto/pairs')
    dump_orthodb_data(groups, '', '', 'auto/groups', args.outputfile, 'a')

    groups = select_ortholog_groups(args.outputfile, 'augo/pairs')
    dump_orthodb_data(groups, '', '', 'augo/groups', args.outputfile, 'a')
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
