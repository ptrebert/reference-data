#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import pandas as pd
import numpy as np
import collections as col
import traceback as trb
import argparse as argp
import fnmatch as fnm


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
                        dest='strategy')
    args = parser.parse_args()
    return args


def tax_ids():
    """
    :return:
    """
    ncbi_ids = [9615, 10090, 9606, 9913, 9823, 9031]
    regular = ['dog', 'mouse', 'human', 'cow', 'pig', 'chicken']
    abbr = ['cfa', 'mmu', 'hsa', 'bta', 'ssc', 'gga']
    lookup = dict()
    for t, n, a in zip(ncbi_ids, regular, abbr):
        lookup[t] = {'name': n, 'code': a}
    return lookup


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


def merge_tables(genes, ogroups, mapping):
    """
    :param genes:
    :param ogroups:
    :param mapping:
    :return:
    """
    taxons = tax_ids()
    indexer = genes['tax_id'].isin(taxons.keys())
    species_genes = (genes.loc[indexer, :]).copy(deep=True)
    indexer = mapping['odb_gene_id'].isin(species_genes['odb_gene_id'])
    species_map = (mapping.loc[indexer, :]).copy(deep=True)
    indexer = ogroups['og_id'].isin(species_map['og_id'])
    species_og = (ogroups.loc[indexer, :]).copy(deep=True)

    merged = species_og.merge(species_map, on='og_id', how='outer')
    merged = merged.merge(species_genes, on='odb_gene_id', how='outer')
    merged['clade_id'] = merged['og_id'].str.slice(start=0, stop=9)
    return merged


def raw_process_input(args):
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
    merged = merge_tables(genes_tab, og_tab, og2genes_tab)
    return merged


def read_gene_model(basedir, species):
    """
    :param basedir:
    :param species:
    :return:
    """
    abbr = SPECIES_MAP[species]
    modelfiles = os.listdir(basedir)
    modelfile = fnm.filter(modelfiles, '{}_*.genes.bed.gz'.format(abbr))
    assert len(modelfile) == 1, 'No gene set identified: {} / {}'.format(species, abbr)
    fpath = os.path.join(basedir, modelfile[0])
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
        rem_ogids = ogid_counts.loc[~match_idx, :].index
        remove_ogids.add(rem_ogids)
        ogid_counts = ogid_counts.loc[, :]
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
            select_genes.add(sub_sub[agname])
            select_genes.add(sub_sub[bgname])
        merged = merged.loc[~merged['og_id'].isin(remove_ogids), :]

    # continue here: select og_ids from A and B, merge and return

    return 'foo'


def reduce_to_subset(df, subsetpath):
    """
    :param df:
    :param subsetpath:
    :return:
    """
    # what follows:
    # the OrthoDB flat files seem not to have a clear info about the "hierarchy"
    # among the clades, so this simple heuristic selects the maximal set of OG IDs
    # that do not contain gene duplicates

    # first 9 characters of OG ID are the clade ID
    clade_ids = set(joined_subset['og_id'].str.slice(start=0, stop=9).unique())
    # make a temp working copy of the dataset
    clade_data = joined_subset.copy()
    clade_counter = col.Counter()
    for cid in clade_ids:
        subs = clade_data.loc[clade_data['og_id'].str.startswith(cid), :]
        clade_counter[cid] = subs.shape[0]
    select_ogids = set()
    remove_ogids = set()
    selected_genes = set()
    for cid, _ in reversed(clade_counter.most_common()):
        subs = clade_data.loc[clade_data['og_id'].str.startswith(cid), :]
        if subs.empty:
            continue
        for ogid in subs['og_id'].unique():
            sub_spec = subs.loc[subs['og_id'] == ogid, 'species'].unique().tolist()
            if len(sub_spec) != 5:
                remove_ogids.add(ogid)
                continue
            sub_genes = subs.loc[subs['og_id'] == ogid, 'gene_name'].tolist()
            if len(set(sub_genes)) != len(sub_genes):
                # no clue what kind of ortholog group could contain the same
                # gene twice or more - what is that? just discard...
                remove_ogids.add(ogid)
                continue
            if any([g in selected_genes for g in sub_genes]):
                remove_ogids.add(ogid)
            else:
                remove_ogids.add(ogid)
                select_ogids.add(ogid)
                selected_genes = selected_genes.union(set(sub_genes))
        idx_remove = clade_data['og_id'].isin(remove_ogids)
        clade_data = clade_data.loc[~idx_remove, :]
        if clade_data.empty:
            # this should happen automatically...
            break
    uniq_idx = joined_subset['og_id'].isin(select_ogids)
    joined_subset = joined_subset.loc[uniq_idx, :]
    # record info about group sizes
    ogid_abund = col.Counter(joined_subset['og_id'])
    group_size = []
    group_bal = []
    for ogid, count in ogid_abund.most_common():
        group_size.append([ogid, count])
        # check if the group is at least balanced
        this_group = col.Counter(joined_subset.loc[joined_subset['og_id'] == ogid, 'species'])
        median = np.median([n for i, n in this_group.most_common()])
        balanced = True
        for i, n in this_group.most_common():
            if (median - 1) <= n <= (median + 1):
                pass
            else:
                balanced = False
                break
        if balanced:
            group_bal.append([ogid, 1])
        else:
            group_bal.append([ogid, 0])
    group_size = pd.DataFrame(group_size, columns=['og_id', 'group_size'])
    group_bal = pd.DataFrame(group_bal, columns=['og_id', 'group_balanced'])
    joined_subset = joined_subset.merge(group_size, on='og_id', how='outer')
    assert joined_subset['group_size'].notnull().all(), 'Failed'
    joined_subset = joined_subset.merge(group_bal, on='og_id', how='outer')
    assert joined_subset['group_balanced'].notnull().all(), 'Failed'
    num_genes = joined_subset['gene_name'].size
    uniq_genes = joined_subset['gene_name'].unique().size
    assert num_genes == uniq_genes, 'Mismatch: {} vs {}'.format(num_genes, uniq_genes)
    joined_subset = joined_subset.astype({'og_id': str, 'og_name': str, 'odb_gene_id': str,
                                          'tax_id': np.int32, 'gene_name': str, 'chrom': str,
                                          'start': np.int32, 'end': np.int32, 'strand': np.int8,
                                          'symbol': str, 'species': str, 'group_size': np.int32,
                                          'group_balanced': np.int8}, copy=True)
    return joined_subset


def main():
    """
    :return:
    """
    args = parse_command_line()
    # create or load cached data
    if args.cache != 'no_cache':
        cache_dir = os.path.dirname(args.inputfiles[0])
        setattr(args, 'cache', os.path.join(cache_dir, args.cache))
    if args.cache == 'no_cache':
        merged = raw_process_input(args)
    elif os.path.isfile(args.cache):
        with pd.HDFStore(args.cache, 'r') as hdf:
            merged = hdf['cache']
    else:
        merged = raw_process_input(args)
        with pd.HDFStore(args.cache, 'w', complib='blosc', complevel=9) as hdf:
            hdf.put('cache', merged, format='table')

    for primary in ['human', 'mouse']:
        prime_genes = read_gene_model(args.annotation, primary)
        for secondary in SPECIES_MAP.keys():
            if primary == secondary:
                continue
            sec_genes = read_gene_model(args.annotation, secondary)
            pairs = select_ortholog_pairs(merged.copy(), prime_genes.copy(), sec_genes.copy(),
                                          primary, secondary, 'chr[0-9XYZW]+$')
            raise
    #
    # indexer = merged['gene_name'].notnull()
    # merged = merged.loc[indexer, :]
    # indexer = merged['og_id'].notnull()
    # raw_merged = merged.loc[indexer, :].copy()
    # with pd.HDFStore(args.outputfile, 'w', complib='blosc', complevel=9) as hdf:
    #     hdf.put('/raw', raw_merged, format='table')
    # species = tax_ids()
    # shared_entries = set()
    # for k, vals in species.items():
    #     this_species = raw_merged.loc[merged['tax_id'] == k, 'og_id']
    #     if not shared_entries:
    #         shared_entries = set(this_species.tolist())
    #     else:
    #         shared_entries = shared_entries.intersection(set(this_species.tolist()))
    # indexer = raw_merged['og_id'].isin(shared_entries)
    # shared_merged = raw_merged.loc[indexer, :].copy()
    # with pd.HDFStore(args.outputfile, 'a', complib='blosc', complevel=9) as hdf:
    #     hdf.put('/shared/raw', shared_merged, format='table')
    # subset = reduce_to_subset(shared_merged, args.subsets)
    # assert not subset.empty, 'Created empty subset of data'
    # with pd.HDFStore(args.outputfile, 'a', complib='blosc', complevel=9) as hdf:
    #     hdf.put('/shared/subset', subset, format='table')
    # md = collect_subset_metadata(subset)
    # with pd.HDFStore(args.outputfile, 'a', complib='blosc', complevel=9) as hdf:
    #     hdf.put('/metadata', md, format='table')
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
