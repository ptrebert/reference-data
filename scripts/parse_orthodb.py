#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import pandas as pd
import numpy as np
import collections as col
import traceback as trb
import argparse as argp


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, nargs='+', dest='inputfiles')
    parser.add_argument('--cache-file', '-cf', type=str, dest='cache', default='no_cache')
    parser.add_argument('--subsets', '-sub', type=str, dest='subsets')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')
    args = parser.parse_args()
    return args


def tax_ids():
    """
    :return:
    """
    ncbi_ids = [9615, 10090, 9606, 9913, 9823]
    regular = ['dog', 'mouse', 'human', 'cow', 'pig']
    abbr = ['cfa', 'mmu', 'hsa', 'bta', 'ssc']
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
    return merged


def raw_process_input(args):
    """
    :param args:
    :return:
    """

    gene_tab, og_tab, map_tab = None, None, None
    for inpf in args.inputfiles:
        assert 'odb9' in inpf, 'Unexpected file: {}'.format(inpf)
        ttype = os.path.basename(inpf).split('.')[0].split('_')[1]
        if ttype == 'genes':
            assert gene_tab is None, 'Duplicate gene: {}'.format(inpf)
            gene_tab = pd.read_csv(inpf, delimiter='\t', header=None,
                                   usecols=['odb_gene_id', 'tax_id', 'gene_name'],
                                   names=odb_header(ttype), skip_blank_lines=True,
                                   dtype={'tax_id': np.int32, 'gene_name': str, 'odb_gene_id': str})
        elif ttype == 'OGs':
            assert og_tab is None, 'Duplicate OG table: {}'.format(inpf)
            og_tab = pd.read_csv(inpf, delimiter='\t', header=None,
                                 names=odb_header(ttype), skip_blank_lines=True,
                                 usecols=['og_id', 'og_name'], dtype={'og_id': str, 'og_name': str})
        elif ttype == 'OG2genes':
            assert map_tab is None, 'Duplicate map table: {}'.format(inpf)
            map_tab = pd.read_csv(inpf, delimiter='\t', header=None,
                                  names=odb_header(ttype), skip_blank_lines=True,
                                  dtype={'og_id': str, 'odb_gene_id': str})
        else:
            raise ValueError('Unexpected table type: {}'.format(inpf))
    merged = merge_tables(gene_tab, og_tab, map_tab)
    return merged


def read_subset(fpath):
    """
    :param fpath:
    :return:
    """
    def conv_strand(token):
        replace = {'+': 1, '-': -1}
        return replace[token]

    subset = pd.read_csv(fpath, delimiter='\t', header=0, skip_blank_lines=True,
                         names=['chrom', 'start', 'end', 'gene_name', 'score', 'strand', 'symbol'],
                         usecols=['chrom', 'start', 'end', 'gene_name', 'strand', 'symbol'],
                         converters={'strand': conv_strand},
                         dtype={'chrom': str, 'start': np.int32, 'end': np.int32,
                                'strand': np.int8, 'symbol': str, 'name': str},
                         na_values={'symbol': ['n/a']})
    return subset


def reduce_to_subset(df, subsetpath):
    """
    :param df:
    :param subsetpath:
    :return:
    """
    lookup = tax_ids()
    mapfiles = os.listdir(subsetpath)
    mapfiles = [os.path.join(subsetpath, fp) for fp in mapfiles if fp.endswith('.genes.bed.gz')]
    subset_shared = set()
    all_subsets = None
    for mf in mapfiles:
        species = os.path.basename(mf).split('_')[0]
        sub_genes = read_subset(mf)
        select_id = 0
        for taxid, vals in lookup.items():
            if species == vals['code']:
                select_id = taxid
                break
        assert select_id != 0, 'No taxon ID selected: {}'.format(species)
        spec_common = lookup[select_id]['name']
        sub_genes['species'] = spec_common
        if all_subsets is None:
            all_subsets = sub_genes.copy()
        else:
            all_subsets = pd.concat([all_subsets, sub_genes], axis=0,
                                    ignore_index=True, join='outer')
        subset = df.loc[df['gene_name'].isin(sub_genes['gene_name']), :]
        assert not subset.empty, 'No genes selected for species: {}'.format(species)
        if not subset_shared:
            subset_shared = set(subset['og_id'].tolist())
        else:
            subset_shared = subset_shared.intersection(set(subset['og_id'].tolist()))
    joined_subset = df.loc[df['og_id'].isin(subset_shared), :].copy()
    joined_subset = joined_subset.merge(all_subsets, on='gene_name', how='outer')
    joined_subset = joined_subset.dropna(axis=0, how='any', inplace=False,
                                         subset=['gene_name', 'og_id', 'species'])
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


def collect_subset_metadata(df):
    """
    :param df:
    :return:
    """
    md = [['subset', 'protein_coding']]
    if df['chrom'].isin(['chrX', 'chrY']).any():
        md.append(['chromosomes', 'whole_genome'])
    else:
        md.append(['chromosomes', 'autosomes'])
    species = df['species'].unique()
    for s in species:
        spec_sub = df.loc[df['species'] == s, 'gene_name']
        spec_sub = spec_sub.unique()
        md.append(['{}_genes'.format(s), str(spec_sub.size)])
    mdf = pd.DataFrame(md, columns=['key', 'value'])
    return mdf


def main():
    """
    :return:
    """
    args = parse_command_line()
    if args.cache != 'no_cache':
        cache_dir = os.path.dirname(args.inputfiles[0])
        setattr(args, 'cache', os.path.join(cache_dir, args.cache))
    if args.cache == 'no_cache':
        merged = raw_process_input(args)
    elif os.path.isfile(args.cache):
        with pd.HDFStore(args.cache, 'r') as hdf:
            merged = hdf['raw']
    else:
        merged = raw_process_input(args)
        with pd.HDFStore(args.cache, 'w', complib='blosc', complevel=9) as hdf:
            hdf.put('raw', merged, format='table')
    indexer = merged['gene_name'].notnull()
    merged = merged.loc[indexer, :]
    indexer = merged['og_id'].notnull()
    raw_merged = merged.loc[indexer, :].copy()
    with pd.HDFStore(args.outputfile, 'w', complib='blosc', complevel=9) as hdf:
        hdf.put('/raw', raw_merged, format='table')
    species = tax_ids()
    shared_entries = set()
    for k, vals in species.items():
        this_species = raw_merged.loc[merged['tax_id'] == k, 'og_id']
        if not shared_entries:
            shared_entries = set(this_species.tolist())
        else:
            shared_entries = shared_entries.intersection(set(this_species.tolist()))
    indexer = raw_merged['og_id'].isin(shared_entries)
    shared_merged = raw_merged.loc[indexer, :].copy()
    with pd.HDFStore(args.outputfile, 'a', complib='blosc', complevel=9) as hdf:
        hdf.put('/shared/raw', shared_merged, format='table')
    subset = reduce_to_subset(shared_merged, args.subsets)
    assert not subset.empty, 'Created empty subset of data'
    with pd.HDFStore(args.outputfile, 'a', complib='blosc', complevel=9) as hdf:
        hdf.put('/shared/subset', subset, format='table')
    md = collect_subset_metadata(subset)
    with pd.HDFStore(args.outputfile, 'a', complib='blosc', complevel=9) as hdf:
        hdf.put('/metadata', md, format='table')
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
