# coding=utf-8

import os as os
import sys as sys
import collections as col
import csv as csv
import operator as op
import fnmatch as fnm
import io as io

from pipelines.auxmod.auxiliary import read_chromsizes, open_comp, check_bounds


def derive_encode_enh_name(names):
    """
    :param names:
    :return:
    """
    all_names = names.split('@')
    all_names = col.Counter([n.rsplit('-', 1)[0] for n in all_names])
    types = {'Proximal-Prediction': 'PP', 'Distal-Prediction': 'DP'}
    return types[all_names.most_common()[0][0]]


def process_merged_encode_enhancer(inputfile, outputfile, boundcheck, nprefix, etype):
    """
    :param inputfile:
    :param outputfile:
    :param boundcheck:
    :param nprefix:
    :param etype:
    :return:
    """
    bounds = read_chromsizes(boundcheck)
    opn, mode, conv = open_comp(inputfile, True)
    regions = []
    enhname = derive_encode_enh_name
    bchk = check_bounds
    with opn(inputfile, 'r') as infile:
        for line in infile:
            line = conv(line)
            if not line:
                continue
            c, s, e, n = line.split()
            stat = bchk(c, s, e, bounds, inputfile)
            if not stat:
                # chromosome not in check file, otherwise raises
                continue
            enhtype = enhname(n)
            if enhtype != etype:
                continue
            regions.append((c, s, e, enhtype))
    assert regions, 'No regions read from file {}'.format(inputfile)
    regions = sorted(regions, key=lambda x: (x[0], int(x[1])))
    outbuffer = io.StringIO()
    for idx, reg in enumerate(regions, start=1):
        outbuffer.write('\t'.join(reg[:3]) + '\t{}{}_{}\n'.format(nprefix, reg[3], idx))
    opn, mode, conv = open_comp(outputfile, False)
    with opn(outputfile, mode) as outf:
        _ = outf.write(conv(outbuffer.getvalue()))
    return outputfile


def process_vista_enhancer(inputfile, outputfile, boundcheck, nprefix, species):
    """
    :param inputfile:
    :param outputfile:
    :param boundcheck:
    :param nprefix:
    :param species:
    :return:
    """
    bounds = read_chromsizes(boundcheck)
    opn, mode, conv = open_comp(inputfile, True)
    regions = []
    bchk = check_bounds
    with opn(inputfile, mode) as infile:
        for line in infile:
            line = conv(line)
            if not line or not line.startswith('>'):
                continue
            header = line.strip('>').split('|')
            if header[0] != species:
                continue
            c, coord = header[1].strip().split(':')
            s, e = coord.split('-')
            stat = bchk(c, s, e, bounds, inputfile)
            if not stat:
                # chromosome not in check file, otherwise raises
                continue
            regname = nprefix + '_' + header[3].strip()[:3] + '_' + header[2].strip('elmnt ')
            regions.append((c, s, e, regname))
    assert regions, 'No regions read from file {}'.format(inputfile)
    regions = sorted(regions, key=lambda x: (x[0], int(x[1]), int(x[2])))
    outbuffer = io.StringIO()
    for reg in regions:
        outbuffer.write('\t'.join(reg) + '\n')
    opn, mode, conv = open_comp(outputfile, False)
    with opn(outputfile, mode) as outf:
        _ = outf.write(conv(outbuffer.getvalue()))
    return outputfile


def merge_genehancer_annotation(inputfiles, outputfile):
    """
    :param inputfiles:
    :param outputfile:
    :return:
    """
    enh_collect = dict()
    with open(inputfiles[0], 'r', encoding='utf-8', newline='\n') as enh_file:
        assert '\r\n' not in enh_file.readline(), 'Linefeed detected in {}'.format(inputfiles[0])
        # _ = enh_file.readline()
        header = ['chrom', 'start', 'end', 'cluster_id',  'GHid', 'enhancer_score', 'is_elite']
        rows = csv.DictReader(enh_file, delimiter='\t', fieldnames=header)
        for row in rows:
            row['chrom'] = 'chr' + row['chrom']
            assert row['cluster_id'] not in enh_collect, 'Cluster duplicate: {}'.format(row)
            enh_collect[row['cluster_id']] = row
    assoc_collect = col.defaultdict(list)
    with open(inputfiles[1], 'r', encoding='utf-8', newline='\n') as assoc:
        assert '\r\n' not in assoc.readline(), 'Linefeed detected in {}'.format(inputfiles[1])
        # _ = assoc.readline()

        header = ['cluster_id', 'gene_id', 'enh_gene_dist', 'assoc_score', 'is_elite']
        rows = csv.DictReader(assoc, delimiter='\t', fieldnames=header)
        for row in rows:
            if row['gene_id'].startswith('ENSG'):
                # could be Ensembl ID
                try:
                    check = row['gene_id']
                    _ = int(check.replace('ENSG', ''))
                    row['name'] = row['gene_id']
                    row['symbol'] = 'n/a'
                except ValueError:
                    row['symbol'] = row['gene_id']
                    row['name'] = 'n/a'
            else:
                row['symbol'] = row['gene_id']
                row['name'] = 'n/a'
            assoc_collect[row['cluster_id']].append(row)
    final = []
    for cluster, enh in enh_collect.items():
        try:
            genes = assoc_collect[cluster]
            vals = [(x['name'], x['symbol'], x['assoc_score'], x['enh_gene_dist'])
                    for x in sorted(genes, key=lambda x: float(x['assoc_score']), reverse=True)]
            names = ','.join([x[0] for x in vals])
            symbols = ','.join([x[1] for x in vals])
            scores = ','.join([x[2] for x in vals])
            dists = ','.join([x[3] for x in vals])
            enh['name'] = names
            enh['symbol'] = symbols
            enh['assoc_score'] = scores
            enh['enh_gene_dist'] = dists
        except KeyError:
            enh['name'] = 'n/a'
            enh['symbol'] = 'n/a'
            enh['assoc_score'] = '0'
            enh['enh_gene_dist'] = '0'
        finally:
            final.append(enh)
    final = sorted(final, key=lambda x: (x['chrom'], int(x['start']), int(x['end'])))
    header = ['chrom', 'start', 'end', 'GHid', 'enhancer_score', 'is_elite', 'cluster_id',
              'name', 'symbol', 'assoc_score', 'enh_gene_dist']
    with open(outputfile, 'w', encoding='utf-8', newline='\n') as out:
        _ = out.write('#')
        writer = csv.DictWriter(out, fieldnames=header, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(final)
    return out


def process_hsa_mapped(inputfile, outputfile, raw_enhancers):
    """
    :param inputfile:
    :param outputfile:
    :return:
    """
    enh_lengths = dict()
    with open(raw_enhancers, 'r') as enh:
        _ = enh.readline()
        for line in enh:
            parts = line.split()
            enh_lengths[parts[4]] = int(parts[2]) - int(parts[1])

    header = ['chrom', 'start', 'end', 'GHid', 'enhancer_score', 'is_elite', 'cluster_id',
              'name', 'symbol', 'assoc_score', 'enh_gene_dist']
    enh_collect = col.defaultdict(list)
    chrom_collect = col.defaultdict(set)
    with open(inputfile, 'r') as mapped:
        rows = csv.DictReader(mapped, delimiter='\t', fieldnames=header)
        for row in rows:
            enh_collect[row['GHid']].append(row)
            chrom_collect[row['GHid']].add(row['chrom'])
    enh_out = []
    for ghid, enhancer in enh_collect.items():
        chroms = chrom_collect[ghid]
        if len(chroms) > 2:
            # 5 enhancer were scattered during the mapping; discard
            continue
        elif len(chroms) == 2:
            # select majority fragments; needed for 18 enhancers
            enhancer = select_map_majority(enhancer)
        else:
            pass
        ghid = enhancer[0]['GHid']
        orig_len = enh_lengths[ghid]
        if len(enhancer) == 1:
            new_len = int(enhancer[0]['end']) - int(enhancer[0]['start'])
            if new_len <= orig_len * 1.1:
                enh_out.append(enhancer[0])
            else:
                raise AssertionError('New enhancer too long: {} {}'.format(orig_len, enhancer[0]))
        else:
            # this is the case for ~1500 enhancers
            new_enh = join_enh_fragments(enhancer, orig_len)
            enh_out.append(new_enh)

    enh_out = sorted(enh_out, key=lambda x: (x['chrom'], int(x['start']), int(x['end'])))

    out_header = header[:-1]
    with open(outputfile, 'w') as out:
        _ = out.write('#')
        writer = csv.DictWriter(out, fieldnames=out_header, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(enh_out)
    return outputfile


def join_enh_fragments(fragments, target_len):
    """
    :param fragments:
    :param target_len:
    :return:
    """
    enh = None
    l = 0
    for f in fragments:
        f['start'] = int(f['start'])
        f['end'] = int(f['end'])
        fl = f['end'] - f['start']
        f['length'] = fl
        f['mid'] = f['start'] + int(f['length'] / 2)
        if fl <= (target_len * 1.1):
            if fl > l:
                enh = f
                l = fl
    lo, hi = enh['mid'] - target_len, enh['mid'] + target_len
    s, e = enh['start'], enh['end']
    for f in fragments:
        if int(f['start']) < hi and int(f['end']) > lo:
            ns, ne = min(s, f['start']), max(e, f['end'])
            nl = ne - ns
            if nl > target_len * 1.1:
                break
            else:
                s, e = ns, ne
                mid = s + int(nl / 2)
                lo, hi = mid - target_len, mid + target_len
    assert (e - s) < target_len * 1.1, 'Joined fragments too long: {}\n{}\n{}'.format((e - s), enh, fragments)
    assert s < e, 'Invalid coordinates: {} - {}'.format(s, e)
    enh['start'] = str(s)
    enh['end'] = str(e)
    return enh


def select_map_majority(fragments):
    """
    :param fragments:
    :return:
    """
    c = col.Counter()
    for f in fragments:
        c[f['chrom']] += 1
    mc = c.most_common(1)
    fragments = [f for f in fragments if f['chrom'] == mc[0][0]]
    return fragments


def build_genehancer_map_params(inputfile, chainfiles, outdir, cmd, jobcall):
    """
    :param inputfile:
    :param chainfiles:
    :param outdir:
    :param cmd:
    :param jobcall:
    :return:
    """
    arglist = []
    for chf in chainfiles:
        fp, fn = os.path.split(chf)
        qry = fn.split('.')[0].split('_to_')[1]
        outname = '{}_enh_genehancer_raw.bed'.format(qry)
        outpath = os.path.join(outdir, outname)
        tmp = cmd.format(**{'chainfile': chf})
        arglist.append([inputfile, outpath, tmp, jobcall])
    if chainfiles and os.path.isfile(inputfile):
        assert arglist, 'No arguments created'
    return arglist


def build_genehancer_dset_params(inputfiles, outdir, seqfiles, genefiles, cmd, jobcall):
    """
    :param inputfiles:
    :param outdir:
    :param seqfiles:
    :param genefiles:
    :param cmd:
    :param jobcall:
    :return:
    """
    arglist = []
    for inpf in inputfiles:
        fp, fn = os.path.split(inpf)
        assm = fn.split('_')[0]
        seq = fnm.filter(seqfiles, '*/' + assm + '*.2bit')
        try:
            assert len(seq) == 1, 'Genomic sequence error: {} - {}'.format(assm, seq)
        except AssertionError:
            if assm in ['mm10', 'hg38', 'ihec38', 'GRCh37']:
                continue
            raise
        seqfile = seq[0]
        genes = fnm.filter(genefiles, '*' + assm + '*reg5p.h5')
        try:
            assert len(genes) == 1, 'Gene file error: {} - {}'.format(assm, genes)
        except AssertionError:
            if assm in ['mm10', 'hg38', 'ihec38', 'GRCh37']:
                continue
            raise
        genefile = genes[0]
        tmp = cmd.format(**{'sequence': seqfile, 'genes': genefile})
        outname = assm + '_gh_dataset.h5'
        outpath = os.path.join(outdir, outname)
        arglist.append([inpf, outpath, tmp, jobcall])
    if inputfiles:
        assert arglist, 'No arguments created'
    return arglist
