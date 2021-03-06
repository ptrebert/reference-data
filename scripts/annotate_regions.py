#!/usr/bin/env python3
# coding=utf-8

import sys as sys
import os as os
import traceback as trb
import argparse as argp
import csv as csv
import multiprocessing as mp
import io as io
import twobitreader as tbr


__VAVOURI__ = 'DOI:10.1186/gb-2012-13-11-r110'
__VAVOURI_WINDOW__ = 500
__VAVOURI_STEPSIZE__ = 50  # not specified in paper
__VAVOURI_OFFSET_UP__ = 300
__VAVOURI_OFFSET_DOWN__ = 200


def get_gc_content(seq):
    total = seq.count('c') + seq.count('g')
    return round((total / len(seq)) * 100., 3)


def get_obs_exp_ratio(seq):
    """
    Obs/Exp CpG = Number of CpG * N / (Number of C * Number of G)
    :param seq:
    :return:
    """
    cpg = seq.count('cg')
    # fudge factors to avoid ZeroDivision
    c = max(1, seq.count('c'))
    g = max(1, seq.count('g'))
    return round((cpg * len(seq)) / (c * g), 3)


def get_cpg_content(seq):
    """
    :param seq:
    :return:
    """
    cpg = seq.count('cg')
    return round((cpg / len(seq) * 100) * 2, 3)


def vavouri_promoter_class(region):
    """
    :param region:
    :return:
    """
    if region['strand'] == '-':
        start = int(region['start']) + __VAVOURI_OFFSET_DOWN__
        end = int(region['end']) - __VAVOURI_OFFSET_UP__
    else:
        start = int(region['start']) + __VAVOURI_OFFSET_UP__
        end = int(region['end']) - __VAVOURI_OFFSET_DOWN__
    assert start < end, 'Invalid region: {} - {} - {}'.format(start, end, region)
    scan_length = end - start - __VAVOURI_WINDOW__
    reg_length = int(region['end']) - int(region['start'])
    seq = region['seq']
    assert len(seq) == reg_length, 'Region - sequence mismatch: {} - {}'.format(len(seq), reg_length)
    lcp = 1
    hcp = False
    # +1 here to make it right-inclusive
    for idx in range(0, scan_length + 1, __VAVOURI_STEPSIZE__):
        subseq = seq[idx:idx+__VAVOURI_WINDOW__].lower()
        gc = get_gc_content(subseq)
        oe = get_obs_exp_ratio(subseq)
        if oe > 0.75 and gc > 55.:
            region['promoter_class'] = 'HCP'
            assert start + idx + __VAVOURI_WINDOW__ <= end, 'Stepped out of region: {}'.format(start + idx + __VAVOURI_WINDOW__)
            region['hcp_win_start'] = start + idx
            region['hcp_win_end'] = start + idx + __VAVOURI_WINDOW__
            hcp = True
            break  # single window enough
        elif oe <= 0.48:
            lcp &= 1
        else:
            lcp &= 0
    if not hcp:
        region['hcp_win_start'] = 'n/a'
        region['hcp_win_end'] = 'n/a'
        if lcp > 0:
            region['promoter_class'] = 'LCP'
        else:
            region['promoter_class'] = 'ICP'
    region['promoter_schema'] = 'vavouri'
    return region


def generic_gc_features(region):
    """
    :param region:
    :return:
    """
    s = region['seq'].lower()
    gc = get_gc_content(s)
    oe = get_obs_exp_ratio(s)
    cpg = get_cpg_content(s)
    region['gc_pct'] = gc
    region['oe_ratio'] = oe
    region['cpg_pct'] = cpg
    return region


def add_region_sequence(allregions, genomefile):
    """
    :param allregions:
    :param genomefile:
    :return:
    """
    genome = tbr.TwoBitFile(genomefile)
    current_chrom = ''
    chrom_seq = ''
    for reg in allregions:
        if reg['#chrom'] != current_chrom:
            current_chrom = reg['#chrom']
            chrom_seq = genome[current_chrom]
        reg['seq'] = chrom_seq[int(reg['start']):int(reg['end'])]
    return allregions


def read_region_file(fpath):
    """
    :param fpath:
    :return:
    """
    regions = []
    with open(fpath, 'r', newline='') as inf:
        rows = csv.DictReader(inf, delimiter='\t')
        fields = rows.fieldnames
        for r in rows:
            regions.append(r)
    return regions, fields


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, dest='input')
    parser.add_argument('--genome', '-g', type=str, dest='genome')
    parser.add_argument('--output', '-o', type=str, dest='output')
    parser.add_argument('--feature', '-f', type=str, nargs='*',
                        choices=['gc', 'vavouri'], dest='feature')
    parser.add_argument('--dump-fasta', action='store_true', default=False, dest='dumpfasta')
    parser.add_argument('--worker', '-w', type=int, default=1, dest='worker')
    args = parser.parse_args()
    return args


def make_fasta_header(region, suffix):
    """
    TODO: this function needs to be adapted to be more
    generic - right now, it makes too many assumptions
    about available header fields
    :param region:
    :return:
    """
    header = '>'
    loc = region['#chrom'] + ':' + str(region['start']) + '-' + str(region['end']) + '|'
    header += loc
    try:
        header += '|'.join([region['gene_name'], region['gene_id'], region['strand'], region['promoter_class']])
    except KeyError:
        header += region['name']
    if suffix:
        header += '|' + suffix
    region['header'] = header
    return region


def dump_regions(regions, fieldnames, output, fasta):
    """
    TODO: this function needs to be simplified in order to not
    contain code specific for Sarvesh's project
    :param regions:
    :param fieldnames:
    :param output:
    :return:
    """
    with open(output, 'w', newline='') as outf:
        writer = csv.DictWriter(outf, delimiter='\t', fieldnames=fieldnames,
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(regions)
    if fasta:
        # for Sarvesh
        if '_cgi' in os.path.basename(output):
            suffix = 'CGI'
        elif '_noncgi' in os.path.basename(output):
            suffix = 'nonCGI'
        elif '_distal' in os.path.basename(output):
            suffix = 'DIST'
        elif '_proximal' in os.path.basename(output):
            suffix = 'PROX'
        else:
            suffix = ''
        regions = [make_fasta_header(r, suffix) for r in regions]
        fabuffer = io.StringIO()
        for reg in regions:
            fabuffer.write(reg['header'] + '\n')
            fabuffer.write(reg['seq'] + '\n\n')
        base = output.rsplit('.', 1)[0]
        fastafile = base + '.fa'
        with open(fastafile, 'w') as outf:
            _ = outf.write(fabuffer.getvalue())
    return


def run_annotate():
    """
    :return:
    """
    exc = 0
    try:
        args = parse_command_line()
        regions, fields = read_region_file(args.input)
        regions = add_region_sequence(regions, args.genome)
        with mp.Pool(args.worker) as pool:
            if 'gc' in args.feature:
                fields.extend(['gc_pct', 'oe_ratio', 'cpg_pct'])
                regions = pool.map(generic_gc_features, regions)
            if 'vavouri' in args.feature:
                fields.extend(['promoter_class', 'promoter_schema', 'hcp_win_start', 'hcp_win_end'])
                regions = pool.map(vavouri_promoter_class, regions)
        dump_regions(regions, fields, args.output, args.dumpfasta)
    except Exception as err:
        trb.print_exc()
        sys.stderr.write('\nError: {}\n'.format(err))
        if isinstance(err.args[0], int):
            exc = err.args[0]
        else:
            exc = 1
    finally:
        return exc


if __name__ == '__main__':
    stat = run_annotate()
    sys.exit(stat)

