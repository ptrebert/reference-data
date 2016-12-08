#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import argparse as argp
import traceback as trb
import re as re
import gzip as gz
import bz2 as bz
import operator as op
import functools as fnt
import json as js
import numpy as np
import numpy.ma as msk
import io as io
import pickle as pck
import collections as col


Chain = col.namedtuple('Chain', ['tname', 'tsize', 'tstrand', 'tstart', 'tend',
                                 'qname', 'qsize', 'qstrand', 'qstart', 'qend',
                                 'chainid', 'score'])
Block = col.namedtuple('Block', ['tchrom', 'tstart', 'tend', 'tstrand',
                                 'qchrom', 'qstart', 'qend', 'qstrand',
                                 'name', 'score', 'size'])


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser(add_help=True)
    parser.add_argument('--task', '-tk', type=str, choices=['symmfilt'], dest='task', required=True)

    #parser.add_argument('--task', '-tk', type=str, choices=['chaintobed', 'qfilter', 'symmext',
    #                                                        'maptobedgraph', 'swap'], dest='task')

    comgroup = parser.add_argument_group('File arguments')
    comgroup.add_argument('--chain-file', '-chf', type=str, dest='chainfile',
                          help='A standard UCSC chain file')
    comgroup.add_argument('--map-file', '-mpf', type=str, dest='mapfile',
                          help='A one-to-one mapping file created from a chain file')
    comgroup.add_argument('--output-file', '-otf', type=str, dest='outputfile')

    comgroup = parser.add_argument_group('Task: Filter chains for symmetrical coverage [symmfilt]')
    comgroup.add_argument('--min-score', '-msc', type=int, dest='minscore', default=0)
    comgroup.add_argument('--min-size', '-msz', type=int, dest='minsize', default=1)
    comgroup.add_argument('--chrom', '-chr', type=str, dest='chrom')

    comgroup = parser.add_argument_group('Task: Dump chain to BED [chaintobed]')
    comgroup.add_argument('--no-buffer', '-nob', action='store_true', default=False, dest='nobuffer')

    comgroup = parser.add_argument_group('Task: Check for mergeable blocks in a mapping [symmext]')
    comgroup.add_argument('--query-sizes', '-qcs', type=str, dest='qsizes')
    comgroup.add_argument('--gap-size', '-gap', type=int, default=200, dest='gapsize')
    comgroup.add_argument('--cov-ratio', '-cov', type=float, default=0.7, dest='covratio')
    comgroup.add_argument('--min-cov', '-mcv', type=int, default=100, dest='mincov')
    comgroup.add_argument('--check-ext', '-cex', type=int, nargs='+', default=[100, 50, 25], dest='checkext')

    comgroup = parser.add_argument_group('Task: Dump mapping to bedGraph [maptobedgraph]')
    comgroup.add_argument('--query', '-qry', action='store_true', default=False,
                          help='Dump (unsorted) query coordinates; otherwise output is target coordinates')

    args = parser.parse_args()
    assert max(args.checkext) < args.gapsize,\
        'Maximal extension has to be shorter than gap size: {} vs {}'.format(max(args.checkext), args.gapsize)
    setattr(args, 'checkext', sorted(args.checkext, reverse=True))

    return args


def text_file_mode(fpath, read=True):
    """
    Naive determination of file type and appropriate selection
    of opening function and mode. Text-mode for compressed
    files requires Python 3.4+

    :param fpath:
    :return:
    """
    if read:
        assert os.path.isfile(fpath), 'Invalid path to file: {}'.format(fpath)
    ext = fpath.split('.')[-1].lower()
    if ext in ['gz', 'gzip']:
        f, m = gz.open, 'rt'  # Python 3.4+
    elif ext in ['bz', 'bz2', 'bzip', 'bzip2']:
        f, m = bz.open, 'rt'
    else:
        f, m = open, 'r'
    if not read:
        m = m.replace('r', 'w')
    return f, m


def identify_file_extension(fpath):
    """
    :param fpath:
    :return:
    """
    fp, fn = os.path.split(fpath)
    comp_ext = fn.rsplit('.', 1)[1]
    is_compressed = comp_ext in ['zip', 'gz', 'gzip', 'bz', 'bz2', 'bzip2']
    if is_compressed:
        ext = '.'.join(fn.rsplit('.', 2)[1:])
    else:
        ext = fn.rsplit('.', 1)[1]
    return ext


def read_chain_header(line):
    """
    Convenience wrapper to parse chain header information
    :param line:
    :return:
    """
    _, score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd, chainid = line.strip().split()
    return Chain(tName, int(tSize), tStrand, int(tStart), int(tEnd), qName, int(qSize), qStrand, int(qStart), int(qEnd), chainid, int(score))


def read_chain_block(line):
    size, dt, dq = line.strip().split()
    return int(size), int(dt), int(dq)


def write_block_to_buffer(block, iobuffer):
    """
    :param block:
    :param iobuffer:
    :return:
    """
    iobuffer.write('\t'.join([block.tchrom, str(block.tstart), str(block.tend), block.tstrand]))
    iobuffer.write('\t' + block.name + '\t')
    iobuffer.write('\t'.join([block.qchrom, str(block.qstart), str(block.qend), block.qstrand]) + '\n')
    return


def check_skip(selector, chrom):
    """
    :param selector:
    :param chrom:
    :return:
    """
    if selector is None:
        return False
    elif selector.match(chrom) is None:
        return True
    else:
        return False


def chromsize_from_chain(chainfile, chrom, target=True):
    """
    :param chainfile:
    :param chrom:
    :return:
    """
    read_head = read_chain_header
    opn, mode = text_file_mode(chainfile)
    chrom_size = 0
    with opn(chainfile, mode=mode, encoding='ascii') as chf:
        for line in chf:
            if line.strip() and line.startswith('chain'):
                chain = read_head(line)
                if chain.tname == chrom and target:
                    chrom_size = chain.tsize
                    break
                elif chain.qname == chrom:
                    chrom_size = chain.qsize
                    break
                else:
                    continue
    assert chrom_size > 0, 'No entry in chain file {} for chromosome: {}'.format(chainfile, chrom)
    return chrom_size


def get_chain_iterator(fobj, tselect=None, qselect=None, read_num=0, min_size=1, min_score=0):
    """
    Returns an iterator over chain files as used by UCSC liftOver tool.
    The design assumes a simple parallelization, i.e. many processes can read
    the same chain file, each one filtering out 1 target chain and many query chains.
    This always assumes that there are some blocks in the chain file,
    otherwise raises AssertionError.

    :param fobj: object supporting line iteration/read
     :type: file-like object opened in text mode
    :param tselect: compiled regex object to check for allowed target chromosomes
     :type: None or re.regex object
    :param qselect: compiled regex object to check for allowed query chromosomes
     :type: None or re.regex object
    :param read_num: read only this many chains
     :type: int, default 0
    :param min_size: skip chains with size smaller than this number
     :type: int, default 0
    :param min_score, default 0
     :type: skip chains with score lower than this number
    :return:
    """
    tchrom, qchrom = None, None
    tstrand, qstrand = '+', None
    trun, qrun = -1, -1
    exp_tend, exp_qend = -1, -1
    # bring reader function to local scope
    read_head = read_chain_header
    read_block = read_chain_block
    tskip = fnt.partial(check_skip, *(tselect,))
    qskip = fnt.partial(check_skip, *(qselect,))
    skip = False
    chain_id = ''
    chain_score = 0
    cc = 0  # count chains
    bc = 0  # count blocks per chain
    for line in fobj:
        if line.strip():
            if line.startswith('chain'):
                assert trun == exp_tend or exp_tend == -1, 'Missed expected block end: {} vs {}'.format(trun, exp_tend)
                assert qrun == exp_qend or exp_qend == -1, 'Missed expected block end: {} vs {}'.format(qrun, exp_qend)
                if 0 < read_num <= cc:
                    # read user specified number of chains,
                    # ignore the rest
                    break
                bc = 0  # reset block count for new chain
                chain = read_head(line)
                assert chain.tstrand == '+', 'Reverse target chain is unexpected: {}'.format(line)
                chain_id = chain.chainid
                chain_score = chain.score
                chain_size = chain.tend - chain.tstart
                if chain_score < min_score or chain_size < min_size:
                    # the size check rests on the understanding that the size of a chain
                    # is the same in both target and query - note that this is the
                    # length including all gaps, NOT the number of aligning bases
                    skip = True
                    continue
                tchrom = chain.tname
                skip = tskip(tchrom)
                if skip:
                    continue
                qchrom = chain.qname
                skip = qskip(qchrom)
                if skip:
                    continue
                cc += 1
                trun = chain.tstart
                exp_tend = chain.tend
                qstrand = chain.qstrand
                if qstrand == '+':
                    qrun = chain.qstart
                    exp_qend = chain.qend
                elif qstrand == '-':
                    qrun = chain.qsize - chain.qend
                    exp_qend = chain.qsize - chain.qstart
                else:
                    raise ValueError('Unknown query strand in chain file: {}'.format(line))
            else:
                if skip:
                    continue
                try:
                    blocksize, tgap, qgap = read_block(line)
                    bc += 1
                    yield Block(tchrom, trun, trun + blocksize, tstrand,
                                qchrom, qrun, qrun + blocksize, qstrand,
                                chain_id + '.' + str(bc), chain_score, blocksize)
                    trun += blocksize + tgap
                    qrun += blocksize + qgap
                except ValueError:
                    # end of chain - last line just contains
                    # size of last block
                    blocksize = int(line)
                    bc += 1
                    yield Block(tchrom, trun, trun + blocksize, tstrand,
                                qchrom, qrun, qrun + blocksize, qstrand,
                                chain_id + '.' + str(bc), chain_score, blocksize)
                    trun += blocksize
                    qrun += blocksize
    try:
        tpat = tselect.pattern
    except AttributeError:
        # in case no pattern was specified, tselect is None
        tpat = 'all'
    try:
        qpat = qselect.pattern
    except AttributeError:
        qpat = 'all'
    assert cc > 0, 'No chain read from file for target/query selectors: {}'.format(tpat, qpat)
    assert bc > 0 or skip, 'No aln. blocks read from chain file for target/query selectors: {} / {}'.format(tpat, qpat)
    return


def assert_valid_region(block, blocktype):
    """
    :param block:
    :param blocktype:
    :return:
    """
    assert block.tstart < block.tend, 'Invalid target region: {}'.format(block)
    assert block.qstart < block.qend, 'Invalid query region: {}'.format(block)
    tlen = block.tend - block.tstart
    qlen = block.qend - block.qstart
    assert tlen == qlen, 'Non-symmetric block detected: {}\n{} vs {} ({})'.format(block, tlen, qlen, blocktype)
    return


def chain_symmetry_filter(args):
    """
    This function is used to filter a chain file that is already single
    coverage for the target species but not for the query. It produces
    symmetric blocks that, in total, represent a symmetric one-to-one
    mapping between target and query assembly.
    The chain file has to be sorted in descending order by chain score.
    In case of overlapping chains in the query, the chain with the lower
    score will be adapted (the block edges rearranged) to cover a
    [close to] maximal set of genomic bases in the query.
    The iteration through the chain file is query-centric, i.e., due to
    the assumption that the chains are already single coverage for the
    target, the iterator selects all target chromosomes that map to the
    same query chromosome. This allows for easy parallel processing
    of the file.
    :param args:
    :return: fixed return value 0
    """
    chromsize = chromsize_from_chain(args.chainfile, args.chrom, target=False)
    # save query regions in form of masked arrays
    qchrom = msk.masked_array(np.arange(chromsize, dtype=np.int32))
    qchrom.mask = False

    # variables to collect some statistics
    blocks_total = 0
    blocks_contained = 0
    blocks_symmetric = 0
    blocks_bridging = 0
    blocks_scattered = 0
    discard_bridge = 0
    coverage_bridge = 0
    used_fragments = 0
    discard_fragments = 0
    coverage_fragments = 0
    # output buffer for symmetrical blocks
    block_buffer = io.StringIO()
    write_block = write_block_to_buffer
    opn, mode = text_file_mode(args.chainfile, read=True)
    valid_block = assert_valid_region
    with opn(args.chainfile, mode) as infile:
        chainit = get_chain_iterator(infile, qselect=re.compile(args.chrom + '$'),
                                     min_size=args.minsize, min_score=args.minscore)
        last_score = sys.maxsize
        for block in chainit:
            assert block.score <= last_score, 'Chains not sorted in descending order for score: {}'.format(block)
            last_score = block.score
            blocks_total += 1
            qs, qe = block.qstart, block.qend
            pos_masked = msk.count_masked(qchrom[qs:qe])
            if pos_masked == block.size:
                # all values masked -> block is contained
                # in block with higher score
                blocks_contained += 1
            elif pos_masked == 0:
                # no values masked -> use block as-is
                qchrom[qs:qe].mask = 1
                write_block(block, block_buffer)
                blocks_symmetric += 1
            else:
                # block overlaps with other blocks in query
                # I expect this to occur mostly when a chain
                # fills a gap between two higher order chains, i.e.,
                # left and right end can be simply adjusted
                unmask_slices = msk.clump_unmasked(qchrom[qs:qe])
                if len(unmask_slices) == 1:
                    # the common case, i.e., block overlaps
                    # at left or right end with higher order blocks
                    blocks_bridging += 1
                    bridge = unmask_slices[0]
                    left_offset = bridge.start
                    bridge_length = bridge.stop - bridge.start
                    if bridge_length < args.minsize:
                        discard_bridge += 1
                        continue
                    block = block._replace(**{'tstart': block.tstart + left_offset,
                                              'tend': block.tstart + left_offset + bridge_length,
                                              'qstart': block.qstart + left_offset,
                                              'qend': block.qstart + left_offset + bridge_length,
                                              'size': -1})
                    valid_block(block, 'Bridge')
                    coverage_bridge += bridge_length
                    pos_masked = msk.count_masked(qchrom[block.qstart:block.qend])
                    assert pos_masked == 0, 'Bridge edge adjustment failed by {}: {}'.format(pos_masked, block)
                    assert block.qstart < block.qend, 'Bridge edge adjustment created invalid region: {}'.format(block)
                    write_block(block, block_buffer)
                else:
                    # rather rare case unless no selection based
                    # on score/size is done and many small chains are
                    # still present in the filtering process
                    blocks_scattered += 1
                    # sort all fragments by size
                    fragments = sorted(unmask_slices, key=lambda x: x.stop - x.start, reverse=True)
                    for idx, f in enumerate(fragments):
                        frag_len = f.stop - f.start
                        if frag_len < args.minsize:
                            discard_fragments += len(fragments[idx:])
                            break
                        used_fragments += 1
                        left_offset = f.start
                        tmp_block = block._replace(**{'tstart': block.tstart + left_offset,
                                                      'tend': block.tstart + left_offset + frag_len,
                                                      'qstart': block.qstart + left_offset,
                                                      'qend': block.qstart + left_offset + frag_len,
                                                      'size': -1})
                        valid_block(tmp_block, 'Fragment')
                        coverage_fragments += frag_len
                        pos_masked = msk.count_masked(qchrom[tmp_block.qstart:tmp_block.qend])
                        assert pos_masked == 0, 'Fragment edge adjustment failed by {}: {}'.format(pos_masked, tmp_block)
                        assert tmp_block.qstart < tmp_block.qend, 'Fragment edge adjustment created invalid region: {}'.format(tmp_block)
                        write_block(tmp_block, block_buffer)
                    continue
    # all the int() because JSON cannot serialize numpy.int*
    chrom_cov = int(qchrom.mask.sum())
    metadata = {'chainfile': os.path.basename(args.chainfile), 'chromosome': args.chrom,
                'chrom_size': int(chromsize), 'chrom_cov': chrom_cov, 'chrom_cov_pct': round(chrom_cov / chromsize, 4),
                'blocks_total': int(blocks_total),
                'raw_symmetric_cov': int(chrom_cov - coverage_bridge - coverage_fragments),
                'raw_symmetric_cov_pct': round((chrom_cov - coverage_bridge - coverage_fragments) / chrom_cov, 4),
                'blocks_contained': int(blocks_contained), 'blocks_contained_pct': round(blocks_contained / blocks_total, 4),
                'blocks_symmetric': int(blocks_symmetric), 'blocks_symmetric_pct': round(blocks_symmetric / blocks_total, 4),
                'blocks_bridging': int(blocks_bridging), 'blocks_bridging_pct': round(blocks_bridging / blocks_total, 4),
                'blocks_scattered': int(blocks_scattered), 'blocks_scattered_pct': round(blocks_scattered / blocks_total, 4),
                'bridges_discarded': int(discard_bridge), 'bridges_discarded_pct': round(discard_bridge / max(1, blocks_bridging), 4),
                'bridge_coverage': int(coverage_bridge),
                'fragments_used': int(used_fragments), 'fragments_discarded': int(discard_fragments),
                'fragments_used_pct': round(used_fragments / max(1, used_fragments + discard_fragments), 4),
                'fragments_discarded_pct': round(discard_fragments / max(1, used_fragments + discard_fragments), 4)}
    ofext = identify_file_extension(args.outputfile)
    stat_out = args.outputfile.replace(ofext, 'json')
    with open(stat_out, 'w') as outf:
        _ = js.dump(metadata, outf, indent=1, sort_keys=True)
    opn, mode = text_file_mode(args.outputfile, read=False)
    with opn(args.outputfile, mode) as outf:
        _ = outf.write(block_buffer.getvalue())
    return 0


# ==========================================
# outdated functions, currently not used
# ==========================================


def build_query_chrom_masks(chromsizes):
    """
    :param chromsizes:
    :return:
    """
    chrommasks = dict()
    selector = re.compile('chr[0-9]+$')
    with open(chromsizes, 'r') as infile:
        for line in infile.readlines():
            c, s = line.strip().split()[:2]
            if selector.match(c) is not None:
                chrommasks[c] = msk.masked_array(np.arange(int(s)), dtype=np.int32)
                chrommasks[c].mask = False
    assert chrommasks, 'No chromosomes were masked'
    return chrommasks


def get_mapfile_iterator(fobj):
    """
    :param fobj:
    :return:
    """
    pml = parse_map_line
    while 1:
        lines = fobj.readlines(32768)
        if not lines:
            break
        for line in lines:
            yield pml(line)
    return


def mask_query_blocks(mapfile, qmasks):
    """
    :param mapfile:
    :param qmasks:
    :return:
    """
    qc, qs, qe = 5, 6, 7
    with gz.open(mapfile, 'rt') as infile:
        lineiter = get_mapfile_iterator(infile)
        for line in lineiter:
            assert msk.count_masked(qmasks[line[qc]][line[qs]:line[qe]]) == 0,\
                'Overlapping query blocks at position: {}'.format(line)
            qmasks[line[qc]][line[qs]:line[qe]].mask = 1
    return qmasks


def parse_map_line(line):
    """
    :param line:
    :return:
    """
    tchrom, tstart, tend, tstrand, blockid, qchrom, qstart, qend, qstrand = line.strip().split()
    return tchrom, int(tstart), int(tend), tstrand, blockid.split('.')[0], qchrom, int(qstart), int(qend), qstrand


def check_block_merge(blocks, chrom):
    """ A test run of this function on hg19 / mm9 indicates that a straightforward
    merge of consecutive blocks is not possible (due to varying gap sizes)

    Hence, this function is deprecated / useless and left here just as a reminder

    :param blocks:
    :param chrom:
    :return:
    """
    if len(blocks) < 2:
        return blocks, False
    outblocks = []
    merged = False
    for idx in range(len(blocks)):
        subset = blocks[idx:]
        if len(subset) == 1:
            outblocks.append(subset[0])
            break
        elif (subset[-1][2] - subset[0][1]) == (subset[-1][7] - subset[0][6]) and \
                sum(s[7] - s[6] for s in subset) == chrom[subset[0][6]:subset[0][7]].mask.sum():
            # same length in target and query and no intervening blocks in query
            outblocks.append((subset[0][0], subset[0][1], subset[-1][2], subset[0][3],
                              subset[0][4],
                              subset[0][5], subset[0][6], subset[-1][7], subset[0][8]))
            chrom[subset[0][6]:subset[-1][7]].mask = 1
            merged = True
            break
        else:
            outblocks.append(blocks[idx])
    assert outblocks, 'No output generated for input of length {}'.format(len(blocks))
    return outblocks, merged


def check_gap_extendable(covratio, mincov, checkext, blocks, chrom):
    """
    :param blocks:
    :param chrom:
    :return:
    """

    ts = blocks[0][1]
    qs = blocks[0][6]
    te = blocks[-1][2]
    qe = blocks[-1][7]
    max_len = max(te - ts, qe - qs)
    tcov = sum([b[2] - b[1] for b in blocks])
    qcov = sum([b[7] - b[6] for b in blocks])
    assert tcov == qcov, 'Coverage not equal ({} vs {}) for blocks: {}'.format(tcov, qcov, blocks)
    if tcov < mincov:
        return blocks
    outblocks = []
    for ext in checkext:
        new_len = max_len // ext * ext + ext
        if chrom[qs:qs+new_len].mask.sum() == qcov:
            # no overlapping blocks from other chains in this gap, extend
            if qcov / (qs+new_len - qs) < covratio:
                continue
            outblocks = [(blocks[0][0], ts, ts+new_len, blocks[0][3],
                          blocks[0][4],
                          blocks[0][5], qs, qs+new_len, blocks[0][8])]
            chrom[qs:qs+new_len].mask = 1
            break
    if not outblocks:
        if chrom[qs:qs+max_len].mask.sum() == qcov:
            if not qcov / (qs+max_len - qs) < covratio:
                outblocks = [(blocks[0][0], ts, ts+max_len, blocks[0][3],
                              blocks[0][4],
                              blocks[0][5], qs, qs+max_len, blocks[0][8])]
                chrom[qs:qs+max_len].mask = 1
            else:
                outblocks = blocks
        else:
            outblocks = blocks
    return outblocks


def check_gap_extend(gapsize, covratio, mincov, checkext, blocks, chrom):
    """
    :param blocks:
    :param chrom:
    :return:
    """
    check_pos = []
    for idx, block in enumerate(blocks):
        try:
            if blocks[idx+1][1] - block[2] >= gapsize and \
                    blocks[idx+1][6] - block[7] >= gapsize:
                check_pos.append(idx)
        except IndexError:
            continue
    if not check_pos:
        return blocks
    outblocks = []
    start_pos = 0
    check_ext = fnt.partial(check_gap_extendable, *(covratio, mincov, checkext))
    for inc_end in check_pos:
        s, e = start_pos, inc_end + 1
        subset = blocks[s:e]
        if len(subset) == 1:
            # isolated singleton
            outblocks.extend(subset)
            start_pos = inc_end + 1
            continue
        assert subset, 'Invalid subset selection of blocks - indices {} to {}'.format(s, e)
        outblocks.extend(check_ext(subset, chrom))
        start_pos = inc_end + 1
    if start_pos < len(blocks):
        outblocks.extend(blocks[start_pos:])
    assert outblocks, 'No output created for input blocks: {}'.format(blocks)
    assert len(outblocks) == len(set(outblocks)), 'Duplicate blocks created: {}'.format(outblocks)
    return outblocks


def symm_extend_blocks(args):
    """
    :param args:
    :return:
    """
    mf_ext = get_file_extension(args, 'mapfile', True)
    maskfile = args.mapfile.replace(mf_ext, '.pck')
    if not os.path.isfile(maskfile):
        cm = build_query_chrom_masks(args.qsizes)
        qmasks = mask_query_blocks(args.mapfile, cm)
        pck.dump(qmasks, open(maskfile, 'wb'))
    else:
        maskfile = args.mapfile.replace(mf_ext, '.pck')
        qmasks = pck.load(open(maskfile, 'rb'))
    raw_cov = summarize_coverage(qmasks)
    out_ext = get_file_extension(args, 'output', True)
    raw_cov_stats = collect_block_statistics(qmasks, args.output.replace(out_ext, '.rawcov.pck'))
    pml = parse_map_line
    ext = 0
    maps = 0
    # empty output file in case there is
    # an incomplete file from a previous run
    with gz.open(args.output, 'wt') as outf:
        _ = outf.write('')
    outbuffer = []
    all_mrg_cov = 0
    check_gap = fnt.partial(check_gap_extend, *(args.gapsize, args.covratio, args.mincov, args.checkext))
    with gz.open(args.mapfile, 'rt') as infile:
        chainbuffer = []
        chain = ''
        chrom = ''
        for line in infile:
            l = pml(line)
            maps += 1
            if l[4] != chain:
                try:
                    chaincov = sum([b[2] - b[1] for b in chainbuffer])
                    outblocks = check_gap(chainbuffer, qmasks[chrom])
                    mrgcov = sum([b[2] - b[1] for b in outblocks])
                    all_mrg_cov += mrgcov
                    assert mrgcov >= chaincov, 'Coverage shrink - lost blocks: {}\n\n{}'.format(chainbuffer, outblocks)
                    outbuffer.extend(outblocks)
                    if len(outblocks) != len(chainbuffer):
                        ext += 1
                except KeyError:
                    pass
                chainbuffer = [l]
                chrom = l[5]
                chain = l[4]
                if len(outbuffer) > 1000000:
                    with gz.open(args.output, 'at') as outf:
                        for b in outbuffer:
                            _ = outf.write('\t'.join(map(str, b)) + '\n')
                    outbuffer = []
            else:
                chainbuffer.append(l)
        chaincov = sum([b[2] - b[1] for b in chainbuffer])
        outblocks = check_gap(chainbuffer, qmasks[chrom])
        mrgcov = sum([b[2] - b[1] for b in outblocks])
        all_mrg_cov += mrgcov
        assert mrgcov >= chaincov, 'Coverage shrink - lost blocks: {}\n\n{}'.format(chainbuffer, outblocks)
        outbuffer.extend(outblocks)
        if len(outblocks) != len(chainbuffer):
            ext += 1
    with gz.open(args.output, 'at') as outf:
        for b in outbuffer:
            _ = outf.write('\t'.join(map(str, b)) + '\n')
        outbuffer = []
    mrg_cov = summarize_coverage(qmasks)
    mrg_cov_stats = collect_block_statistics(qmasks, args.output.replace(out_ext, '.mrgcov.pck'))
    assert mrg_cov == all_mrg_cov, 'Coverage mismatch between mask and regions: {} - {}'.format(mrg_cov, all_mrg_cov)
    md = dict()
    md['num_maps'] = maps
    md['num_ext_blocks'] = ext
    md['raw_cov_bp'] = int(raw_cov)
    md['mrg_cov_bp'] = int(mrg_cov)
    md['mrg_ratio'] = round(mrg_cov / raw_cov, 5)
    md['input_mapfile'] = args.mapfile
    md['output_mapfile'] = args.output
    md['raw_cov_stats'] = raw_cov_stats
    md['mrg_cov_stats'] = mrg_cov_stats
    md['param_min_cov'] = args.mincov
    md['param_cov_ratio'] = args.covratio
    md['param_gap_size'] = args.gapsize
    md['param_check_ext'] = args.checkext
    with open(args.output.replace('.tsv.gz', '.json'), 'w') as metadata:
        js.dump(md, metadata, indent=1)
    assert not outbuffer, 'Out buffer not empty'
    return 0


def collect_block_statistics(qmasks, outfile):
    """
    :param qmasks:
    :return:
    """
    gapdist = []
    sizedist = []
    for chrom, masked in qmasks.items():
        sizes = msk.clump_masked(masked)
        sizedist.append([s.stop - s.start for s in sizes])
        tmp = msk.masked_array(masked.data, msk.logical_not(masked.mask))
        leftend, rightend = msk.flatnotmasked_edges(tmp)
        # ignore the regions to the chromosome boundaries
        gaps = msk.clump_masked(tmp[leftend:rightend+1])
        gapdist.append([g.stop - g.start for g in gaps])
    gapdist = np.concatenate(gapdist)
    sizedist = np.concatenate(sizedist)
    stats = {'gaps': gapdist, 'blocks': sizedist}
    with open(outfile, 'wb') as outf:
        pck.dump(stats, outf)
    return outfile


def summarize_coverage(chroms):
    """
    :param chroms:
    :return:
    """
    cov = 0
    for k, v in chroms.items():
        cov += v.mask.sum()
    return cov


def _block_to_bed_line(block):
    """ tchrom, trun, trun + size, tstrand, qchrom, qrun, qrun + size, qstrand, chain_id + '.' + str(bc), chain_score
    :param block:
    :return:
    """
    return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(*block)


def chain_to_bed_unbuffered(args):
    """
    :param args:
    :return:
    """
    chainit = get_chain_iterator(gz.open(args.chainfile, 'rt'), min_score=1000, min_size=args.mincov)
    bbl = _block_to_bed_line
    for block in chainit:
        sys.stdout.write(bbl(block))
    return 0


def chain_to_bed(args):
    """
    :param args:
    :return:
    """
    outbuffer = io.StringIO()
    chars = 0
    chainit = get_chain_iterator(gz.open(args.chainfile, 'rt'), min_score=1000, min_size=args.mincov)
    bbl = _block_to_bed_line
    for block in chainit:
        chars += outbuffer.write(bbl(block))
        if chars > 1000000:
            sys.stdout.write(outbuffer.getvalue())
            outbuffer = io.StringIO()
    sys.stdout.write(outbuffer.getvalue())
    return 0


def _target_map_to_bedgraph(line):
    """
    :param line:
    :return:
    """
    cols = line.strip().split()
    return '{}\t{}\t{}\t1\n'.format(cols[0], cols[1], cols[2])


def _query_map_to_bedgraph(line):
    """ Note to self: even though the strand can be minus,
    all coordinates are given in plus orientation by construction
    :param line:
    :return:
    """
    cols = line.strip().split()
    return '{}\t{}\t{}\t1\n'.format(cols[5], cols[6], cols[7])


def map_to_bedgraph(mapfile, query):
    """
    :param mapfile:
    :param query:
    :return:
    """
    if not query:
        readline = _target_map_to_bedgraph
    else:
        readline = _query_map_to_bedgraph
    with gz.open(mapfile, 'rt') as infile:
        for line in infile:
            sys.stdout.write(readline(line))
    return 0


def swap_map(mapfile):
    """
    :param mapfile:
    :return:
    """
    inverse = op.itemgetter(*(5, 6, 7, 3, 4, 0, 1, 2, 8))
    with gz.open(mapfile, 'rt') as infile:
        for line in infile:
            if not line:
                continue
            cols = line.strip().split()
            swl = '\t'.join(inverse(cols)) + '\n'
            sys.stdout.write(swl)
    return


def get_file_extension(args, which, compressed):
    """
    :param args:
    :param which:
    :param compressed:
    :return:
    """
    fpath = getattr(args, which)
    if compressed:
        if not fpath.endswith('.gz') or fpath.endswith('.gzip'):
            fpath += '.gz'
        fileext = '.' + '.'.join(fpath.split('.')[-2:])
    else:
        fileext = '.' + fpath.split('.')[-1]
    setattr(args, which, fpath)
    return fileext


if __name__ == '__main__':
    try:
        args = parse_command_line()
        cmd_select = {'symmfilt': chain_symmetry_filter}
        run_cmd = cmd_select[args.task]
        _ = run_cmd(args)
    except Exception:
        trb.print_exc()
        sys.exit(1)
    else:
        sys.exit(0)
