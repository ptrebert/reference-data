#!/usr/bin/env python3
# coding=utf-8

import sys as sys
import io as io
import gzip as gz
import argparse as argp
import csv as csv
import operator as op
import traceback as trb
import twobitreader as tbr


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


def get_seqrc(sequence):
    """
    :param sequence:
    :return:
    """
    rcmap = {'A': 'T', 'a': 't',
             'C': 'G', 'c': 'g',
             'G': 'C', 'g': 'c',
             'T': 'A', 't': 'a',
             'N': 'N', 'n': 'n'}
    return ''.join(map(lambda x: rcmap[x], reversed(list(sequence))))


def check_need_rc(seq, is_minus):
    """
    :param seq:
    :param is_minus:
    :return:
    """
    if is_minus:
        seq = get_seqrc(seq)
    return seq.upper()


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, dest='input', required=True)
    parser.add_argument('--seq-file', type=str, dest='seqfile', required=True)
    parser.add_argument('--output', '-o', type=str, dest='output', required=True)
    parser.add_argument('--coord-base', type=int, default=0, choices=[0, 1], dest='coordbase')
    args = parser.parse_args()
    return args


def read_transcript_bed(fp, offset):
    """
    :param fp:
    :param offset:
    :return:
    """
    opn, mode = det_open_mode(fp, True)
    coords = []
    get_coords = op.itemgetter(*('#chrom', 'start', 'end', 'name', 'strand'))
    with opn(fp, mode) as infile:
        rows = csv.DictReader(infile, delimiter='\t')
        for r in rows:
            trans = list(get_coords(r))
            trans[1] = int(trans[1])
            trans[2] = int(trans[2]) + offset
            coords.append(tuple(trans))
    return coords


def build_transcriptome(transcripts, seqfile, outfile):
    """
    :param transcripts:
    :param seqfile:
    :param outfile:
    :return:
    """
    outbuffer = io.StringIO('')
    current_chrom = None
    chrom_seq = ''
    wg_seq = tbr.TwoBitFile(seqfile)
    opn, mode = det_open_mode(outfile, False)
    for t in transcripts:
        if t[0] != current_chrom:
            current_chrom = t[0]
            chrom_seq = wg_seq[current_chrom]
            with opn(outfile, mode) as outf:
                _ = outf.write(outbuffer.getvalue())
            mode = mode.replace('w', 'a')
            outbuffer = io.StringIO()
        minus = t[4] == '-'
        trans_seq = check_need_rc(chrom_seq[t[1]:t[2]], minus)
        outbuffer.write('>{}\n'.format(t[3]))
        for idx in range(0, len(trans_seq), 80):
            outbuffer.write(trans_seq[idx:idx+80] + '\n')
    return 0


if __name__ == '__main__':
    try:
        args = parse_command_line()
        transcripts = read_transcript_bed(args.input, args.coordbase)
        _ = build_transcriptome(transcripts, args.seqfile, args.output)
        sys.exit(0)
    except Exception as err:
        trb.print_exc()
        sys.stderr.write('\nError: {}\n'.format(err))
        sys.exit(1)
