#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import gzip as gz

import pandas as pd


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--chain-file', '-c', type=str, dest='chainfile')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')

    args = parser.parse_args()
    return args


def parse_chain_header(line):
    """
    Reminder from genome.ucsc.edu/goldenpath/help/chain.html :

    chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id

    The alignment start and end positions are represented as zero-based half-open intervals.

    When the strand value is "-", position coordinates are listed in terms of the reverse-complemented sequence.

    :param line:
    :return:
    """
    _, score, tchrom, tsize, tstrand, tstart, tend, qchrom, qsize, qstrand, qstart, qend, chainid = line.split()
    assert tstrand == '+', 'Reverse target strand not expected'
    if qstrand == '-':
        tmpstart = int(qsize) - int(qend)
        qend = int(qsize) - int(qstart)
        qstart = tmpstart
        assert qstart < qend, 'Malformed reverse query block: {} / {} from {}'.format(qstart, qend, line)
    tstart, tend, qstart, qend = int(tstart), int(tend), int(qstart), int(qend)
    tlen = tend - tstart
    qlen = qend - qstart
    res = int(score), tchrom, tsize, '+', tstart, tend,\
        qchrom, int(qsize), qstrand, qstart, qend, int(chainid), tlen, qlen
    return res


def main():
    """
    :return:
    """
    chain_columns = ['score', 'tchrom', 'tsize', 'tstrand', 'tstart', 'tend']
    chain_columns += ['qchrom', 'qsize', 'qstrand', 'qstart', 'qend', 'chainid', 'tlength', 'qlength']
    args = parse_command_line()
    chains = []
    with gz.open(args.chainfile, 'rt') as chf:
        for line in chf:
            if line.startswith('chain'):
                chain = parse_chain_header(line)
                chains.append(chain)
    chainlut = pd.DataFrame(chains, columns=chain_columns)
    with pd.HDFStore(args.outputfile, 'w', complib='blosc', complevel=9) as hdf:
        hdf.put('chainlut', chainlut, format='fixed')
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
