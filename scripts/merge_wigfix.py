#!/usr/bin/env python
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import gzip as gz
import multiprocessing as mp

import numpy as np
import pandas as pd


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, nargs='+', dest='inputfiles')
    parser.add_argument('--chrom', '-c', type=str, dest='chromsizes')
    parser.add_argument('--workers', '-w', type=int, default=4, dest='workers')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')
    args = parser.parse_args()
    return args


def parse_chromsizes(fpath):
    """
    :param fpath:
    :return:
    """
    chroms = dict()
    with open(fpath, 'r') as infile:
        for line in infile:
            c, s = line.strip().split()
            chroms[c] = int(s)
    return chroms


def parse_wiggle_header(header):
    """
    :param header:
    :return:
    """
    _, chrom, start, step = header.strip().split()
    foo, chrom = chrom.split('=')
    assert foo == 'chrom', 'Unexpected wiggle header: {}'.format(header)
    bar, start = start.split('=')
    assert bar == 'start', 'Unexpected wiggle header: {}'.format(header)
    loo, step = step.split('=')
    assert loo == 'step', 'Unexpected wiggle header: {}'.format(header)
    # NB: wiggle coordinates start at 1, make 0-based here
    return chrom, int(start) - 1, int(step)


def process_fixed_wiggle_track(params):
    """
    :param params:
    :return:
    """
    fpath, chroms, outputfile, lock = params
    with gz.open(fpath, 'rt') as infile:
        header = infile.readline()
        chrom, start, step = parse_wiggle_header(header)
        carray = np.zeros(chroms[chrom], dtype=np.float32)
        for line in infile:
            try:
                carray[start:start+step] = float(line)
                start += step
            except ValueError:
                if not line.strip():
                    continue
                _, start, step = parse_wiggle_header(line)
                continue
    with lock:
        with pd.HDFStore(outputfile, 'a', complib='blosc', complevel=9) as hdf:
            hdf.put(chrom, pd.Series(carray, dtype=np.float32), format='fixed')
            hdf.flush()
    return chrom


def main(args):
    """
    :param args:
    :return:
    """
    chroms = parse_chromsizes(args.chromsizes)
    with pd.HDFStore(args.outputfile, 'w', complib='blosc', complevel=9) as hdf:
        pass

    with mp.Manager() as mng:
        file_lock = mng.Lock()
        params = [(fp, chroms, args.outputfile, file_lock) for fp in args.inputfiles]
        workers = mp.Pool(args.workers)
        res = workers.imap_unordered(process_fixed_wiggle_track, params)
        for chrom in res:
            #print('{} done'.format(chrom))
            pass
    return


if __name__ == '__main__':
    try:
        args = parse_command_line()
        main(args)
    except Exception as err:
        trb.print_exc()
        raise err
    else:
        sys.exit(0)
