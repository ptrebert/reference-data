#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import collections as col
import gzip as gz
import csv as csv
import pickle as pck
import numpy as np


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, nargs='+', dest='inputfiles')
    parser.add_argument('--model', '-m', type=str, dest='model')
    parser.add_argument('--output', '-o', type=str, dest='output')
    args = parser.parse_args()
    return args


def get_model_header(fpath):
    """
    :param fpath:
    :return:
    """
    with gz.open(fpath, 'rt') as modelfile:
        header = modelfile.readline().strip().strip('#').split()
    return header


def process_overlap_files(inputfiles, header):
    """
    :param inputfiles:
    :param header:
    :return:
    """
    cage_cov = {'+': dict(), '-': dict()}
    feat_count = col.Counter()
    skip_partial = 0
    for inpf in inputfiles:
        with open(inpf, 'r', newline='') as infile:
            rows = csv.DictReader(infile, fieldnames=header, delimiter='\t')
            for idx, row in enumerate(rows):
                assert row['cage_strand'] == row['strand'], 'Feature overlap has to be strand specific: {}'.format(row)
                feature_type = row['feature']
                cage_length = int(row['cage_end']) - int(row['cage_start'])
                if cage_length >= 1000:
                    continue  # this is much less than 1%
                cage_strand = row['cage_strand']
                cs = int(row['cage_start'])
                ce = int(row['cage_end'])
                feat_count[feature_type] += 1
                try:
                    if cage_strand == '-':
                        ft5p = int(row['end'])
                        offset = ft5p - 500
                        new_cs = cs - offset
                        new_ce = ce - offset
                    else:
                        ft5p = int(row['start'])
                        offset = ft5p - 500
                        new_cs = cs - offset
                        new_ce = ce - offset
                    assert new_cs + cage_length == new_ce, 'Cage length mismatch: {} / {} --- {}'.format(new_cs, new_ce, row)
                    if new_cs < 0 or new_cs + cage_length > 1000:
                        skip_partial += 1
                        continue
                    cage_cov[cage_strand][feature_type][new_cs:new_ce] += 1
                except KeyError:
                    cage_cov[cage_strand][feature_type] = np.zeros(1001, dtype=np.int32)
                    cage_cov[cage_strand][feature_type][new_cs:new_ce] += 1
    for feat, region in cage_cov['-'].items():
        rev_region = region[::-1]
        cage_cov['+'][feat] += rev_region
    return cage_cov['+'], feat_count


def main():
    """
    :return:
    """
    args = parse_command_line()
    model_header = get_model_header(args.model)
    cage_header = ['cage_chrom', 'cage_start', 'cage_end', 'cage_strand']
    full_header = cage_header + model_header
    cov, ftcount = process_overlap_files(args.inputfiles, full_header)
    with open(args.output, 'wb') as dump:
        pck.dump({'values': cov, 'counts': ftcount}, dump)
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
