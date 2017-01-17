#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--map-1', '-m1', type=str, dest='map1')
    parser.add_argument('--map-2', '-m2', type=str, dest='map2')
    parser.add_argument('--map-ref', '-mr', type=str, choices=['target', 'query'], dest='mapref')
    parser.add_argument('--output', '-o', type=str, dest='output')
    args = parser.parse_args()
    return args


def compare_chrom_mask(c1, c2):
    """
    Masks are given as bool arrays with
    True/1/masked = not conserved
    False/0/unmasked = conserved
    :param c1: np.array
    :param c2: np.array
    :return:
    """
    assert c1.size == c2.size, 'Expected chromosome masks of same size: {} vs {}'.format(c1.size, c2.size)
    raw_cov_c1 = c1.size - c1.sum()
    raw_cov_c2 = c2.size - c2.sum()
    joined_cov = c1.size - (c1 | c2).sum()
    return raw_cov_c1, raw_cov_c2, joined_cov


def main():
    """
    :return:
    """
    args = parse_command_line()


if __name__ == '__main__':
    try:
        main()
    except Exception as err:
        trb.print_exc(file=sys.stderr)
        sys.stderr.write('\nError: {}\n'.format(str(err)))
        sys.exit(1)
    else:
        sys.exit(0)
