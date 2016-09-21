# coding=utf-8

import os as os
import gzip as gz
import numpy as np
import fnmatch as fnm


def collect_full_paths(rootdir, pattern):
    """
    :param rootdir:
    :param pattern:
    :return:
    """
    all_files = []
    for root, dirs, files in os.walk(rootdir):
        if files:
            filt = fnm.filter(files, pattern)
            for f in filt:
                all_files.append(os.path.join(root, f))
    return all_files


def read_chromsizes(fp):
    """
    :param fp:
    :return:
    """
    bounds = dict()
    with open(fp, 'r') as infile:
        for line in infile:
            if not line.strip():
                continue
            name, size = line.strip().split()
            bounds[name] = int(size)
    assert bounds, 'No boundaries read from file {}'.format(fp)
    return bounds


def open_comp(fp, read=True):
    """ unnecessary in Python 3.5+, open in 'rt' mode
    :param fp:
    :param read:
    :return:
    """
    if fp.endswith('.gz') or fp.endswith('.gzip'):
        if read:
            return gz.open, 'rb', lambda x: x.decode('ascii').strip()
        else:
            return gz.open, 'wb', lambda x: x.encode('ascii')
    else:
        if read:
            return open, 'r', lambda x: x.strip()
        else:
            return open, 'w', lambda x: x


def check_bounds(chrom, start, end, bounds, fname):
    """
    :param chrom:
    :param start:
    :param end:
    :param bounds:
    :return:
    """
    good = False
    try:
        assert int(end) <= bounds[chrom],\
            'Region {} out of bounds: {} - {} - {} in file {}'.format(chrom, start, end, fname)
        _ = int(start)
        good = True
    except ValueError:
        raise ValueError('Non-integer coordinates {} - {} - {} in file {}'.format(chrom, start, end, fname))
    except KeyError:
        # chromosome not in check file, skip
        pass
    return good


def collect_region_stats(inputfile, outputfile):
    """
    :param inputfile:
    :param outputfile:
    :return:
    """
    opn, mode, enc = open_comp(inputfile, True)
    lens = []
    with opn(inputfile, mode) as infile:
        for line in infile:
            line = enc(line)
            if not line or line.startswith('#'):
                continue
            cols = line.split()
            lens.append(int(cols[2]) - int(cols[1]))
    assert lens, 'No regions read from file {}'.format(inputfile)
    percentiles = [0, 5, 25, 50, 75, 95, 100]
    header = ['name', 'num_regions']
    for pct in percentiles:
        header.append('pct_' + str(pct).zfill(3))
    scores = np.percentile(lens, percentiles)
    infos = [os.path.basename(inputfile), '{}'.format(len(lens))] + [str(int(s)) for s in scores]
    opn, mode, enc = open_comp(outputfile, False)
    with opn(outputfile, mode) as outfile:
        _ = outfile.write(enc('\t'.join(header) + '\n'))
        _ = outfile.write(enc('\t'.join(infos) + '\n'))
    return outputfile
