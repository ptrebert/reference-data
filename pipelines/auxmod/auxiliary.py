# coding=utf-8

import os as os
import gzip as gz
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

