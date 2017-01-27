# coding=utf-8

import os as os
import io as io
import re as re
import gzip as gz
import functools as fnt

from pipelines.auxmod.auxiliary import read_chromsizes, check_bounds, open_comp


def adjust_coordinates(regtype, start, end, strand):
    """
    :param regtype:
    :param start:
    :param end:
    :param strand:
    :return:
    """
    norm = 1 if strand in ['+', '.'] else -1
    if regtype in ['win5p', 'win3p']:
        refloc = start if regtype == 'win5p' else end
        s, e = refloc - 250, refloc + 250
    elif regtype == 'body':
        s, e = start, end
    elif regtype == 'reg5p':
        if norm > 0:
            s, e = start - 2250, start - 250
        else:
            s, e = end + 250, end + 2250
    else:
        raise ValueError('Unknown region type: {}'.format(regtype))
    return s, e


def make_bed_roi(inputfile, outputfile, chromfile, regtype):
    """
    :param inputfile:
    :param outputfile:
    :param chromfile:
    :param regtype:
    :return:
    """
    assert os.path.isfile(inputfile), 'Invalid input file: {}'.format(inputfile)
    chrom_bounds = read_chromsizes(chromfile)
    opn, mode, dec = open_comp(inputfile, read=True)
    regions = []
    adjust = fnt.partial(adjust_coordinates, *(regtype,))
    with opn(inputfile, mode) as infile:
        header = dec(infile.readline())
        if header.startswith('#'):
            regions.append('\t'.join(header.split()))  # assert \t separation
        else:
            infile.seek(0)
        for line in infile:
            region = dec(line).split()
            start, end = adjust(int(region[1]), int(region[2]), region[5])
            try:
                _ = check_bounds(region[0], start, end, chrom_bounds, inputfile)
            except AssertionError:
                max_end = chrom_bounds[region[0]]
                assert end > max_end, 'Weirdly malformed region: {} - {} - {} - {}'.format(region, start, end, max_end)
                end = max_end
            assert start < end, 'Invalid region created: {} - {}'.format(start, end)
            region[1] = str(start)
            region[2] = str(end)
            regions.append('\t'.join(region))
    opn, mode, enc = open_comp(outputfile, read=False)
    regions.append('')  # assert newline at end of output file
    with opn(outputfile, mode) as outfile:
        _ = outfile.write(enc('\n'.join(regions)))
    return outputfile


def normalize_join(inputfiles, outputfiles):
    """
    :param inputfiles:
    :param outputfiles:
    :return:
    """
    buffer_ncbi = []
    buffer_ens = []
    prefixed_chroms = re.compile('^[0-9XYM]')
    for fp in inputfiles:
        if fp.endswith('.gz'):
            opn, mode = gz.open, 'rt'
            to_strip = '.bed.gz'
        else:
            opn, mode = open, 'r'
            to_strip = '.bed'
        with opn(fp, mode) as infile:
            for idx, line in enumerate(infile, start=1):
                if not line.strip() or line.startswith('#'):
                    continue
                cols = line.strip().split()
                region = cols[:3]
                try:
                    name = cols[3]
                except IndexError:
                    name = os.path.basename(fp).rstrip(to_strip) + '_' + str(idx)
                try:
                    _ = int(name)
                    name = os.path.basename(fp).rstrip(to_strip) + '_' + str(idx)
                except ValueError:
                    pass
                region.append(name)
                if region[0].startswith('chr'):
                    buffer_ncbi.append(region)
                    region = [region[0][3:], region[1], region[2], region[3]]
                    buffer_ens.append(region)
                else:
                    if prefixed_chroms.match(region[0]):
                        buffer_ens.append(region)
                        region = ['chr' + region[0], region[1], region[2], region[3]]
                        buffer_ncbi.append(region)
                    else:
                        buffer_ncbi.append(region)
                        buffer_ens.append(region)
    buffer_ens = sorted(buffer_ens, key=lambda x: (x[0], int(x[1]), int(x[2])))
    buffer_ncbi = sorted(buffer_ncbi, key=lambda x: (x[0], int(x[1]), int(x[2])))
    ncbi_output = outputfiles[0]
    with open(ncbi_output, 'w') as outf:
        _ = outf.write('\n'.join(['\t'.join(reg) for reg in buffer_ncbi]) + '\n')
    ensembl_output = outputfiles[1]
    with open(ensembl_output, 'w') as outf:
        _ = outf.write('\n'.join(['\t'.join(reg) for reg in buffer_ens]) + '\n')
    return ncbi_output, ensembl_output
