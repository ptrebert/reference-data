# coding=utf-8

import csv as csv

from pipelines.auxmod.auxiliary import read_chromsizes, open_comp, check_bounds


def process_ucsc_cgi(inputfile, outputfile, boundcheck, nprefix):
    """
    :param inputfile:
    :param outputfile:
    :param boundcheck:
    :param nprefix:
    :return:
    """
    bounds = read_chromsizes(boundcheck)
    opn, mode, conv = open_comp(inputfile, True)
    regions = []
    bchk = check_bounds
    with opn(inputfile, mode) as infile:
        for line in infile:
            line = conv(line)
            if not line or line.startswith('#'):
                continue
            c, s, e = line.split()[1:4]
            stat = bchk(c, s, e, bounds, inputfile)
            if not stat:
                # chromosome not in check file, otherwise raises
                continue
            regions.append(line.strip().split('\t')[1:])
    assert regions, 'No regions read from file {}'.format(inputfile)
    regions = sorted(regions, key=lambda x: (x[0], int(x[1]), int(x[2])))
    rowbuffer = []
    for idx, reg in enumerate(regions, start=1):
        reg[3] = '{}_{}'.format(nprefix, idx)
        rowbuffer.append(reg)
    opn, mode, conv = open_comp(outputfile, False)
    cgi_header = ['#chrom', 'start', 'end', 'name', 'length',
                  'cpgNum', 'gcNum', 'perCpg', 'perGc', 'obsExp']
    with opn(outputfile, mode) as outf:
        writer = csv.writer(outf, delimiter='\t')
        writer.writerow(cgi_header)
        for row in rowbuffer:
            writer.writerow(row)
    return outputfile
