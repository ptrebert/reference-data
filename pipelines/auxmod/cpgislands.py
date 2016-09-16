# coding=utf-8

import io as io

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
            regions.append((c, s, e))
    assert regions, 'No regions read from file {}'.format(inputfile)
    regions = sorted(regions, key=lambda x: (x[0], int(x[1]), int(x[2])))
    outbuffer = io.StringIO()
    for idx, reg in enumerate(regions, start=1):
        outbuffer.write('\t'.join(reg) + '\t{}_{}\n'.format(nprefix, idx))
    opn, mode, conv = open_comp(outputfile, False)
    with opn(outputfile, mode) as outf:
        _ = outf.write(conv(outbuffer.getvalue()))
    return outputfile
