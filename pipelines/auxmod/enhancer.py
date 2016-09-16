# coding=utf-8

import collections as col
import io as io

from pipelines.auxmod.auxiliary import read_chromsizes, open_comp, check_bounds


def derive_encode_enh_name(names):
    """
    :param names:
    :return:
    """
    all_names = names.split('@')
    all_names = col.Counter([n.rsplit('-', 1)[0] for n in all_names])
    types = {'Proximal-Prediction': 'PP', 'Distal-Prediction': 'DP'}
    return types[all_names.most_common()[0][0]]


def process_merged_encode_enhancer(inputfile, outputfile, boundcheck, nprefix, etype):
    """
    :param inputfile:
    :param outputfile:
    :param boundcheck:
    :param nprefix:
    :param etype:
    :return:
    """
    bounds = read_chromsizes(boundcheck)
    opn, mode, conv = open_comp(inputfile, True)
    regions = []
    enhname = derive_encode_enh_name
    bchk = check_bounds
    with opn(inputfile, 'r') as infile:
        for line in infile:
            line = conv(line)
            if not line:
                continue
            c, s, e, n = line.split()
            stat = bchk(c, s, e, bounds, inputfile)
            if not stat:
                # chromosome not in check file, otherwise raises
                continue
            enhtype = enhname(n)
            if enhtype != etype:
                continue
            regions.append((c, s, e, enhtype))
    assert regions, 'No regions read from file {}'.format(inputfile)
    regions = sorted(regions, key=lambda x: (x[0], int(x[1])))
    outbuffer = io.StringIO()
    for idx, reg in enumerate(regions, start=1):
        outbuffer.write('\t'.join(reg[:3]) + '\t{}{}_{}\n'.format(nprefix, reg[3], idx))
    opn, mode, conv = open_comp(outputfile, False)
    with opn(outputfile, mode) as outf:
        _ = outf.write(conv(outbuffer.getvalue()))
    return outputfile


def process_vista_enhancer(inputfile, outputfile, boundcheck, nprefix, species):
    """
    :param inputfile:
    :param outputfile:
    :param boundcheck:
    :param nprefix:
    :param species:
    :return:
    """
    bounds = read_chromsizes(boundcheck)
    opn, mode, conv = open_comp(inputfile, True)
    regions = []
    bchk = check_bounds
    with opn(inputfile, mode) as infile:
        for line in infile:
            line = conv(line)
            if not line or not line.startswith('>'):
                continue
            header = line.strip('>').split('|')
            if header[0] != species:
                continue
            c, coord = header[1].strip().split(':')
            s, e = coord.split('-')
            stat = bchk(c, s, e, bounds, inputfile)
            if not stat:
                # chromosome not in check file, otherwise raises
                continue
            regname = nprefix + '_' + header[3].strip()[:3] + '_' + header[2].strip('elmnt ')
            regions.append((c, s, e, regname))
    assert regions, 'No regions read from file {}'.format(inputfile)
    regions = sorted(regions, key=lambda x: (x[0], int(x[1]), int(x[2])))
    outbuffer = io.StringIO()
    for reg in regions:
        outbuffer.write('\t'.join(reg) + '\n')
    opn, mode, conv = open_comp(outputfile, False)
    with opn(outputfile, mode) as outf:
        _ = outf.write(conv(outbuffer.getvalue()))
    return outputfile
