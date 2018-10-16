# coding=utf-8

import os as os
import csv as csv
import collections as col
import fnmatch as fnm
import io as io
import operator as op


def split_by_file(datafile, indexfile, outpath):
    """
    :param datafile:
    :param indexfile:
    :param outpath:
    :return:
    """
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'signal',
               'pv', 'qv', 'filename', 'factor', 'source', 'cell', 'sample']
    get_values = op.itemgetter(*('chrom', 'start', 'end', 'filename',
                                 'factor', 'signal', 'cell', 'sample', 'source'))
    file_buffer = col.defaultdict(list)
    with open(datafile, 'r') as dump:
        rows = csv.DictReader(dump, delimiter='\t', fieldnames=columns)
        for r in rows:
            file_buffer[r['filename']].append(get_values(r))
    index_buffer = io.StringIO()
    index_buffer.write('\t'.join(['filename', 'cellType', 'antibody', 'dataSource']) + '\n')
    for k in sorted(file_buffer.keys()):
        regions = sorted(file_buffer[k], key=lambda x: (x[0], int(x[1]), int(x[2])))
        _, _, _, _, factor, _, cell, sample, _ = regions[0]
        cell = cell.replace(' ', '-')
        new_filename = '_'.join([k, cell, factor, sample]) + '.bed'
        index_entry = '\t'.join([new_filename, cell, factor, 'ENCODE-DeepBlue-' + sample])
        index_buffer.write(index_entry + '\n')
        with open(os.path.join(outpath, new_filename), 'w') as dump:
            writer = csv.writer(dump, delimiter='\t')
            writer.writerows(regions)
    with open(indexfile, 'w') as out:
        _ = out.write(index_buffer.getvalue())
    return indexfile


def make_lola_isect_params(deepblue_check_files, base_out, fisect, fcopy, jobcall):
    """
    :param deepblue_check_files:
    :param base_out:
    :param fisect:
    :param fcopy:
    :param jobcall:
    :return:
    """
    bedfiles = [fp.replace('.dl.chk', '.bed.gz') for fp in deepblue_check_files]
    done = set()
    arglist = []
    processed = 0
    for path in bedfiles:
        if not os.path.isfile(path):
            processed += 1
            continue
        fp, fn = os.path.split(path)
        subp = fp.split('/')
        basename = fn.split('.')[0]
        if basename in done:
            continue
        replicates = fnm.filter(bedfiles, '*/' + basename + '*')
        outname = basename + '.bed'
        outdir = os.path.join(base_out, subp[-2], subp[-1], 'regions')
        if not os.path.isdir(outdir):
            os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, outname)
        if len(replicates) == 1:
            tmp = str(fcopy)
            arglist.append([[path], outpath, tmp, jobcall])
            done.add(basename)
        elif len(replicates) > 1:
            num = len(replicates)
            tmp = fisect.format(**{'NUM': num})
            arglist.append([replicates, outpath, tmp, jobcall])
            done.add(basename)
        else:
            raise AssertionError('No files identified for base name: {}'.format(basename))
    if deepblue_check_files and processed < len(deepblue_check_files):
        assert arglist, 'No argument list created'
    return arglist


def make_lola_index(inputfiles, outputfile):
    """
    :param inputfiles:
    :param outputfile:
    :return:
    """
    indexfile = os.path.join(os.path.dirname(outputfile), 'index.txt')
    if 'DNaseI' in inputfiles[0]:
        idx_header = 'filename\tcellType\ttissue\tsample\tdataSource'
    else:
        idx_header = 'filename\tcellType\ttissue\tantibody\tsample\tdataSource'
    with open(indexfile, 'w') as index:
        _ = index.write(idx_header + '\n')
        for inp in inputfiles:
            fn = os.path.basename(inp)
            proj, sample, assm, tissue, cell, assay = fn.split('.')[0].split('_')
            if assay.startswith('DNase'):
                _ = index.write('\t'.join([fn, cell, tissue, sample, 'DeepBlue']) + '\n')
            else:
                _ = index.write('\t'.join([fn, cell, tissue, assay, sample, 'DeepBlue']) + '\n')
    with open(outputfile, 'w') as check:
        _ = check.write('\n'.join(inputfiles))
    return outputfile
