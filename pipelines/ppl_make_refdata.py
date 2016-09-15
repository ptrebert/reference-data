# coding=utf-8

import os as os
import fnmatch as fnm
import datetime as dt
import gzip as gz

from ruffus import *

from pipelines.auxmod.chromsizes import filter_chromosomes


def touch_checkfile(inputfiles, outputfile):
    """
    :param inputfiles:
    :param outputfile:
    :return:
    """
    timestr = dt.datetime.now().strftime('%Y-%m-%d@%A@%H:%M:%S')
    with open(outputfile, 'w') as outf:
        _ = outf.write(timestr + '\n')
    return outputfile


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
    """ unnecessary in Python 3.5+
    :param fp:
    :return:
    """
    if fp.endswith('.gz'):
        if read:
            return gz.open, 'rb', lambda x: x.decode('ascii').strip()
        else:
            return gz.open, 'wb', lambda x: x.encode('ascii')
    else:
        if read:
            return open, 'r', lambda x: x.strip()
        else:
            return open, 'w', lambda x: x


def process_std_bedfile(inputfile, outputfile, boundcheck):
    """
    :param inputfile:
    :param outputfile:
    :param boundcheck:
    :return:
    """
    bounds = read_chromsizes(boundcheck)
    opn, mode, conv = open_comp(inputfile, True)
    regions = []
    with opn(inputfile, mode) as infile:
        for line in infile:
            line = conv(line)
            if not line:
                continue
            c, s, e, n = line.split()[:4]
            try:
                assert int(e) < bounds[c], 'Region {} out of bounds ({}) in file {}'.format(line, bounds[c], inputfile)
            except KeyError:
                # chromosome not in check file, skip over region
                continue






def build_pipeline(args, config, sci_obj):
    """
    :param args:
    :param config:
    :param sci_obj:
    :return:
    """

    pipe = Pipeline(name=config.get('Pipeline', 'name'))

    # sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    # if args.gridmode:
    #     jobcall = sci_obj.ruffus_gridjob()
    # else:
    #     jobcall = sci_obj.ruffus_localjob()

    workdir = config.get('Pipeline', 'workdir')

    # ============================
    # Major task: chromosome sizes
    #
    dir_task_chromsizes = os.path.join(workdir, 'chromsizes')

    csz_rawdata = os.path.join(dir_task_chromsizes, 'rawdata')
    csz_init = pipe.originate(task_func=lambda x: x,
                              name='csz_init',
                              output=collect_full_paths(csz_rawdata, '*.sizes'))

    outdir = os.path.join(dir_task_chromsizes, 'chrom_primary')
    csz_primary = pipe.subdivide(task_func=filter_chromosomes,
                                 name='csz_primary',
                                 input=output_from(csz_init),
                                 filter=formatter('(?P<ASSM>\w+)\.chrom\.sizes'),
                                 output=[os.path.join(outdir, '{ASSM[0]}_sizes_primary.tsv'),
                                         os.path.join(outdir, '{ASSM[0]}_sizes_primary.bed')],
                                 extras=['chr[0-9A-Z]+$']).mkdir(outdir).jobs_limit(2)

    outdir = os.path.join(dir_task_chromsizes, 'chrom_auto')
    csz_auto = pipe.subdivide(task_func=filter_chromosomes,
                              name='csz_auto',
                              input=output_from(csz_init),
                              filter=formatter('(?P<ASSM>\w+)\.chrom\.sizes'),
                              output=[os.path.join(outdir, '{ASSM[0]}_sizes_autosomes.tsv'),
                                      os.path.join(outdir, '{ASSM[0]}_sizes_autosomes.bed')],
                              extras=['chr[0-9][0-9A-Z]?$']).mkdir(outdir).jobs_limit(2)

    run_task_csz = pipe.merge(task_func=touch_checkfile,
                              name='run_task_csz',
                              input=output_from(csz_primary, csz_auto),
                              output=os.path.join(dir_task_chromsizes, 'run_task_csz.chk'))

    #
    # End: major task chromosome sizes
    # ================================

    # ================================
    # Major task: enhancer
    #
    dir_task_enhancer = os.path.join(workdir, 'enhancer')

    enh_rawdata = os.path.join(dir_task_enhancer, 'rawdata')
    enh_init = pipe.originate(task_func=lambda x: x,
                              name='enh_init',
                              output=collect_full_paths(enh_rawdata, '*'))

    #
    # End: major task chromosome sizes
    # ================================

    run_all = pipe.merge(task_func=touch_checkfile,
                         name='run_all',
                         input=output_from(run_task_csz),
                         output=os.path.join(workdir, 'run_all_refdata.chk'))

    return pipe