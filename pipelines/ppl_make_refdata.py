# coding=utf-8

import os as os
import fnmatch as fnm
import datetime as dt

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
                              output=collect_full_paths(csz_rawdata, '*.sizes'),
                              name='csz_init')

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

    return pipe