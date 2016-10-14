# coding=utf-8

import os as os
import io as io
import datetime as dt

from ruffus import *

from pipelines.auxmod.auxiliary import read_chromsizes, open_comp, collect_full_paths,\
    collect_region_stats
from pipelines.auxmod.enhancer import process_merged_encode_enhancer, process_vista_enhancer
from pipelines.auxmod.cpgislands import process_ucsc_cgi


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
                _ = int(s)
            except ValueError:
                raise ValueError('Non-integer coordinates in BED file {} in line {}'.format(inputfile, line))
            except KeyError:
                # chromosome not in check file, skip over region
                continue
            regions.append((c, s, e, n))
    assert regions, 'No regions read from file {}'.format(inputfile)
    regions = sorted(regions, key=lambda x: (x[0], int(x[1]), int(x[2])))
    outbuffer = io.StringIO()
    for reg in regions:
        outbuffer.write('\t'.join(reg) + '\n')
    opn, mode, conv = open_comp(outputfile, False)
    with opn(outputfile, mode) as outf:
        _ = outf.write(conv(outbuffer.getvalue()))
    return outputfile


def process_noname_bedfile(inputfile, outputfile, boundcheck, nprefix):
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
    with opn(inputfile, mode) as infile:
        for line in infile:
            line = conv(line)
            if not line:
                continue
            c, s, e = line.split()[:3]
            try:
                assert int(e) < bounds[c], 'Region {} out of bounds ({}) in file {}'.format(line, bounds[c], inputfile)
                _ = int(s)
            except ValueError:
                raise ValueError('Non-integer coordinates in BED file {} in line {}'.format(inputfile, line))
            except KeyError:
                # chromosome not in check file, skip over region
                continue
            regions.append((c, s, e))
    assert regions, 'No regions read from file {}'.format(inputfile)
    regions = sorted(regions, key=lambda x: (x[0], int(x[1]), int(x[2])))
    outbuffer = io.StringIO()
    for idx, reg in enumerate(regions, start=1):
        outbuffer.write('\t'.join(reg) + '\t' + nprefix + str(idx) + '\n')
    opn, mode, conv = open_comp(outputfile, False)
    with opn(outputfile, mode) as outf:
        _ = outf.write(conv(outbuffer.getvalue()))
    return outputfile


def build_pipeline(args, config, sci_obj):
    """
    :param args:
    :param config:
    :param sci_obj:
    :return:
    """

    pipe = Pipeline(name=config.get('Pipeline', 'name'))

    workdir = config.get('Pipeline', 'workdir')

    # ==========================
    # Major task: handle genomes
    #

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    dir_task_genomes = os.path.join(workdir, 'genomes')

    dir_gen_rawdata = os.path.join(dir_task_genomes, 'rawdata')
    genomes_raw_init = pipe.originate(task_func=lambda x: x,
                                      name='genomes_raw_init',
                                      output=collect_full_paths(dir_gen_rawdata, '*'))

    dir_gen_2bit = os.path.join(dir_task_genomes, 'wg_2bit')
    cmd = config.get('Pipeline', 'to2bit')
    gen2bit = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='gen2bit',
                             input=output_from(genomes_raw_init),
                             filter=formatter('(?P<ASSEMBLY>[a-zA-Z0-9]+)\.fa\.gz$'),
                             output=os.path.join(dir_gen_2bit, '{ASSEMBLY[0]}.2bit'),
                             extras=[cmd, jobcall])

    run_task_genomes = pipe.merge(task_func=touch_checkfile,
                                  name='task_genomes',
                                  input=output_from(genomes_raw_init, gen2bit),
                                  output=os.path.join(dir_task_genomes, 'run_task_genomes.chk'))
    #
    # End: major task chromosome sizes
    # ================================

    # ============================
    # Major task: chromosome sizes
    #
    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    dir_task_chromsizes = os.path.join(workdir, 'chromsizes')

    gen2bit_files = pipe.originate(task_func=lambda x: x,
                                   name='gen2bit_files',
                                   output=collect_full_paths(dir_gen_2bit, '*.2bit'))

    dir_out_chromwg = os.path.join(dir_task_chromsizes, 'chrom_wg')
    cmd = config.get('Pipeline', 'chromwg')
    csz_chromwg = pipe.subdivide(task_func=sci_obj.get_jobf('in_outpair'),
                                 name='csz_chromwg',
                                 input=output_from(gen2bit_files),
                                 filter=formatter('(?P<ASSM>\w+)\.2bit'),
                                 output=[os.path.join(dir_out_chromwg, '{ASSM[0]}_chrom_wg.tsv'),
                                         os.path.join(dir_out_chromwg, '{ASSM[0]}_chrom_wg.bed')],
                                 extras=[cmd, jobcall]).mkdir(dir_out_chromwg)

    dir_out_chromauto = os.path.join(dir_task_chromsizes, 'chrom_auto')
    cmd = config.get('Pipeline', 'chromauto')
    csz_chromauto = pipe.subdivide(task_func=sci_obj.get_jobf('in_outpair'),
                                   name='csz_chromauto',
                                   input=output_from(gen2bit_files),
                                   filter=formatter('(?P<ASSM>\w+)\.2bit'),
                                   output=[os.path.join(dir_out_chromauto, '{ASSM[0]}_chrom_auto.tsv'),
                                           os.path.join(dir_out_chromauto, '{ASSM[0]}_chrom_auto.bed')],
                                   extras=[cmd, jobcall]).mkdir(dir_out_chromauto)

    run_task_csz = pipe.merge(task_func=touch_checkfile,
                              name='task_csz',
                              input=output_from(csz_chromwg, csz_chromauto),
                              output=os.path.join(dir_task_chromsizes, 'run_task_csz.chk'))

    #
    # End: major task chromosome sizes
    # ================================

    # ================================
    # Major task: enhancer
    #
    dir_task_enhancer = os.path.join(workdir, 'enhancer')
    temp_task_enhancer = os.path.join(dir_task_enhancer, 'temp')

    enh_rawdata = os.path.join(dir_task_enhancer, 'rawdata')
    enh_init = pipe.originate(task_func=lambda x: x,
                              name='enh_init',
                              output=collect_full_paths(enh_rawdata, '*'))

    outdir = os.path.join(dir_task_enhancer, 'bed_format')
    enh_super = pipe.transform(task_func=process_std_bedfile,
                               name='enh_super',
                               input=output_from(enh_init),
                               filter=formatter('all_(?P<ASSM>\w+)_bed\.bed'),
                               output=os.path.join(outdir, '{ASSM[0]}_enh_super.bed'),
                               extras=[os.path.join(dir_task_chromsizes,
                                                    'chrom_complete',
                                                    '{ASSM[0]}_sizes_chroms.tsv')]).mkdir(outdir)

    enh_hsa_fantom = pipe.transform(task_func=process_noname_bedfile,
                                    name='enh_hsa_fantom',
                                    input=os.path.join(enh_rawdata, 'human_permissive_enhancers_phase_1_and_2.bed.gz'),
                                    filter=formatter(),
                                    output=os.path.join(outdir, 'hg19_enh_fantom_p1p2.bed'),
                                    extras=[os.path.join(dir_task_chromsizes,
                                                         'chrom_complete',
                                                         'hg19_sizes_chroms.tsv'),
                                            'hFE_']).mkdir(outdir)

    enh_mmu_fantom = pipe.transform(task_func=process_noname_bedfile,
                                    name='enh_mmu_fantom',
                                    input=os.path.join(enh_rawdata, 'mouse_permissive_enhancers_phase_1_and_2.bed.gz'),
                                    filter=formatter(),
                                    output=os.path.join(outdir, 'mm9_enh_fantom_p1p2.bed'),
                                    extras=[os.path.join(dir_task_chromsizes,
                                                         'chrom_complete',
                                                         'mm9_sizes_chroms.tsv'),
                                            'mFE_']).mkdir(outdir)

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'catgzbed')
    enh_cat_encode = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                                  name='enh_cat_encode',
                                  input=output_from(enh_init),
                                  filter=formatter('(?P<SAMPLE>\w+)_(?P<TYPE>predictions)\.bed\.gz$'),
                                  output=os.path.join(temp_task_enhancer, 'hg19_encode_{TYPE[0]}_union.bed'),
                                  extras=[cmd, jobcall]).mkdir(temp_task_enhancer)

    cmd = config.get('Pipeline', 'mrgbed')
    enh_mrg_encode = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                    name='enh_mrg_encode',
                                    input=output_from(enh_cat_encode),
                                    filter=suffix('union.bed'),
                                    output='merged.bed',
                                    output_dir=temp_task_enhancer,
                                    extras=[cmd, jobcall])

    enh_hsa_encpp = pipe.transform(task_func=process_merged_encode_enhancer,
                                   name='enh_hsa_encpp',
                                   input=output_from(enh_mrg_encode),
                                   filter=suffix('encode_predictions_merged.bed'),
                                   output='enh_encode_prox.bed',
                                   output_dir=outdir,
                                   extras=[os.path.join(dir_task_chromsizes,
                                                        'chrom_complete',
                                                        'hg19_sizes_chroms.tsv'),
                                           'hEE', 'PP']).mkdir(outdir)

    enh_hsa_encdp = pipe.transform(task_func=process_merged_encode_enhancer,
                                   name='enh_hsa_encdp',
                                   input=output_from(enh_mrg_encode),
                                   filter=suffix('encode_predictions_merged.bed'),
                                   output='enh_encode_dist.bed',
                                   output_dir=outdir,
                                   extras=[os.path.join(dir_task_chromsizes,
                                                        'chrom_complete',
                                                        'hg19_sizes_chroms.tsv'),
                                           'hEE', 'DP']).mkdir(outdir)

    enh_hsa_vista = pipe.transform(task_func=process_vista_enhancer,
                                   name='enh_hsa_vista',
                                   input=os.path.join(enh_rawdata, 'vista_enhancers_hg19_mm9.fa'),
                                   filter=formatter(),
                                   output=os.path.join(outdir, 'hg19_enh_vista.bed'),
                                   extras=[os.path.join(dir_task_chromsizes,
                                                        'chrom_complete',
                                                        'hg19_sizes_chroms.tsv'),
                                           'hVE', 'Human']).mkdir(outdir)

    enh_mmu_vista = pipe.transform(task_func=process_vista_enhancer,
                                   name='enh_mmu_vista',
                                   input=os.path.join(enh_rawdata, 'vista_enhancers_hg19_mm9.fa'),
                                   filter=formatter(),
                                   output=os.path.join(outdir, 'mm9_enh_vista.bed'),
                                   extras=[os.path.join(dir_task_chromsizes,
                                                        'chrom_complete',
                                                        'mm9_sizes_chroms.tsv'),
                                           'mVE', 'Mouse']).mkdir(outdir)

    outdir = os.path.join(workdir, 'stats')
    enh_stats = pipe.transform(task_func=collect_region_stats,
                               name='enh_stats',
                               input=output_from(enh_mmu_vista, enh_hsa_vista,
                                                 enh_hsa_fantom, enh_mmu_fantom,
                                                 enh_hsa_encdp, enh_hsa_encpp,
                                                 enh_super),
                               filter=suffix('.bed'),
                               output='.tsv',
                               output_dir=outdir).mkdir(outdir).jobs_limit(2)

    run_task_enh = pipe.merge(task_func=touch_checkfile,
                              name='run_task_enh',
                              input=output_from(enh_super,
                                                enh_hsa_fantom, enh_mmu_fantom,
                                                enh_cat_encode, enh_mrg_encode, enh_hsa_encpp, enh_hsa_encdp,
                                                enh_hsa_vista, enh_mmu_vista,
                                                enh_stats),
                              output=os.path.join(dir_task_enhancer, 'run_task_enh.chk'))
    #
    # End: major task enhancer
    # ================================

    # ================================
    # Major task: CpG islands
    #
    dir_task_cpgislands = os.path.join(workdir, 'cpgislands')

    cgi_rawdata = os.path.join(dir_task_cpgislands, 'rawdata')
    cgi_init = pipe.originate(task_func=lambda x: x,
                              name='cgi_init',
                              output=collect_full_paths(cgi_rawdata, '*'))

    outdir = os.path.join(dir_task_cpgislands, 'bed_format')
    cgi_ucsc = pipe.transform(task_func=process_ucsc_cgi,
                              name='cgi_ucsc',
                              input=output_from(cgi_init),
                              filter=formatter('(?P<ASSM>\w+)_ucsc_cpgislands\.bed\.gz$'),
                              output=os.path.join(outdir, '{ASSM[0]}_cgi_ucsc.bed'),
                              extras=[os.path.join(dir_task_chromsizes,
                                                   'chrom_complete',
                                                   '{ASSM[0]}_sizes_chroms.tsv'),
                                      'CGI']).mkdir(outdir).jobs_limit(2)

    run_task_cgi = pipe.merge(task_func=touch_checkfile,
                              name='run_task_cgi',
                              input=output_from(cgi_ucsc),
                              output=os.path.join(dir_task_cpgislands, 'run_task_cgi.chk'))

    #
    # End: major task CpG islands
    # ================================

    run_all = pipe.merge(task_func=touch_checkfile,
                         name='run_all',
                         input=output_from(run_task_csz, run_task_enh, run_task_cgi),
                         output=os.path.join(workdir, 'run_all_refdata.chk'))

    return pipe