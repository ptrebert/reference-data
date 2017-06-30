# coding=utf-8

import os as os
import io as io
import datetime as dt

from ruffus import *

from pipelines.auxmod.auxiliary import read_chromsizes, open_comp, collect_full_paths
from pipelines.auxmod.enhancer import process_merged_encode_enhancer, \
    process_vista_enhancer, merge_genehancer_annotation, process_hsa_mapped, \
    build_genehancer_map_params
from pipelines.auxmod.cpgislands import process_ucsc_cgi
from pipelines.auxmod.chainfiles import build_chain_filter_commands, \
    build_symm_filter_commands, check_lifted_blocks, filter_rbest_net
from pipelines.auxmod.bedroi import make_bed_roi, normalize_join
from pipelines.auxmod.orthologs import match_ortholog_files
from pipelines.auxmod.promoter import split_promoter_files, score_promoter_regions
from pipelines.auxmod.loladb import split_by_file, make_lola_isect_params

from pipelines.auxmod.project_sarvesh import make_5p_window


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

    cmd = config.get('Pipeline', 'genstdfagz')
    dir_gen_stdfagz = os.path.join(dir_task_genomes, 'rawnorm')
    genstdfagz = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='genstdfagz',
                                input=output_from(genomes_raw_init),
                                filter=formatter('(?P<ASSM>\w+[0-9])\.fa\.gz'),
                                output=os.path.join(dir_gen_stdfagz, '{ASSM[0]}.fa.gz'),
                                extras=[cmd, jobcall]).mkdir(dir_gen_stdfagz)

    dir_gen_2bit = os.path.join(dir_task_genomes, 'wg_2bit')

    cmd = config.get('Pipeline', 'to2bit')
    gen2bit = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='gen2bit',
                             input=output_from(genstdfagz),
                             filter=formatter('(?P<ASSM>[a-zA-Z0-9]+)\.fa\.gz'),
                             output=os.path.join(dir_gen_2bit, '{ASSM[0]}.2bit'),
                             extras=[cmd, jobcall]).mkdir(dir_gen_2bit)

    gen2bitren = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='gen2bitren',
                                input=output_from(genomes_raw_init),
                                filter=formatter('Home_sapiens.+\.fa\.gz'),
                                output=os.path.join(dir_gen_2bit, 'GRCh37.2bit'),
                                extras=[cmd, jobcall])

    genomes_2bit_init = pipe.originate(task_func=lambda x: x,
                                       name='genomes_2bit_init',
                                       output=collect_full_paths(dir_gen_2bit, '*.2bit')).follows(gen2bit).follows(gen2bitren)

    dir_gen_fagz = os.path.join(dir_task_genomes, 'wg_fagz')
    cmd = config.get('Pipeline', 'tofagz')
    genfagz = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='tofagz',
                             input=output_from(genomes_2bit_init),
                             filter=suffix('.2bit'),
                             output='.fa.gz',
                             output_dir=dir_gen_fagz,
                             extras=[cmd, jobcall]).mkdir(dir_gen_fagz)

    dir_gen_fabgz = os.path.join(dir_task_genomes, 'wg_fabgz')
    cmd = config.get('Pipeline', 'tofabgz').replace('\n', ' ')
    genfabgz = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='tofabgz',
                              input=output_from(genomes_2bit_init),
                              filter=formatter('(?P<ASSM>[0-9A-Za-z]+)\.2bit'),
                              output=os.path.join(dir_gen_fabgz, '{ASSM[0]}.fa.bgz'),
                              extras=[cmd, jobcall]).mkdir(dir_gen_fabgz)

    dir_gen_fasta = os.path.join(dir_task_genomes, 'wg_fasta')
    cmd = config.get('Pipeline', 'tofasta')
    genfasta = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='genfasta',
                              input=output_from(genomes_2bit_init),
                              filter=suffix('.2bit'),
                              output='.fa',
                              output_dir=dir_gen_fasta,
                              extras=[cmd, jobcall]).mkdir(dir_gen_fasta)

    run_task_genomes = pipe.merge(task_func=touch_checkfile,
                                  name='task_genomes',
                                  input=output_from(genomes_raw_init, genstdfagz, gen2bit,
                                                    genomes_2bit_init, genfagz, genfabgz,
                                                    gen2bitren, genfasta),
                                  output=os.path.join(dir_task_genomes, 'run_task_genomes.chk'))
    #
    # End: major task genomes
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

    dir_out_chromaugo = os.path.join(dir_task_chromsizes, 'chrom_augo')
    cmd = config.get('Pipeline', 'chromaugo')
    csz_chromaugo = pipe.subdivide(task_func=sci_obj.get_jobf('in_outpair'),
                                   name='csz_chromaugo',
                                   input=output_from(gen2bit_files),
                                   filter=formatter('(?P<ASSM>\w+)\.2bit'),
                                   output=[os.path.join(dir_out_chromaugo, '{ASSM[0]}_chrom_augo.tsv'),
                                           os.path.join(dir_out_chromaugo, '{ASSM[0]}_chrom_augo.bed')],
                                   extras=[cmd, jobcall]).mkdir(dir_out_chromaugo)

    run_task_csz = pipe.merge(task_func=touch_checkfile,
                              name='task_csz',
                              input=output_from(csz_chromwg, csz_chromauto, csz_chromaugo),
                              output=os.path.join(dir_task_chromsizes, 'run_task_csz.chk'))

    #
    # End: major task chromosome sizes
    # ================================

    # ================================
    # Major task: bowtie2 indices
    #

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()
    dir_task_srmidx = os.path.join(workdir, 'srmidx')

    dir_bowtie2_base = os.path.join(dir_task_srmidx, 'bowtie2')
    cmd = config.get('Pipeline', 'bt2idx')
    bt2_idx = pipe.subdivide(task_func=sci_obj.get_jobf('in_pat'),
                             name='bt2_idx',
                             input=output_from(genfasta),
                             filter=formatter('(?P<ASSM>mm9)\.fa'),
                             output=os.path.join(dir_bowtie2_base, '{ASSM[0]}*'),
                             extras=[dir_bowtie2_base, '{ASSM[0]}*', cmd, jobcall])
    bt2_idx = bt2_idx.mkdir(dir_bowtie2_base)

    run_task_srmidx = pipe.merge(task_func=touch_checkfile,
                                 name='task_srmidx',
                                 input=output_from(bt2_idx),
                                 output=os.path.join(dir_task_srmidx, 'run_task_srmidx.chk'))

    #
    # End of: bowtie2 indices
    # ================================

    # ================================
    # Major task: CpG islands
    #
    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()
    dir_task_cpgislands = os.path.join(workdir, 'cpgislands')

    cgi_rawdata = os.path.join(dir_task_cpgislands, 'rawdata')
    cgi_init = pipe.originate(task_func=lambda x: x,
                              name='cgi_init',
                              output=collect_full_paths(cgi_rawdata, '*'))

    outdir = os.path.join(dir_task_cpgislands, 'bed_format')
    cgi_ucsc_wg = pipe.transform(task_func=process_ucsc_cgi,
                                 name='cgi_ucsc_wg',
                                 input=output_from(cgi_init),
                                 filter=formatter('(?P<ASSM>\w+)_ucsc_cpgislands\.bed\.gz$'),
                                 output=os.path.join(outdir, '{ASSM[0]}_cgi_ucsc_wg.bed'),
                                 extras=[os.path.join(dir_task_chromsizes,
                                                      'chrom_wg',
                                                      '{ASSM[0]}_chrom_wg.tsv'),
                                         'CGI']).mkdir(outdir).jobs_limit(2)

    outdir = os.path.join(dir_task_cpgislands, 'bed_format')
    cgi_ucsc_auto = pipe.transform(task_func=process_ucsc_cgi,
                                   name='cgi_ucsc_auto',
                                   input=output_from(cgi_init),
                                   filter=formatter('(?P<ASSM>\w+)_ucsc_cpgislands\.bed\.gz$'),
                                   output=os.path.join(outdir, '{ASSM[0]}_cgi_ucsc_auto.bed'),
                                   extras=[os.path.join(dir_task_chromsizes,
                                                        'chrom_auto',
                                                        '{ASSM[0]}_chrom_auto.tsv'),
                                           'CGI']).mkdir(outdir).jobs_limit(2)

    cmd = config.get('Pipeline', 'hg19toh37hd')
    cgi_ucsc_h37 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                  name='cgi_ucsc_h37',
                                  input=os.path.join(outdir, 'hg19_cgi_ucsc_wg.bed'),
                                  filter=formatter(),
                                  output=os.path.join(outdir, 'h37_cgi_ucsc_augo.bed'),
                                  extras=[cmd, jobcall])

    run_task_cgi = pipe.merge(task_func=touch_checkfile,
                              name='task_cgi',
                              input=output_from(cgi_ucsc_wg, cgi_ucsc_auto, cgi_ucsc_h37),
                              output=os.path.join(dir_task_cpgislands, 'task_cgi.chk'))

    #
    # End: major task CpG islands
    # ================================

    # ===============================
    # Major task: gene model
    #
    dir_task_genemodel = os.path.join(workdir, 'genemodel')
    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    gm_rawdata = os.path.join(dir_task_genemodel, 'rawdata')
    gm_init = pipe.originate(task_func=lambda x: x,
                             name='gm_init',
                             output=collect_full_paths(gm_rawdata, '*.gz'))

    bedout = os.path.join(dir_task_genemodel, 'bed_format')
    cmd = config.get('Pipeline', 'gtftobed')
    gm_gtftobed = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='gm_gtftobed',
                                 input=output_from(gm_init),
                                 filter=suffix('.gtf.gz'),
                                 output='.bed.gz',
                                 output_dir=bedout,
                                 extras=[cmd, jobcall]).mkdir(bedout)

    gm_gfftobed = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='gm_gfftobed',
                                 input=output_from(gm_init),
                                 filter=suffix('.gff3.gz'),
                                 output='.bed.gz',
                                 output_dir=bedout,
                                 extras=[cmd, jobcall]).mkdir(bedout)

    cmd = config.get('Pipeline', 'enstobed').replace('\n', ' ')
    gm_enstobed = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='gm_enstobed',
                                 input=output_from(gm_init),
                                 filter=formatter('ensGene_(?P<ENSVER>v[0-9]+)_UCSC_(?P<SPECIES>[a-z]+)_(?P<ASSM>\w+)\.txt\.gz'),
                                 output=os.path.join(bedout, '{SPECIES[0]}_{ASSM[0]}_ensembl_{ENSVER[0]}.bed.gz'),
                                 extras=[cmd, jobcall])

    subout = os.path.join(dir_task_genemodel, 'subsets')

    pc_subset = os.path.join(subout, 'protein_coding')
    cmd = config.get('Pipeline', 'subpchg19').replace('\n', ' ')
    gm_pc_hg19 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='gm_pc_hg19',
                                input=output_from(gm_gtftobed),
                                filter=formatter('gencode\.(?P<GENCVER>v19)\.annotation\.bed\.gz'),
                                output=os.path.join(pc_subset, 'hsa_hg19_gencode_{GENCVER[0]}.genes.bed.gz'),
                                extras=[cmd, jobcall]).mkdir(pc_subset)

    cmd = config.get('Pipeline', 'subpcmm9').replace('\n', ' ')
    gm_pc_mm9 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                               name='gm_pc_mm9',
                               input=output_from(gm_gtftobed),
                               filter=formatter('gencode\.(?P<GENCVER>vM1)\.annotation\.bed\.gz'),
                               output=os.path.join(pc_subset, 'mmu_mm9_gencode_{GENCVER[0]}.genes.bed.gz'),
                               extras=[cmd, jobcall]).mkdir(pc_subset)

    cmd = config.get('Pipeline', 'subpcbta').replace('\n', ' ')
    gm_pc_bta7 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='gm_pc_bta7',
                                input=output_from(gm_gfftobed),
                                filter=formatter('Ensembl(?P<ENSVER>75)_.+\.bed\.gz'),
                                output=os.path.join(pc_subset, 'bta_bosTau7_ensembl_v{ENSVER[0]}.genes.bed.gz'),
                                extras=[cmd, jobcall]).mkdir(pc_subset)

    cmd = config.get('Pipeline', 'subpcens').replace('\n', ' ')
    gm_pc_ens = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                               name='gm_pc_ens',
                               input=output_from(gm_enstobed),
                               filter=formatter('(?P<ANNOTID>\w+)\.bed\.gz'),
                               output=os.path.join(pc_subset, '{ANNOTID[0]}.genes.bed.gz'),
                               extras=[cmd, jobcall]).mkdir(pc_subset)

    # dump various regions of interest from the subset of protein coding genes
    pc_roi_dump = os.path.join(pc_subset, 'roi_bed')
    gm_pc_reg5p = pipe.transform(task_func=make_bed_roi,
                                 name='gm_pc_reg5p',
                                 input=output_from(gm_pc_hg19, gm_pc_mm9, gm_pc_bta7, gm_pc_ens),
                                 filter=formatter('(?P<ANNOTID>[a-z]{3}_(?P<ASSM>[a-zA-Z0-9]+)_\w+)\.genes\.bed\.gz'),
                                 output=os.path.join(pc_roi_dump, '{ANNOTID[0]}.reg5p.bed.gz'),
                                 extras=[os.path.join(dir_out_chromaugo, '{ASSM[0]}_chrom_augo.tsv'),
                                         'reg5p']).mkdir(pc_roi_dump).jobs_limit(4)

    gm_pc_uprr = pipe.transform(task_func=make_bed_roi,
                                name='gm_pc_uprr',
                                input=output_from(gm_pc_hg19, gm_pc_mm9, gm_pc_bta7, gm_pc_ens),
                                filter=formatter('(?P<ANNOTID>[a-z]{3}_(?P<ASSM>[a-zA-Z0-9]+)_\w+)\.genes\.bed\.gz'),
                                output=os.path.join(pc_roi_dump, '{ANNOTID[0]}.uprr.bed.gz'),
                                extras=[os.path.join(dir_out_chromaugo, '{ASSM[0]}_chrom_augo.tsv'),
                                        'uprr']).mkdir(pc_roi_dump).jobs_limit(4)

    gm_pc_body = pipe.transform(task_func=make_bed_roi,
                                name='gm_pc_body',
                                input=output_from(gm_pc_hg19, gm_pc_mm9, gm_pc_bta7, gm_pc_ens),
                                filter=formatter('(?P<ANNOTID>[a-z]{3}_(?P<ASSM>[a-zA-Z0-9]+)_\w+)\.genes\.bed\.gz'),
                                output=os.path.join(pc_roi_dump, '{ANNOTID[0]}.body.bed.gz'),
                                extras=[os.path.join(dir_out_chromaugo, '{ASSM[0]}_chrom_augo.tsv'),
                                        'body']).mkdir(pc_roi_dump).jobs_limit(4)

    cmd = config.get('Pipeline', 'subpchdf').replace('\n', ' ')
    pc_roi_conv_hdf = os.path.join(pc_subset, 'roi_hdf')
    gm_pc_hdf = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                               name='gm_pc_hdf',
                               input=output_from(gm_pc_reg5p, gm_pc_uprr, gm_pc_body),
                               filter=formatter('(?P<SPEC>\w+)_(?P<ASSM>\w+)_(?P<ANNOT>\w+)_'
                                                '(?P<VER>\w+)\.(?P<REGTYPE>\w+)\.bed\.gz'),
                               output=os.path.join(pc_roi_conv_hdf, '{SPEC[0]}_{ASSM[0]}_{ANNOT[0]}_'
                                                                    '{VER[0]}.{REGTYPE[0]}.h5'),
                               extras=[cmd, jobcall]).mkdir(pc_roi_conv_hdf)

    run_task_genemodel = pipe.merge(task_func=touch_checkfile,
                                    name='task_genemodel',
                                    input=output_from(gm_gtftobed, gm_gfftobed, gm_enstobed,
                                                      gm_pc_hg19, gm_pc_mm9, gm_pc_bta7, gm_pc_ens,
                                                      gm_pc_uprr, gm_pc_reg5p, gm_pc_body,
                                                      gm_pc_hdf),
                                    output=os.path.join(dir_task_genemodel, 'run_task_genemodel.chk'))

    #
    # End of major task: gene model
    # ================================

    # ================================
    # Major task: transcript model, Salmon indices
    #
    dir_task_transmodel = os.path.join(workdir, 'transcriptome')
    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    tm_init = pipe.originate(task_func=lambda x: x,
                             name='tm_init',
                             output=collect_full_paths(pc_subset, '*.transcripts.bed.gz'))

    dir_tm_fasta = os.path.join(dir_task_transmodel, 'fasta')
    cmd = config.get('Pipeline', 'tmfasta')
    tm_fasta = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='tm_fasta',
                              input=output_from(tm_init),
                              filter=formatter('(?P<SPECIES>[a-z]+)_(?P<ASSM>[a-z0-9A-Z]+)_(?P<ANNVER>\w+)\.transcripts\.bed\.gz'),
                              output=os.path.join(dir_tm_fasta, '{SPECIES[0]}_{ASSM[0]}_{ANNVER[0]}.transcripts.fa.gz'),
                              extras=[cmd, jobcall]).mkdir(dir_tm_fasta)

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    dir_tm_qindex = os.path.join(dir_task_transmodel, 'qindex')

    # k=31 - read length ~75 and above; Salmon does not support higher values for k
    # k=31 is the developer's recommendation for read length of 75
    cmd = config.get('Pipeline', 'genidx31').replace('\n', ' ')
    tm_idxgen31 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='tm_idxgen31',
                                 input=output_from(tm_fasta),
                                 filter=formatter('(?P<TRANSCRIPTOME>(hsa|mmu)_\w+)\.transcripts\.fa\.gz'),
                                 output='{subpath[0][1]}/qindex/{TRANSCRIPTOME[0]}.k31.idx/rsd.bin',
                                 extras=[cmd, jobcall])
    tm_idxgen31 = tm_idxgen31.mkdir(dir_tm_qindex)

    # k=19 - read length ~50
    cmd = config.get('Pipeline', 'genidx19').replace('\n', ' ')
    tm_idxgen19 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='tm_idxgen19',
                                 input=output_from(tm_fasta),
                                 filter=formatter('(?P<TRANSCRIPTOME>(cfa|bta|mmu|ssc|gga)_\w+)\.transcripts\.fa\.gz'),
                                 output='{subpath[0][1]}/qindex/{TRANSCRIPTOME[0]}.k19.idx/rsd.bin',
                                 extras=[cmd, jobcall])
    tm_idxgen19 = tm_idxgen19.mkdir(dir_tm_qindex)

    run_task_transmodel = pipe.merge(task_func=touch_checkfile,
                                     name='task_transmodel',
                                     input=output_from(tm_fasta, tm_idxgen31, tm_idxgen19),
                                     output=os.path.join(dir_task_transmodel, 'run_task_transmodel.chk')).mkdir(dir_task_transmodel)
    #
    # End of major task: transcript model
    # ====================================

    # ================================
    # Major task: orthologs
    #
    dir_task_orthologs = os.path.join(workdir, 'orthologs')

    ortho_rawdata = os.path.join(dir_task_orthologs, 'rawdata')
    ortho_files = collect_full_paths(ortho_rawdata, '*.gz')
    ortho_init = pipe.originate(task_func=lambda x: x,
                                name='ortho_init',
                                output=ortho_files)

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    ortho_hdf_out = os.path.join(dir_task_orthologs, 'hdf')
    cmd = config.get('Pipeline', 'proc_hcop').replace('\n', ' ')
    proc_hcop = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                             name='proc_hcop',
                             input=output_from(ortho_init),
                             filter=formatter('human_[a-z]+_(?P<SOURCE>hcop)_six_column\.txt\.gz'),
                             output=os.path.join(ortho_hdf_out, '{SOURCE[0]}_6species.h5'),
                             extras=[cmd, jobcall])
    proc_hcop = proc_hcop.mkdir(ortho_hdf_out)

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'proc_orthodb').replace('\n', ' ')
    proc_orthodb = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                                name='proc_orthodb',
                                input=output_from(ortho_init),
                                filter=formatter('(?P<SOURCE>odb9)_\w+\.tab\.gz$'),
                                output=os.path.join(ortho_hdf_out, '{SOURCE[0]}_6species.h5'),
                                extras=[cmd, jobcall])
    proc_orthodb = proc_orthodb.mkdir(ortho_hdf_out)
    proc_orthodb = proc_orthodb.follows(run_task_genemodel)

    run_task_orthologs = pipe.merge(task_func=touch_checkfile,
                                    name='task_orthologs',
                                    input=output_from(proc_hcop, proc_orthodb),
                                    output=os.path.join(dir_task_orthologs, 'run_task_orthologs.chk'))

    #
    # End of major task: orthologs
    # ================================

    # ================================
    # Major task promoter
    #
    dir_task_promoter = os.path.join(workdir, 'promoter')
    encode_expcell_map = os.path.join(workdir, 'encode_exp-to-cell.json')

    prom_rawdata = os.path.join(dir_task_promoter, 'rawdata')
    prom_init_hg19 = pipe.originate(task_func=lambda x: x,
                                    name='prom_init_hg19',
                                    output=collect_full_paths(os.path.join(prom_rawdata, '20161130-062049-promoter-like-hg19-BothDNaseAndH3K4me3.v3'),
                                                              '*.bed.gz'))

    dir_prom_temp = os.path.join(dir_task_promoter, 'temp')
    prom_prox_hg19 = pipe.merge(task_func=split_promoter_files,
                                name='prom_prox_hg19',
                                input=output_from(prom_init_hg19),
                                output=os.path.join(dir_prom_temp, 'ENCODE_hg19_proximal.bed'),
                                extras=[encode_expcell_map, 'Proximal']).mkdir(dir_prom_temp)

    prom_dist_hg19 = pipe.merge(task_func=split_promoter_files,
                                name='prom_dist_hg19',
                                input=output_from(prom_init_hg19),
                                output=os.path.join(dir_prom_temp, 'ENCODE_hg19_distal.bed'),
                                extras=[encode_expcell_map, 'Distal']).mkdir(dir_prom_temp)

    cmd = config.get('Pipeline', 'prommerge')
    prom_mrg = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='prom_mrg',
                              input=output_from(prom_prox_hg19, prom_dist_hg19),
                              filter=suffix('.bed'),
                              output='.mrg.bed',
                              output_dir=dir_prom_temp,
                              extras=[cmd, jobcall])

    dir_prom_bed = os.path.join(dir_task_promoter, 'bed_format')
    prom_bed = pipe.transform(task_func=score_promoter_regions,
                              name='prom_bed',
                              input=output_from(prom_mrg),
                              filter=suffix('.mrg.bed'),
                              output='_prom.bed',
                              output_dir=dir_prom_bed).mkdir(dir_prom_bed)


    #
    # End of major task: promoter
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

    dir_task_enh_gh = os.path.join(enh_rawdata, 'GeneHancerDump_v4-4-2')
    enh_mrg_genehancer = pipe.merge(task_func=merge_genehancer_annotation,
                                    name='merge_genehancer',
                                    input=[os.path.join(dir_task_enh_gh, 'enhancers.txt'),
                                           os.path.join(dir_task_enh_gh, 'enhancer_gene_associations.txt')],
                                    output=os.path.join(outdir, 'hg38_enh_genehancer.bed'))
    enh_mrg_genehancer = enh_mrg_genehancer.mkdir(outdir)
    enh_mrg_genehancer = enh_mrg_genehancer.follows(enh_init)

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('CrossmapEnv')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'hg19map')
    enh_gh_hg19map = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                    name='enh_gh_hg19map',
                                    input=output_from(enh_mrg_genehancer),
                                    filter=formatter(),
                                    output=os.path.join(temp_task_enhancer, 'hg19_enh_genehancer_raw.bed'),
                                    extras=[cmd, jobcall])

    enh_gh_hg19proc = pipe.transform(task_func=process_hsa_mapped,
                                     name='enh_gh_hg19proc',
                                     input=output_from(enh_gh_hg19map),
                                     filter=suffix('raw.bed'),
                                     output='mrg.bed',
                                     extras=[os.path.join(dir_task_enh_gh, 'enhancers.txt')])

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'hg19final')
    enh_gh_hg19final = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                      name='enh_gh_hg19final',
                                      input=output_from(enh_gh_hg19proc),
                                      filter=suffix('_mrg.bed'),
                                      output='.bed',
                                      output_dir=outdir,
                                      extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'hg19mldata').replace('\n', ' ')
    enh_gh_mldata = pipe.subdivide(task_func=sci_obj.get_jobf('in_pat'),
                                   name='enh_gh_mldata',
                                   input=output_from(enh_gh_hg19final),
                                   filter=formatter(),
                                   output=os.path.join(temp_task_enhancer, 'hg19_genehancer_mldata_cvidx*'),
                                   extras=[temp_task_enhancer, 'hg19*.pck', cmd, jobcall])

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'hg19tomap')
    enh_gh_tomap = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                  name='enh_gh_tomap',
                                  input=output_from(enh_gh_hg19final),
                                  filter=suffix('.bed'),
                                  output='_map.bed',
                                  output_dir=temp_task_enhancer,
                                  extras=[cmd, jobcall])

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('CrossmapEnv')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'ghmap')
    enh_gh_map_params = build_genehancer_map_params(os.path.join(temp_task_enhancer, 'hg19_enh_genehancer_map.bed'),
                                                    collect_full_paths(os.path.join(workdir, 'chainfiles', 'rbest_chain_filt'), 'hg19_*'),
                                                    temp_task_enhancer, cmd, jobcall)
    enh_gh_map = pipe.files(sci_obj.get_jobf('in_out'),
                            enh_gh_map_params,
                            name='enh_gh_map')
    enh_gh_map = enh_gh_map.mkdir(temp_task_enhancer)
    enh_gh_map = enh_gh_map.follows(enh_gh_tomap)

    run_task_enh = pipe.merge(task_func=touch_checkfile,
                              name='run_task_enh',
                              input=output_from(enh_super,
                                                enh_hsa_fantom, enh_mmu_fantom,
                                                enh_cat_encode, enh_mrg_encode, enh_hsa_encpp, enh_hsa_encdp,
                                                enh_hsa_vista, enh_mmu_vista, enh_mrg_genehancer, enh_gh_hg19map,
                                                enh_gh_hg19proc, enh_gh_hg19final, enh_gh_tomap, enh_gh_map),
                              output=os.path.join(dir_task_enhancer, 'run_task_enh.chk'))
    #
    # End: major task enhancer
    # ================================

    # ================================
    # Major task: chain files
    #
    dir_task_chainfiles = os.path.join(workdir, 'chainfiles')
    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    simon_says = config.getboolean('Pipeline', 'liftover_task')
    jacques_says = config.getboolean('Pipeline', 'liftover_check')

    chf_rawdata = os.path.join(dir_task_chainfiles, 'rawdata')
    raw_chainfiles = collect_full_paths(chf_rawdata, '*over.chain.gz')
    chf_init = pipe.originate(task_func=lambda x: x,
                              name='chf_init',
                              output=raw_chainfiles)

    prefilter_chains = os.path.join(dir_task_chainfiles, 'prefilter')
    cmd = config.get('Pipeline', 'chfprefilt').replace('\n', ' ')
    chf_prefilter = pipe.files(sci_obj.get_jobf('in_out'),
                               build_chain_filter_commands(raw_chainfiles, dir_out_chromauto,
                                                           prefilter_chains, cmd, jobcall),
                               name='chf_prefilter')
    chf_prefilter = chf_prefilter.mkdir(prefilter_chains)

    chain_re = '(?P<TARGET>\w+)_to_(?P<QUERY>\w+)(?P<EXT>\.[\w\.]+)'

    cmd = config.get('Pipeline', 'chfswap')
    chf_swap_dir = os.path.join(dir_task_chainfiles, 'tmp_swap')
    chf_swap = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='chf_swap',
                              input=output_from(chf_prefilter),
                              filter=formatter(chain_re),
                              output=os.path.join(chf_swap_dir, '{QUERY[0]}_to_{TARGET[0]}.tbest.chain.gz'),
                              extras=[cmd, jobcall])
    chf_swap = chf_swap.mkdir(chf_swap_dir)

    cmd = config.get('Pipeline', 'qrybnet')
    chf_rbest_net = os.path.join(dir_task_chainfiles, 'rbest_net')
    qrybnet = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='qrybnet',
                             input=output_from(chf_swap),
                             filter=formatter(chain_re),
                             output=os.path.join(chf_rbest_net, '{TARGET[0]}_to_{QUERY[0]}.rbest.net.gz'),
                             extras=[cmd, jobcall])
    qrybnet = qrybnet.mkdir(chf_rbest_net)

    cmd = config.get('Pipeline', 'qrybchain')
    chf_rbest_chain = os.path.join(dir_task_chainfiles, 'rbest_chain')
    qrybchain = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                               name='qrybchain',
                               input=output_from(qrybnet),
                               filter=formatter(chain_re),
                               output=os.path.join(chf_rbest_chain, '{TARGET[0]}_to_{QUERY[0]}.rbest.chain.gz'),
                               extras=[cmd, jobcall])
    qrybchain = qrybchain.mkdir(chf_rbest_chain)

    cmd = config.get('Pipeline', 'trgbchain')
    trgbchain = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                               name='trgbchain',
                               input=output_from(qrybchain),
                               filter=formatter(chain_re),
                               output=os.path.join(chf_rbest_chain, '{QUERY[0]}_to_{TARGET[0]}.rbest.chain.gz'),
                               extras=[cmd, jobcall]).mkdir(chf_rbest_chain)
    trgbchain = trgbchain.mkdir(chf_rbest_chain)

    # The netting process with subsequent conversion back to chain files
    # somehow re-introduces low scoring/small chains into the chain file
    # despite the fact that the initial filtering (step 1) removed all of those
    # Second, netToBed does not print the strand information, just the query
    # chromosome name. Thus, have to use chainfiles and chainToAxt -bed
    # to get target regions with strand information in query that can then
    # be lifted to the actual query
    rbest_re = '(?P<TARGET>(hg19|mm9))_to_(?P<QUERY>(mm9|mm10|bosTau7|canFam3|susScr2|galGal3))(?P<EXT>\.[\w\.]+)'

    cmd = config.get('Pipeline', 'trgbfilt').replace('\n', ' ')
    chf_rbest_filter = os.path.join(dir_task_chainfiles, 'rbest_chain_filt')
    trgbfilt = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='trgbfilt',
                              input=output_from(trgbchain),
                              filter=formatter(rbest_re),
                              output=os.path.join(chf_rbest_filter, '{TARGET[0]}_to_{QUERY[0]}.rbest.chain.filt.gz'),
                              extras=[cmd, jobcall])
    trgbfilt = trgbfilt.mkdir(chf_rbest_filter)

    cmd = config.get('Pipeline', 'trgbbed').replace('\n', ' ')
    chf_rbest_bed = os.path.join(dir_task_chainfiles, 'rbest_chain_bed')
    trgbbed = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='trgbbed',
                             input=output_from(trgbfilt),
                             filter=formatter(rbest_re),
                             output=os.path.join(chf_rbest_bed, '{TARGET[0]}_to_{QUERY[0]}.rbest.bed'),
                             extras=[cmd, jobcall])
    trgbbed = trgbbed.mkdir(chf_rbest_bed)

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('CrossmapEnv')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'crossmap_time').replace('\n', ' ')
    chf_qry_crossmap = os.path.join(dir_task_chainfiles, 'qry_mapped', 'crossmap')
    crossmap = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='crossmap',
                              input=output_from(trgbbed),
                              filter=formatter(rbest_re),
                              output=os.path.join(chf_qry_crossmap, '{TARGET[0]}_to_{QUERY[0]}.mapped.bed'),
                              extras=[cmd, jobcall])
    crossmap = crossmap.mkdir(chf_qry_crossmap)

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'liftover_time').replace('\n', ' ')
    chf_qry_liftover = os.path.join(dir_task_chainfiles, 'qry_mapped', 'liftover')
    liftover = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='liftover',
                              input=output_from(trgbbed),
                              filter=formatter(rbest_re),
                              output=os.path.join(chf_qry_liftover, '{TARGET[0]}_to_{QUERY[0]}.lifted.bed'),
                              extras=[cmd, jobcall])
    liftover = liftover.mkdir(chf_qry_liftover)
    liftover = liftover.active_if(simon_says)

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    dir_chf_ovlchecks = os.path.join(dir_task_chainfiles, 'ovl_checks')
    cmd = config.get('Pipeline', 'bedselfovl')
    trgselfovl = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='trgselfovl',
                                input=output_from(trgbbed),
                                filter=suffix('.rbest.bed'),
                                output='.trg.self-ovl.txt',
                                output_dir=dir_chf_ovlchecks,
                                extras=[cmd, jobcall])
    trgselfovl = trgselfovl.mkdir(dir_chf_ovlchecks)
    trgselfovl = trgselfovl.active_if(jacques_says)

    cmd = config.get('Pipeline', 'bedselfovl')
    cmqselfovl = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='cmqselfovl',
                                input=output_from(crossmap),
                                filter=suffix('.mapped.bed'),
                                output='.cmq.self-ovl.txt',
                                output_dir=dir_chf_ovlchecks,
                                extras=[cmd, jobcall])
    cmqselfovl = cmqselfovl.mkdir(dir_chf_ovlchecks)
    cmqselfovl = cmqselfovl.active_if(jacques_says)

    lfqselfovl = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='lfqselfovl',
                                input=output_from(liftover),
                                filter=suffix('.lifted.bed'),
                                output='.lfq.self-ovl.txt',
                                output_dir=dir_chf_ovlchecks,
                                extras=[cmd, jobcall])
    lfqselfovl = lfqselfovl.mkdir(dir_chf_ovlchecks)
    lfqselfovl = lfqselfovl.active_if(jacques_says)

    chf_mrgblocks = os.path.join(dir_task_chainfiles, 'tsv_map')
    cmd = config.get('Pipeline', 'mrgblocks').replace('\n', ' ')
    merge_blocks = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                  name='merge_blocks',
                                  input=output_from(trgbbed),
                                  filter=formatter(rbest_re),
                                  output=os.path.join(chf_mrgblocks, '{TARGET[0]}_to_{QUERY[0]}.aln.tsv.gz'),
                                  extras=[cmd, jobcall])
    merge_blocks = merge_blocks.mkdir(chf_mrgblocks)
    merge_blocks = merge_blocks.follows(crossmap)

    cmd = config.get('Pipeline', 'switchblocks').replace('\n', ' ')
    switch_blocks = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                   name='switch_blocks',
                                   input=output_from(trgbbed),
                                   filter=formatter('(?P<TARGET>hg19)_to_(?P<QUERY>mm9)(?P<EXT>\.[\w\.]+)'),
                                   output=os.path.join(chf_mrgblocks, '{QUERY[0]}_to_{TARGET[0]}.aln.tsv.gz'),
                                   extras=[cmd, jobcall])
    switch_blocks = switch_blocks.follows(merge_blocks)

    alntsv_re = '(?P<TARGET>\w+)_to_(?P<QUERY>\w+)(?P<EXT>\.[\w\.]+)'

    cmd = config.get('Pipeline', 'mapidx').replace('\n', ' ')
    hdfmap_dir = os.path.join(dir_task_chainfiles, 'hdf_map')
    mapidx = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='mapidx',
                            input=output_from(merge_blocks, switch_blocks),
                            filter=formatter(alntsv_re),
                            output=os.path.join(hdfmap_dir, '{TARGET[0]}_to_{QUERY[0]}.idx.h5'),
                            extras=[cmd, jobcall])
    mapidx = mapidx.mkdir(hdfmap_dir)

    run_task_chains = pipe.merge(task_func=touch_checkfile,
                                 name='task_chains',
                                 input=output_from(chf_init, chf_prefilter, chf_swap,
                                                   qrybnet, qrybchain, trgbchain, trgbfilt,
                                                   trgbbed, crossmap, liftover,
                                                   trgselfovl, cmqselfovl, lfqselfovl,
                                                   merge_blocks, switch_blocks, mapidx),
                                 output=os.path.join(dir_task_chainfiles, 'run_task_chainfiles.chk'))

    #
    # End of major task: chain files
    # ================================

    # ================================
    # Major task: blacklist regions for NGS data
    #
    dir_task_blacklists = os.path.join(workdir, 'blacklists')
    blacklists_rawdata = os.path.join(dir_task_blacklists, 'rawdata')

    dir_blnormjoin = os.path.join(dir_task_blacklists, 'normjoin')
    blnj_grch37 = pipe.collate(task_func=normalize_join,
                               name='blnjgrch37',
                               input=[os.path.join(blacklists_rawdata, fn) for fn in ['hs37d5_extENCODE_blacklist_DEEP.bed',
                                                                                      'wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz']],
                               filter=formatter(),
                               output=[os.path.join(dir_blnormjoin, 'hg19_ngs_blacklist.bed'),
                                       os.path.join(dir_blnormjoin, 'GRCh37_ngs_blacklist.bed')]).mkdir(dir_blnormjoin)

    blnj_grch38 = pipe.collate(task_func=normalize_join,
                               name='blnjgrch38',
                               input=[os.path.join(blacklists_rawdata, fn) for fn in ['hg38.blacklist.bed.gz']],
                               filter=formatter(),
                               output=[os.path.join(dir_blnormjoin, 'hg38_ngs_blacklist.bed'),
                                       os.path.join(dir_blnormjoin, 'GRCh38_ngs_blacklist.bed')]).mkdir(dir_blnormjoin)

    blnj_grcm38 = pipe.collate(task_func=normalize_join,
                               name='blnjgrcm38',
                               input=[os.path.join(blacklists_rawdata, fn) for fn in ['GRCm38_ENCODE-DKFZ_blacklist_DEEP.bed',
                                                                                      'GRCm38_alt-SHeyne_blacklist_DEEP.bed',
                                                                                      'mm10.blacklist.bed.gz']],
                               filter=formatter(),
                               output=[os.path.join(dir_blnormjoin, 'mm10_ngs_blacklist.bed'),
                                       os.path.join(dir_blnormjoin, 'GRCm38_ngs_blacklist.bed')]).mkdir(dir_blnormjoin)

    blnj_mgsc37 = pipe.collate(task_func=normalize_join,
                               name='blnjmgsc37',
                               input=[os.path.join(blacklists_rawdata, fn) for fn in ['mm9-blacklist.bed.gz']],
                               filter=formatter(),
                               output=[os.path.join(dir_blnormjoin, 'mm9_ngs_blacklist.bed'),
                                       os.path.join(dir_blnormjoin, 'MGSCv37_ngs_blacklist.bed')]).mkdir(dir_blnormjoin)

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    blacklists_init = pipe.originate(task_func=lambda x: x,
                                     name='blacklists_init',
                                     output=collect_full_paths(dir_blnormjoin, '*.bed'))

    dir_blmerge = os.path.join(dir_task_blacklists, 'merged')
    cmd = config.get('Pipeline', 'blmerge')
    blmerge = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='blmerge',
                             input=output_from(blacklists_init),
                             filter=suffix('.bed'),
                             output='.mrg.bed',
                             output_dir=dir_blmerge,
                             extras=[cmd, jobcall]).mkdir(dir_blmerge)

    run_task_blacklists = pipe.merge(task_func=touch_checkfile,
                                     name='task_blacklists',
                                     input=output_from(blacklists_init, blnj_grch37, blnj_mgsc37,
                                                       blnj_grcm38, blnj_grch38, blmerge),
                                     output=os.path.join(dir_task_blacklists, 'run_task_blacklists.chk'))

    #
    # End of major task: chain files
    # ================================

    # ================================
    # Major task: conservation tracks
    #
    dir_task_conservation = os.path.join(workdir, 'conservation')
    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()
    cons_rawdata = os.path.join(dir_task_conservation, 'rawdata')
    raw_consfiles = collect_full_paths(cons_rawdata, '*.gz')
    cons_init = pipe.originate(task_func=lambda x: x,
                               name='cons_init',
                               output=raw_consfiles)

    cmd = config.get('Pipeline', 'elembed').replace('\n', ' ')
    dir_elemtobed = os.path.join(dir_task_conservation, 'phastcons', 'elements')
    cons_elemtobed = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                    name='elemtobed',
                                    input=output_from(cons_init),
                                    filter=suffix('.tsv.gz'),
                                    output='.bed.gz',
                                    output_dir=dir_elemtobed,
                                    extras=[cmd, jobcall]).mkdir(dir_elemtobed)

    sci_obj.set_config_env(dict(config.items('NodeJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    dir_phylop_scores = os.path.join(dir_task_conservation, 'phylop', 'scores')
    cmd = config.get('Pipeline', 'phylopscores')
    phylop_scores = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                                 name='phylop_scores',
                                 input=output_from(cons_init),
                                 filter=formatter('chr[0-9XY]+\.(?P<SOURCE>phyloP46way)\.wigFix\.gz'),
                                 output=os.path.join(dir_phylop_scores, 'hg19_{SOURCE[0]}_scores.h5'),
                                 extras=[cmd, jobcall])
    phylop_scores = phylop_scores.mkdir(dir_phylop_scores)

    run_task_conservation = pipe.merge(task_func=touch_checkfile,
                                       name='task_cons',
                                       input=output_from(cons_init, cons_elemtobed,
                                                         phylop_scores),
                                       output=os.path.join(dir_task_conservation, 'run_task_conservation.chk'))

    #
    # End of major task: conservation tracks
    # ======================================

    # ======================================
    # Major task: TSS modeling
    #
    dir_task_tssmodel = os.path.join(workdir, 'tss')
    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    tss_cage_rawfiles = collect_full_paths(os.path.join(dir_task_tssmodel, 'rawdata'), '*.bed')
    tss_init = pipe.originate(task_func=lambda x: x,
                              name='tss_init',
                              output=tss_cage_rawfiles)

    dir_tss_isect = os.path.join(dir_task_tssmodel, 'isectfeat')
    cmd = config.get('Pipeline', 'tssisect')
    tss_isect = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                               name='tss_isect',
                               input=output_from(tss_init),
                               filter=suffix('.bed'),
                               output='_closes_feat.tsv',
                               output_dir=dir_tss_isect,
                               extras=[cmd, jobcall]).mkdir(dir_tss_isect)

    #
    # End of major task: TSS
    # =======================================

    # =======================================
    # Major task: LOLA database
    #

    dir_task_loladb = os.path.join(workdir, 'lola')
    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    deepblue_check_files = collect_full_paths(os.path.join(dir_task_loladb, 'rawdata', 'download'), '*dl.chk')

    lola_init = pipe.originate(task_func=lambda x: x,
                               output=deepblue_check_files,
                               name='lola_init')

    cmd_isect = config.get('Pipeline', 'lola_isect')
    cmd_copy = config.get('Pipeline', 'lola_copy')
    lola_isect_params = make_lola_isect_params(deepblue_check_files,
                                               os.path.join(dir_task_loladb, 'LOLACustom'),
                                               cmd_isect, cmd_copy, jobcall)
    lola_isect = pipe.files(sci_obj.get_jobf('ins_out'),
                            lola_isect_params,
                            name='lola_isect')
    lola_isect = lola_isect.mkdir(os.path.join(dir_task_loladb, 'LOLACustom', 'hg19'))
    lola_isect = lola_isect.mkdir(os.path.join(dir_task_loladb, 'LOLACustom', 'mm9'))

    #
    # End of major task: LOLA DB
    # =======================================


    run_all = pipe.merge(task_func=touch_checkfile,
                         name='run_all',
                         input=output_from(run_task_csz, run_task_enh, run_task_cgi,
                                           run_task_genemodel, run_task_transmodel),
                                           #run_task_chains),
                         output=os.path.join(workdir, 'run_all_refdata.chk'))

    # =========================================================================
    # =========================================================================
    # Extra section to create some reference annotation for particular projects
    # =========================================================================
    # =========================================================================

    # Sarvesh: targeted vs dispersed transcription initiation
    #
    #
    dir_project_sarvesh = os.path.join(workdir, 'projects', 'sarvesh')
    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    proj_srv_init = pipe.originate(task_func=lambda x: x,
                                   name='proj_srv_init',
                                   output=collect_full_paths(gm_rawdata, '*v19*.gtf.gz'))

    dir_pcgenes = os.path.join(dir_project_sarvesh, 'pcgenes')
    cmd = config.get('Pipeline', 'srv_pcgenes')
    srv_pcgenes = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='srv_pcgenes',
                                 input=output_from(proj_srv_init),
                                 filter=suffix('.gtf.gz'),
                                 output='.bed.gz',
                                 output_dir=dir_pcgenes,
                                 extras=[cmd, jobcall]).mkdir(dir_pcgenes)

    srv_2kbwin = pipe.transform(task_func=make_5p_window,
                                name='srv_2kbwin',
                                input=output_from(srv_pcgenes),
                                filter=suffix('.bed.gz'),
                                output='_5p_2kb.bed.gz',
                                output_dir=dir_pcgenes,
                                extras=[1500, 500])

    cmd = config.get('Pipeline', 'srv_cgiprom').replace('\n', ' ')
    srv_cgiprom = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='srv_cgiprom',
                                 input=output_from(srv_2kbwin),
                                 filter=suffix('.bed.gz'),
                                 output='_cgi.bed',
                                 output_dir=dir_pcgenes,
                                 extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'srv_noncgi').replace('\n', ' ')
    srv_noncgi = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='srv_noncgi',
                                input=output_from(srv_2kbwin),
                                filter=suffix('.bed.gz'),
                                output='_noncgi.bed',
                                output_dir=dir_pcgenes,
                                extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'srv_annreg')
    srv_annreg = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='srv_annreg',
                                input=output_from(srv_cgiprom, srv_noncgi),
                                filter=suffix('.bed'),
                                output='_feat.bed',
                                output_dir=dir_pcgenes,
                                extras=[cmd, jobcall])

    # different gene set: all protein coding

    dir_allgenes = os.path.join(dir_project_sarvesh, 'allgenes')
    cmd = config.get('Pipeline', 'srv_allgenes')
    srv_allgenes = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                  name='srv_allgenes',
                                  input=output_from(proj_srv_init),
                                  filter=suffix('.gtf.gz'),
                                  output='.bed.gz',
                                  output_dir=dir_allgenes,
                                  extras=[cmd, jobcall]).mkdir(dir_allgenes)

    cmd = config.get('Pipeline', 'srv_extgenes')
    srv_extgenes = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                  name='srv_extgenes',
                                  input=output_from(srv_allgenes),
                                  filter=suffix('.bed.gz'),
                                  output='_ext1kb.bed.gz',
                                  output_dir=dir_allgenes,
                                  extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'srv_cgidist')
    srv_cgidist = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='srv_cgidist',
                                 input=output_from(srv_extgenes),
                                 filter=formatter(),
                                 output=os.path.join(dir_allgenes, 'hg19_CGI_gene_distal.bed'),
                                 extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'srv_cgiprox')
    srv_cgiprox = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='srv_cgiprox',
                                 input=output_from(srv_extgenes),
                                 filter=formatter(),
                                 output=os.path.join(dir_allgenes, 'hg19_CGI_gene_proximal.bed'),
                                 extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'srv_anncgi')
    srv_anncgi = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='srv_anncgi',
                                input=output_from(srv_cgidist, srv_cgiprox),
                                filter=suffix('.bed'),
                                output='_feat.bed',
                                output_dir=dir_allgenes,
                                extras=[cmd, jobcall])

    run_project_sarvesh = pipe.merge(task_func=touch_checkfile,
                                     name='project_sarvesh',
                                     input=output_from(proj_srv_init, srv_pcgenes, srv_2kbwin,
                                                       srv_cgiprom, srv_noncgi, srv_annreg,
                                                       srv_allgenes, srv_extgenes, srv_cgidist, srv_cgiprox,
                                                       srv_anncgi),
                                     output=os.path.join(dir_project_sarvesh, 'run_project_sarvesh.chk'))

    #
    # End of project Sarvesh
    # ======================

    # ======================
    # Pairwise alignment information
    # between human (hg19) and mouse (mm10)
    # including gonosomes for DEEP

    dir_project_deep = os.path.join(workdir, 'projects', 'deep')
    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    deep_chainfile = os.path.join(workdir, 'chainfiles', 'rawdata', 'hg19ToMm10.over.chain.gz')

    deep_aln_init = pipe.originate(task_func=lambda x: x,
                                   name='deep_aln_init',
                                   output=[deep_chainfile])

    deep_prefiltchain_dir = os.path.join(dir_project_deep, 'temp', 'prefilter')
    cmd = config.get('Pipeline', 'deep_chfprefilt').replace('\n', ' ')
    deep_chainfilter_params = build_chain_filter_commands([deep_chainfile], dir_out_chromaugo,
                                                          deep_prefiltchain_dir, cmd, jobcall)
    deep_prefilter_chain = pipe.files(sci_obj.get_jobf('in_out'),
                                      deep_chainfilter_params,
                                      name='deep_prefilter_chain')
    deep_prefilter_chain = deep_prefilter_chain.mkdir(deep_prefiltchain_dir)
    deep_prefilter_chain = deep_prefilter_chain.follows(deep_aln_init)

    chain_re = '(?P<TARGET>\w+)_to_(?P<QUERY>\w+)(?P<EXT>\.[\w\.]+)'

    cmd = config.get('Pipeline', 'chfswap')
    deep_swap_dir = os.path.join(dir_project_deep, 'temp', 'swap')
    deep_chf_swap = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                   name='deep_chf_swap',
                                   input=output_from(deep_prefilter_chain),
                                   filter=formatter(chain_re),
                                   output=os.path.join(deep_swap_dir, '{QUERY[0]}_to_{TARGET[0]}.tbest.chain.gz'),
                                   extras=[cmd, jobcall])
    deep_chf_swap = deep_chf_swap.mkdir(deep_swap_dir)

    cmd = config.get('Pipeline', 'deep_qrybnet')
    deep_rbestnet_dir = os.path.join(dir_project_deep, 'temp', 'rbest_net')
    deep_qrybnet = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                  name='deep_qrybnet',
                                  input=output_from(deep_chf_swap),
                                  filter=formatter(chain_re),
                                  output=os.path.join(deep_rbestnet_dir, '{TARGET[0]}_to_{QUERY[0]}.rbest.net.gz'),
                                  extras=[cmd, jobcall])
    deep_qrybnet = deep_qrybnet.mkdir(deep_rbestnet_dir)

    cmd = config.get('Pipeline', 'deep_qrybchain')
    deep_rbestchain_dir = os.path.join(dir_project_deep, 'temp', 'rbest_chain')
    deep_qrybchain = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                    name='deep_qrybchain',
                                    input=output_from(deep_qrybnet),
                                    filter=formatter(chain_re),
                                    output=os.path.join(deep_rbestchain_dir, '{TARGET[0]}_to_{QUERY[0]}.rbest.chain.gz'),
                                    extras=[cmd, jobcall])
    deep_qrybchain = deep_qrybchain.mkdir(deep_rbestchain_dir)

    cmd = config.get('Pipeline', 'trgbchain')
    deep_trgbchain = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                    name='deep_trgbchain',
                                    input=output_from(deep_qrybchain),
                                    filter=formatter(chain_re),
                                    output=os.path.join(deep_rbestchain_dir, '{QUERY[0]}_to_{TARGET[0]}.rbest.chain.gz'),
                                    extras=[cmd, jobcall])
    deep_trgbchain = deep_trgbchain.mkdir(deep_rbestchain_dir)

    cmd = config.get('Pipeline', 'deep_trgbfilt').replace('\n', ' ')
    deep_rbestfilter_dir = os.path.join(dir_project_deep, 'temp', 'rbest_chain_filt')
    deep_trgbfilt = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                   name='deep_trgbfilt',
                                   input=output_from(deep_trgbchain),
                                   filter=formatter(chain_re),
                                   output=os.path.join(deep_rbestfilter_dir, '{TARGET[0]}_to_{QUERY[0]}.rbest.chain.filt.gz'),
                                   extras=[cmd, jobcall])
    deep_trgbfilt = deep_trgbfilt.mkdir(deep_rbestfilter_dir)

    cmd = config.get('Pipeline', 'trgbbed').replace('\n', ' ')
    deep_rbestbed_dir = os.path.join(dir_project_deep, 'rbest_chain_bed')
    deep_trgbbed = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                  name='deep_trgbbed',
                                  input=output_from(deep_trgbfilt),
                                  filter=formatter(chain_re),
                                  output=os.path.join(deep_rbestbed_dir, '{TARGET[0]}_to_{QUERY[0]}.rbest.bed'),
                                  extras=[cmd, jobcall])
    deep_trgbbed = deep_trgbbed.mkdir(deep_rbestbed_dir)

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('CrossmapEnv')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'deep_crossmap').replace('\n', ' ')
    deep_crossmap = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                   name='deep_crossmap',
                                   input=output_from(deep_trgbbed),
                                   filter=formatter(chain_re),
                                   output=os.path.join(deep_rbestbed_dir, '{TARGET[0]}_to_{QUERY[0]}.mapped.bed'),
                                   extras=[cmd, jobcall])

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'deep_mrgblocks').replace('\n', ' ')
    deep_pwaln_dir = os.path.join(dir_project_deep, 'aln')
    deep_mrgblocks = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                    name='deep_mrgblocks',
                                    input=output_from(deep_crossmap),
                                    filter=formatter(chain_re),
                                    output=os.path.join(deep_pwaln_dir, '{TARGET[0]}_to_{QUERY[0]}.aln.tsv.gz'),
                                    extras=[cmd, jobcall])
    deep_mrgblocks = deep_mrgblocks.mkdir(deep_pwaln_dir)

    cmd = config.get('Pipeline', 'deep_trgbed').replace('\n', ' ')
    deep_trgbed_dir = os.path.join(dir_project_deep, 'temp', 'alnbed')
    deep_trgbed = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='deep_trgbed',
                                 input=output_from(deep_mrgblocks),
                                 filter=suffix('aln.tsv.gz'),
                                 output='hg19.bed.gz',
                                 output_dir=deep_trgbed_dir,
                                 extras=[cmd, jobcall])
    deep_trgbed = deep_trgbed.mkdir(deep_trgbed_dir)

    cmd = config.get('Pipeline', 'deep_mapaln').replace('\n', ' ')
    deep_mapaln = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='deep_mapaln',
                                 input=output_from(deep_trgbed),
                                 filter=formatter(chain_re),
                                 output=os.path.join(deep_trgbed_dir, '{TARGET[0]}_to_{QUERY[0]}.map.bed.gz'),
                                 extras=[cmd, jobcall])

    #
    # End of: alignment hg19 - mm10
    # ======================


    dir_projects = os.path.join(workdir, 'projects')
    run_projects = pipe.merge(task_func=touch_checkfile,
                              name='run_projects',
                              input=output_from(run_project_sarvesh),
                              output=os.path.join(dir_projects, 'run_all_projects.chk'))

    return pipe
