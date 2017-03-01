# coding=utf-8

import os as os
import io as io
import datetime as dt

from ruffus import *

from pipelines.auxmod.auxiliary import read_chromsizes, open_comp, collect_full_paths
from pipelines.auxmod.enhancer import process_merged_encode_enhancer, process_vista_enhancer
from pipelines.auxmod.cpgislands import process_ucsc_cgi
from pipelines.auxmod.chainfiles import build_chain_filter_commands, build_symm_filter_commands
from pipelines.auxmod.bedroi import make_bed_roi, normalize_join
from pipelines.auxmod.orthologs import match_ortholog_files
from pipelines.auxmod.promoter import split_promoter_files, score_promoter_regions

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
    gm_pc_win5p = pipe.transform(task_func=make_bed_roi,
                                 name='gm_pc_win5p',
                                 input=output_from(gm_pc_hg19, gm_pc_mm9, gm_pc_bta7, gm_pc_ens),
                                 filter=formatter('(?P<ANNOTID>[a-z]{3}_(?P<ASSM>[a-zA-Z0-9]+)_\w+)\.genes\.bed\.gz'),
                                 output=os.path.join(pc_roi_dump, '{ANNOTID[0]}.win5p.bed.gz'),
                                 extras=[os.path.join(dir_out_chromaugo, '{ASSM[0]}_chrom_augo.tsv'),
                                         'win5p']).mkdir(pc_roi_dump).jobs_limit(4)

    gm_pc_reg5p = pipe.transform(task_func=make_bed_roi,
                                 name='gm_pc_reg5p',
                                 input=output_from(gm_pc_hg19, gm_pc_mm9, gm_pc_bta7, gm_pc_ens),
                                 filter=formatter('(?P<ANNOTID>[a-z]{3}_(?P<ASSM>[a-zA-Z0-9]+)_\w+)\.genes\.bed\.gz'),
                                 output=os.path.join(pc_roi_dump, '{ANNOTID[0]}.reg5p.bed.gz'),
                                 extras=[os.path.join(dir_out_chromaugo, '{ASSM[0]}_chrom_augo.tsv'),
                                         'reg5p']).mkdir(pc_roi_dump).jobs_limit(4)

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
                               input=output_from(gm_pc_win5p, gm_pc_reg5p, gm_pc_body),
                               filter=formatter('(?P<SPEC>\w+)_(?P<ASSM>\w+)_(?P<ANNOT>\w+)_'
                                                '(?P<VER>\w+)\.(?P<REGTYPE>\w+)\.bed\.gz'),
                               output=os.path.join(pc_roi_conv_hdf, '{SPEC[0]}_{ASSM[0]}_{ANNOT[0]}_'
                                                                    '{VER[0]}.{REGTYPE[0]}.h5'),
                               extras=[cmd, jobcall]).mkdir(pc_roi_conv_hdf)

    run_task_genemodel = pipe.merge(task_func=touch_checkfile,
                                    name='task_genemodel',
                                    input=output_from(gm_gtftobed, gm_gfftobed, gm_enstobed,
                                                      gm_pc_hg19, gm_pc_mm9, gm_pc_bta7, gm_pc_ens,
                                                      gm_pc_win5p, gm_pc_reg5p, gm_pc_body,
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

    dir_tm_mapindex = os.path.join(dir_task_transmodel, 'mapindex')

    # k=31 - read length ~75 and above; Salmon does not support higher values for k
    # k=31 is the developer's recommendation for read length of 75
    cmd = config.get('Pipeline', 'genidx31').replace('\n', ' ')
    tm_idxgen31 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='tm_idxgen31',
                                 input=output_from(tm_fasta),
                                 filter=formatter('(?P<TRANSCRIPTOME>(hsa|mmu)_\w+)\.transcripts\.fa\.gz'),
                                 output='{subpath[0][1]}/qindex/{TRANSCRIPTOME[0]}.k31.idx/rsd.bin',
                                 extras=[cmd, jobcall]).mkdir(dir_tm_mapindex)

    # k=19 - read length ~50
    cmd = config.get('Pipeline', 'genidx19').replace('\n', ' ')
    tm_idxgen19 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='tm_idxgen19',
                                 input=output_from(tm_fasta),
                                 filter=formatter('(?P<TRANSCRIPTOME>(cfa|bta|mmu|ssc)_\w+)\.transcripts\.fa\.gz'),
                                 output='{subpath[0][1]}/qindex/{TRANSCRIPTOME[0]}.k19.idx/rsd.bin',
                                 extras=[cmd, jobcall]).mkdir(dir_tm_mapindex)

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
    ortho_files = collect_full_paths(ortho_rawdata, '*.txt.gz')
    ortho_init = pipe.originate(task_func=lambda x: x,
                                name='ortho_init',
                                output=ortho_files)

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    subset_files = collect_full_paths(pc_subset, '*.genes.*')
    genemodel_dir = bedout  # NB: if ever: should be renamed to gm_bedout for clarity
    cmd = config.get('Pipeline', 'orthconv').replace('\n', ' ')
    outdir = os.path.join(dir_task_orthologs, 'hdf')
    ortho_conv = pipe.files(sci_obj.get_jobf('in_out'),
                            match_ortholog_files(ortho_files, genemodel_dir, subset_files, outdir, cmd, jobcall),
                            name='ortho_conv').mkdir(outdir)

    orthomam_init = pipe.originate(task_func=lambda x: x,
                                   name='orthomam_init',
                                   output=os.path.join(ortho_rawdata, 'orthomam_v9_5vert_orthologs.zip'))

    cmd = config.get('Pipeline', 'orthomamraw').replace('\n', ' ')
    orthomam_raw = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                  name='orthomam_raw',
                                  input=output_from(orthomam_init),
                                  filter=formatter(),
                                  output=os.path.join(ortho_rawdata, 'orthomam_v9_5vert_orthologs.tsv'),
                                  extras=[cmd, jobcall])

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

    run_task_enh = pipe.merge(task_func=touch_checkfile,
                              name='run_task_enh',
                              input=output_from(enh_super,
                                                enh_hsa_fantom, enh_mmu_fantom,
                                                enh_cat_encode, enh_mrg_encode, enh_hsa_encpp, enh_hsa_encdp,
                                                enh_hsa_vista, enh_mmu_vista),
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

    chf_rawdata = os.path.join(dir_task_chainfiles, 'rawdata')
    raw_chainfiles = collect_full_paths(chf_rawdata, '*.chain.gz')
    chf_init = pipe.originate(task_func=lambda x: x,
                              name='chf_init',
                              output=raw_chainfiles)

    filtered_chains = os.path.join(dir_task_chainfiles, 'filtered')
    cmd = config.get('Pipeline', 'chfilt').replace('\n', ' ')
    chf_filter = pipe.files(sci_obj.get_jobf('in_out'),
                            build_chain_filter_commands(raw_chainfiles, dir_out_chromauto,
                                                        filtered_chains, cmd, jobcall),
                            name='chf_filter').mkdir(filtered_chains)

    chain_re = '(?P<TARGET>\w+)_to_(?P<QUERY>\w+)(?P<EXT>\.[\w\.]+)'

    cmd = config.get('Pipeline', 'chfswap')
    chf_swap_dir = os.path.join(dir_task_chainfiles, 'tmp_swap')
    chf_swap = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='chf_swap',
                              input=output_from(chf_filter),
                              filter=formatter(chain_re),
                              output=os.path.join(chf_swap_dir, '{QUERY[0]}_to_{TARGET[0]}.tbest.chain.gz'),
                              extras=[cmd, jobcall]).mkdir(chf_swap_dir)

    cmd = config.get('Pipeline', 'qrybnet')
    chf_rbest_net = os.path.join(dir_task_chainfiles, 'rbest_net')
    qrybnet = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='qrybnet',
                             input=output_from(chf_swap),
                             filter=formatter(chain_re),
                             output=os.path.join(chf_rbest_net, '{TARGET[0]}_to_{QUERY[0]}.rbest.net.gz'),
                             extras=[cmd, jobcall]).mkdir(chf_rbest_net)

    cmd = config.get('Pipeline', 'qrybchain')
    chf_rbest_chain = os.path.join(dir_task_chainfiles, 'rbest_chain')
    qrybchain = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                               name='qrybchain',
                               input=output_from(qrybnet),
                               filter=formatter(chain_re),
                               output=os.path.join(chf_rbest_chain, '{TARGET[0]}_to_{QUERY[0]}.rbest.chain.gz'),
                               extras=[cmd, jobcall]).mkdir(chf_rbest_chain)

    cmd = config.get('Pipeline', 'trgbchain')
    trgbchain = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                               name='trgbchain',
                               input=output_from(qrybchain),
                               filter=formatter(chain_re),
                               output=os.path.join(chf_rbest_chain, '{QUERY[0]}_to_{TARGET[0]}.rbest.chain.gz'),
                               extras=[cmd, jobcall]).mkdir(chf_rbest_chain)

    cmd = config.get('Pipeline', 'trgbnet')
    trgbnet = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='trgbnet',
                             input=output_from(trgbchain),
                             filter=formatter(chain_re),
                             output=os.path.join(chf_rbest_net, '{TARGET[0]}_to_{QUERY[0]}.rbest.net.gz'),
                             extras=[cmd, jobcall]).mkdir(chf_rbest_net)

    cmd = config.get('Pipeline', 'bednet')
    chf_rbest_bed = os.path.join(dir_task_chainfiles, 'rbest_bed')
    chf_bednet = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='chf_bednet',
                                input=output_from(trgbnet, qrybnet),
                                filter=suffix('.rbest.net.gz'),
                                output='.rbest.bed.gz',
                                output_dir=chf_rbest_bed,
                                extras=[cmd, jobcall]).mkdir(chf_rbest_bed)

    cmd = config.get('Pipeline', 'qrysymm')
    symmfilt_chains = os.path.join(dir_task_chainfiles, 'rbest_symm')
    qrysymm = pipe.files(sci_obj.get_jobf('in_out'),
                         build_symm_filter_commands(collect_full_paths(chf_rbest_chain, '*.chain.gz'),
                                                    dir_out_chromauto, symmfilt_chains, cmd, jobcall),
                         name='qrysymm').mkdir(symmfilt_chains)

    cmd = config.get('Pipeline', 'mrgblocks')
    mrgblocks = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                             name='mrgblocks',
                             input=output_from(qrysymm),
                             filter=formatter('(?P<MAPPING>\w+)\.(?P<CHROM>\w+)\.symmap\.tsv\.gz'),
                             output=os.path.join(symmfilt_chains, '{MAPPING[0]}.wg.symmap.tsv.gz'),
                             extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'normblocks')
    normblocks_dir = os.path.join(dir_task_chainfiles, 'rbest_norm')
    normblocks = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='normblocks',
                                input=output_from(mrgblocks),
                                filter=suffix('wg.symmap.tsv.gz'),
                                output='norm.symmmap.tsv.gz',
                                output_dir=normblocks_dir,
                                extras=[cmd, jobcall]).mkdir(normblocks_dir)

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'mapidx').replace('\n', ' ')
    hdfmap_dir = os.path.join(dir_task_chainfiles, 'hdf_map')
    mapidx = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='mapidx',
                            input=output_from(normblocks),
                            filter=formatter(chain_re),
                            output=os.path.join(hdfmap_dir, '{TARGET[0]}_to_{QUERY[0]}.map.h5'),
                            extras=[cmd, jobcall]).mkdir(hdfmap_dir)

    run_task_chains = pipe.merge(task_func=touch_checkfile,
                                 name='task_chains',
                                 input=output_from(chf_init, chf_filter, chf_swap,
                                                   qrybnet, qrybchain, trgbchain, trgbnet,
                                                   chf_bednet, qrysymm, mrgblocks, normblocks,
                                                   mapidx),
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
    raw_consfiles = collect_full_paths(cons_rawdata, '*Elem*.gz')
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

    run_task_conservation = pipe.merge(task_func=touch_checkfile,
                                       name='task_cons',
                                       input=output_from(cons_init, cons_elemtobed),
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


    run_all = pipe.merge(task_func=touch_checkfile,
                         name='run_all',
                         input=output_from(run_task_csz, run_task_enh, run_task_cgi,
                                           run_task_genemodel, run_task_transmodel,
                                           run_task_chains),
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

    dir_projects = os.path.join(workdir, 'projects')
    run_projects = pipe.merge(task_func=touch_checkfile,
                              name='run_projects',
                              input=output_from(run_project_sarvesh),
                              output=os.path.join(dir_projects, 'run_all_projects.chk'))

    return pipe