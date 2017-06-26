#!/usr/bin/env python
# coding=utf-8

import os as os
import sys as sys
import pickle as pck
import traceback as trb
import argparse as argp
import itertools as itt
import gzip as gz
import collections as col
import re as re
import xmlrpc.client as cli

__CACHE_DIR__ = '/TL/deep/fhgfs/projects/pebert/thesis/refdata/lola/rawdata'
__DOWNLAOD_BASE__ = '/TL/deep/fhgfs/projects/pebert/thesis/refdata/lola/rawdata/download'
__DEEPBLUE_USERKEY__ = '/home/pebert/.ssh/id_deepblue.private'


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--user-key', '-uk', type=str, dest='userkey', default=__DEEPBLUE_USERKEY__)
    parser.add_argument('--cache-dir', '-cd', type=str, dest='cachedir', default=__CACHE_DIR__)
    parser.add_argument('--dl-base', '-dl', type=str, dest='dlbase', default=__DOWNLAOD_BASE__)
    parser.add_argument('--genomes', '-g', type=str, nargs='+', dest='genomes',
                        default=['hg19', 'mm9'])
    parser.add_argument('--projects', '-p', type=str, nargs='+', dest='projects',
                        default=['ENCODE', 'Roadmap Epigenomics'])
    parser.add_argument('--techniques', '-t', type=str, nargs='+', dest='techniques',
                        default=['ChIP-seq', 'DNase-seq'])
    parser.add_argument('--force', '-f', action='store_true', default=False, dest='force')
    args = parser.parse_args()
    return args


def init_deepblue_conn(userkey):
    """
    :param userkey:
    :return:
    """
    srv = cli.ServerProxy('http://deepblue.mpi-inf.mpg.de/xmlrpc',
                          encoding='utf-8', allow_none=True)
    s, r = srv.echo(userkey)
    assert s == 'okay', 'Server conn init failed: {}'.format(r)
    return srv


def query_experiments(dbsrv, userkey, cache_file, args):
    """
    :param dbsrv:
    :param userkey:
    :param cache_file:
    :param args:
    :return:
    """
    if os.path.isfile(cache_file):
        with open(cache_file, 'rb') as dump:
            cache = pck.load(dump)
        return cache
    cache = dict()
    for g, t, p in itt.product(args.genomes, args.techniques, args.projects):
        s, r = dbsrv.list_experiments(g, 'peaks', '', '', '', t, p, userkey)
        if s != 'okay':
            sys.stderr.write('\nRequest {} - {} - {} failed: {}\n'.format(g, t, p, r))
            continue
        cache[(g, t, p)] = r
    with open(cache_file, 'wb') as out:
        pck.dump(cache, out)
    return cache


def query_exp_info(dbsrv, userkey, cache_file, experiments):
    """
    :param dbsrv:
    :param userkey:
    :param cache_file:
    :param experiments:
    :return:
    """
    if os.path.isfile(cache_file):
        with open(cache_file, 'rb') as dump:
            cache = pck.load(dump)
    else:
        cache = dict()
    for k, v in experiments.items():
        for eid, fn in v:
            try:
                _ = cache[eid]
            except KeyError:
                s, r = dbsrv.info(eid, userkey)
                if s != 'okay':
                    sys.stderr.write('\nInfo request {} failed: {}\n'.format(eid, r))
                cache[eid] = r[0]
    with open(cache_file, 'wb') as out:
        pck.dump(cache, out)
    return cache


def check_sample_filters(sample, project):
    """
    :param sample:
    :param project:
    :return:
    """
    if project == 'ENCODE':
        bt = sample['biosample_type']
        if bt.startswith('immortal') or bt.startswith('in vitro'):
            return False
        btn = (sample['biosample_term_name']).lower()
        if btn == 'naive b cell':
            return False
        if any([x in btn for x in ['right', 'left', 'female', 'fat',
                                   'induced', 'inferior', 'superior',
                                   'uterus', 'putamen', 'l1-s8', 'gland',
                                   'testis', 'ovary', 'prostate', 'psoas',
                                   'horn', 'pons', 'gyrus', 'nucleus', 'pedal',
                                   'abdomen', 'amnion', 'chorion', 'amniotic',
                                   'mononuclear', 'langerhans', 'cortex', 'cortical',
                                   'occipital', 'pallidus', 'olfactory', 'hematopoietic',
                                   'telencephalon']]):
            return False
        if 'fibroblast' in btn and 'product_id' in sample:
            # numerous fibroblast cell lines
            return False
        if 'life_stage' in sample:
            ls = (sample['life_stage']).lower()
            if any([x in ls for x in ['adult', 'unknown']]):
                return True
        return False
    elif project == 'Roadmap Epigenomics':
        bt = sample['TYPE']
        if bt in ['CellLine', 'ESCDerived']:
            return False
        if 'AGE' in sample:
            age = (sample['AGE']).lower()
            if any([x in age for x in ['fetus', 'child', 'postnatal',
                                       '0y', 'newborn', 'fetal']]):
                return False
        edn = (sample['EDACC_NAME']).lower()
        if any([x in edn for x in ['fetal', 'stimulated', 'gm12878',
                                   'k562', 'induced']]):
            return False
        std = (sample['STD_NAME']).lower()
        if any([x in std for x in ['ips', 'neurospheres', 'derived',
                                   'variant', 'left', 'right', 'nuclei',
                                   'psoas', 'gyrus', 'pedal', 'imr90',
                                   'abdomen', 'amnion', 'chorion']]):
            return False
        return True
    else:
        raise ValueError('Unknown project: {}'.format(project))


def determine_assay_group(tech, mark):
    """
    :param tech:
    :param mark:
    :return:
    """
    hist_match = '^(H3K|H2A|H2B|H4K|H3ac)'
    if tech == 'DNase-seq':
        return 'dnase'
    elif mark.startswith('eGFP') or mark == 'H3T11ph':
        # the histone mark exists only once for H9, ignore
        return 'skip'
    elif re.match(hist_match, mark) is not None:
        return 'histone'
    else:
        return 'tf'


def check_filetype_filters(extras, project):
    """
    :param extras:
    :param project:
    :return:
    """
    if project == 'ENCODE':
        ft = extras['file_type']
    elif project == 'Roadmap Epigenomics':
        ft = extras['type']
    else:
        raise ValueError('Unknown project: {}'.format(project))
    if any([x in ft for x in ['bigBed', 'bed9', 'gappedPeak']]):
        return True
    if ft == 'macs2':
        if 'narrowPeak' not in extras['url']:
            return True
    if ft == 'hotspot':
        if 'all.peaks' in extras['url']:
            return True
    return False


def filter_samples(expinfo):
    """
    :param expinfo:
    :return:
    """
    keepers = dict()
    group_counts = col.Counter()
    for eid, infos in expinfo.items():
        project = infos['project']
        skip = check_filetype_filters(infos['extra_metadata'], project)
        if skip:
            continue
        if 'treated' in infos['description'].lower():
            continue
        t, m = infos['technique'], infos['epigenetic_mark']
        if m == 'H2A.Z':
            m = 'H2AZ'
        grp = determine_assay_group(t, m)
        if grp == 'skip':
            continue
        sid, smp_infos = infos['sample_id'], infos['sample_info']
        try:
            keep = check_sample_filters(smp_infos, project)
            if keep:
                group_key = infos['genome'], project, grp
                infos['assay_group'] = grp
                infos['local_group'] = group_key
                keepers[eid] = infos
                group_counts[group_key] += 1
        except KeyError as ke:
            print(project)
            print('=====')
            for k, v in sorted(smp_infos.items()):
                print(k, ' - ', v)
            raise ke
    return keepers


def make_filenames(expinfo):
    """
    :param expinfo:
    :return:
    """
    dup_counter = col.Counter()
    for eid, infos in expinfo.items():
        project, assm, mark = infos['project'], infos['genome'], infos['epigenetic_mark']
        if mark == 'H2A.Z':
            mark = 'H2AZ'
        if project == 'Roadmap Epigenomics':
            project = 'REMC'
        project = project
        sid, smp = infos['sample_id'], infos['sample_info']
        infos['local_path'] = os.path.join(assm, (project + '_' + infos['assay_group']).lower())
        if project == 'REMC':
            bioname = smp['MNEMONIC']
            bioname = bioname.replace('.', '-')
            remc_id = smp['EID']
            infos['cell_type'] = bioname.split('-', 1)[-1].lower()
            infos['tissue'] = smp['ANATOMY'].lower().replace('_', '-')
            ct, ti = infos['cell_type'], infos['tissue']
            new_name = '_'.join([project, remc_id, assm, ti, ct, mark])
            dup_counter[new_name] += 1
            rep = dup_counter[new_name]
            infos['new_filename'] = new_name + '.r{}.bed.gz'.format(rep)
            infos['checkfile'] = new_name + '.r{}.dl.chk'.format(rep)
        else:
            enc_tissues = ['small intestine', 'large intestine', 'sigmoid colon', 'germinal center',
                           'bone marrow', 'transverse colon', 'esophagus mucosa']
            enc_lut = {'splenic B cell': ('bcell', 'spleen'), 'B cell': ('bcell', 'cord-blood'),
                       'cardiac muscle cell': ('muscle', 'heart'),
                       'bone marrow macrophage': ('macrophage', 'bone-marrow'),
                       'embryonic fibroblast': ('fibroblast', 'embryonic'),
                       'stromal cell bone marrow': ('stromal', 'bone-marrow'),
                       'epidermal melanocyte': ('melanocyte', 'epidermis'),
                       'hepatic stellate cell': ('stellate', 'liver'), 'brain pericyte': ('pericyte', 'brain'),
                       'skeletal muscle myoblast': ('myoblast', 'skeletal-muscle'),
                       'CD1cp myeloid dendritic cell': ('cd1cp-dendritic', 'blood'),
                       'smooth muscle cell brain vasculature': ('smooth-muscle', 'brain'),
                       'regulatory T-lymphocyte': ('Treg', 'blood'), 'CD4p helper T cell': ('cd4p-Th', 'blood'),
                       'CD14p monocyte': ('cd14p-monocyte', 'blood'),
                       'CD4p CD25p alpha-beta regulatory T cell': ('cd4p-cd25p-ab-Treg', 'blood'),
                       'naive thymus-derived CD4p alpha-beta T cell': ('cd4p-ab-Tnaive', 'thymus'),
                       'erythroblast': ('erythroblast', 'blood')}
            btn = smp['biosample_term_name']
            btn = btn.replace(' of the ', ' ')
            btn = btn.replace(' of ', ' ')
            btn = btn.replace('-positive', 'p')
            btn = btn.replace('-negative', 'n')
            btn = btn.replace(',', ' ')
            btn = btn.replace('  ', ' ')

            if btn in enc_lut:
                ct, ti = enc_lut[btn]
                infos['cell_type'] = ct
                infos['tissue'] = ti
            elif len(btn.split()) == 1:
                infos['tissue'] = btn
                infos['cell_type'] = 'whole-sample'
            else:
                if btn.endswith('epithelial cell'):
                    cell_type = 'epithelial'
                    tissue = btn.rsplit(' ', 2)[0].replace(' ', '-')
                    infos['tissue'] = tissue
                    infos['cell_type'] = cell_type
                elif btn.startswith('epithelial cell'):
                    cell_type = 'epithelial'
                    tissue = btn.split(' ', 2)[-1].replace(' ', '-')
                    infos['tissue'] = tissue
                    infos['cell_type'] = cell_type
                elif btn.endswith('endothelial cell'):
                    cell_type = 'endothelial'
                    tissue = btn.rsplit(' ', 2)[0].replace(' ', '-')
                    infos['tissue'] = tissue
                    infos['cell_type'] = cell_type
                elif 'tissue' in btn:
                    infos['cell_type'] = 'whole-sample'
                    infos['tissue'] = btn.replace(' ', '-')
                elif btn in enc_tissues:
                    cell_type = 'whole-sample'
                    tissue = btn.replace(' ', '-')
                    infos['cell_type'] = cell_type
                    infos['tissue'] = tissue
                elif 'astrocyte' in btn:
                    infos['cell_type'] = 'astrocyte'
                    tissue = btn.replace('astrocyte ', '').replace(' ', '-')
                    infos['tissue'] = tissue
                elif btn.startswith('T-helper'):
                    infos['tissue'] = 'blood'
                    cell_type = 'Th{}-Tcell'.format(btn.split()[1])
                    infos['cell_type'] = cell_type
                else:
                    raise ValueError('Unknown biotype: {}'.format(btn))
            ct, ti = infos['cell_type'], infos['tissue']
            new_name = '_'.join([project, sid, assm, ti, ct, mark])
            dup_counter[new_name] += 1
            rep = dup_counter[new_name]
            infos['new_filename'] = new_name + '.r{}.bed.gz'.format(rep)
            infos['checkfile'] = new_name + '.r{}.dl.chk'.format(rep)
    for eid, infos in expinfo.items():
        assert 'cell_type' in infos, 'Cell type missing: {}'.format(eid)
        assert 'tissue' in infos, 'Tissue missing: {}'.format(eid)
    return expinfo


def download_regions(srv, userkey, expinfo, dlbase):
    """
    :param srv:
    :param userkey:
    :param expinfo:
    :return:
    """
    for eid, infos in expinfo.items():
        sub_path = os.path.join(dlbase, infos['local_path'])
        if not os.path.isdir(sub_path):
            os.makedirs(sub_path, exist_ok=True)
        checkfile = os.path.join(sub_path, infos['checkfile'])
        if os.path.isfile(checkfile):
            continue
        datafile = os.path.join(sub_path, infos['new_filename'])
        stat, res = srv.select_experiments(eid, '', 0, 300000000, userkey)
        if stat == 'okay':
            qid = res
            infos['query_id'] = qid
            stat, res = srv.get_regions(qid, 'CHROMOSOME,START,END,@NAME,@SAMPLE_ID,@EPIGENETIC_MARK', userkey)
            if 'error' in stat:
                infos['request_id'] = 'error - {}'.format(res)
                continue
            else:
                rid = res
                proc_stat = 'running'
                while proc_stat != 'done':
                    stat, res = srv.info(rid, userkey)
                    proc_stat = res[0]['state']
                infos['request_id'] = rid
                stat, res = srv.get_request_data(rid, userkey)
                if stat != 'okay':
                    infos['request_dl'] = 'error - {}'.format(res)
                    continue
                else:
                    infos['request_dl'] = 'okay'
                    with gz.open(datafile, 'wt') as dump:
                        _ = dump.write(res + '\n')
                    print('Dumped: ', datafile)
                    with open(checkfile, 'w') as check:
                        _ = check.write('')
        else:
            infos['query_id'] = 'error - {}'.format(res)
    return expinfo


def main(args, db, uk):
    """
    :param args:
    :return:
    """
    cd = args.cachedir
    dl = args.dlbase

    cache_explist = os.path.join(cd, 'db_cache_explist.pck')
    if args.force and os.path.isfile(cache_explist):
        os.unlink(cache_explist)
    exp_collection = query_experiments(db, uk, cache_explist, args)
    cache_expinfo = os.path.join(cd, 'db_cache_expinfo.pck')
    if args.force and os.path.isfile(cache_expinfo):
        os.unlink(cache_expinfo)
    exp_info = query_exp_info(db, uk, cache_expinfo, exp_collection)
    exp_info = filter_samples(exp_info)
    exp_info = make_filenames(exp_info)
    cache_prep_info = os.path.join(cd, 'db_cache_final.pck')
    with open(cache_prep_info, 'wb') as dump:
        pck.dump(exp_info, dump)
    exp_info = download_regions(db, uk, exp_info, dl)
    cache_done = os.path.join(cd, 'db_cache_complete.pck')
    with open(cache_done, 'wb') as dump:
        pck.dump(exp_info, dump)
    return


if __name__ == '__main__':
    try:
        args = parse_command_line()
        assert os.path.isdir(args.cachedir), 'Invalid path to caching folder: {}'.format(args.cachedir)
        if os.path.isfile(args.userkey):
            with open(args.userkey, 'r') as keyfile:
                uk = keyfile.read().strip()
        else:
            uk = args.userkey
        srv = init_deepblue_conn(uk)
        main(args, srv, uk)
    except Exception as err:
        trb.print_exc()
        sys.stderr.write('\nError: {}\n'.format(err))
        sys.exit(1)
    else:
        sys.exit(0)
