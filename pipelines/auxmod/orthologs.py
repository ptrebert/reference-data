# coding=utf-8

import os as os
import fnmatch as fnm


def match_ortholog_files(orthos, genemodels, subsets, outdir, cmd, jobcall):
    """
    :param orthos:
    :param genemodels:
    :param subsets:
    :return:
    """
    spec_assm = {'human': 'hg19', 'mouse': 'mm9', 'pig': 'susScr2',
                 'cow': 'bosTau7', 'dog': 'canFam3'}
    spec_model = {'human': 'gencode.v19.annotation.bed.gz', 'mouse': 'gencode.vM1.annotation.bed.gz',
                  'pig': 'ssc_susScr2_ensembl_v63.bed.gz', 'dog': 'cfa_canFam3_ensembl_v81.bed.gz',
                  'cow': 'Ensembl75_liftOver_Btau_4.6.1_genes.bed.gz'}
    args = []
    allow_skip = False
    for orth in orthos:
        fp, fn = os.path.split(orth)
        s1, s2 = fn.split('_')[:2]
        if s1 in spec_assm and s2 in spec_assm:
            gm1 = os.path.join(genemodels, spec_model[s1])
            gm2 = os.path.join(genemodels, spec_model[s2])
            assm1 = spec_assm[s1]
            sub1 = fnm.filter(subsets, '*_{}_*genes.bed.gz'.format(assm1))
            if len(sub1) < 1:
                allow_skip = True
                continue
            assert len(sub1) == 1, 'Ambig. subset file for {}: {}'.format(s1, sub1)
            sub1 = sub1[0]
            assm2 = spec_assm[s2]
            sub2 = fnm.filter(subsets, '*_{}_*genes.bed.gz'.format(assm2))
            if len(sub2) < 1:
                allow_skip = True
                continue
            assert len(sub2) == 1, 'Ambig. subset file for {}: {}'.format(s2, sub2)
            sub2 = sub2[0]
            outfile = '_'.join([assm1, assm2, 'hcop', 'orthologs']) + '.h5'
            tmp = cmd.format(**{'genemodel1': gm1, 'genemodel2': gm2,
                                'subset1': sub1, 'subset2': sub2})
            outpath = os.path.join(outdir, outfile)
            args.append([orth, outpath, tmp, jobcall])
    if orthos and not allow_skip:
        assert args, 'No run arguments created for ortholog matching'
    return args
