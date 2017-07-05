#!/usr/bin/env python
# coding=utf-8

import os as os
import sys as sys
import argparse as argp
import pickle as pck
import traceback as trb
import multiprocessing as mp
import collections as col

import pandas as pd
import numpy as np
import sklearn as skl
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--index', '-i', type=str, dest='index')
    parser.add_argument('--data', '-d', type=str, dest='data')
    mgrp = parser.add_mutually_exclusive_group()
    mgrp.add_argument('--tune', action='store_true', dest='tune')
    mgrp.add_argument('--apply', action='store_true', dest='apply')
    mgrp.add_argument('--extract', action='store_true', dest='extract')
    parser.add_argument('--models', '-m', type=str, dest='models',
                        default='/TL/deep/fhgfs/projects/pebert/thesis/refdata/enhancer/temp')
    parser.add_argument('--workers', '-w', type=int, default=15, dest='workers')
    parser.add_argument('--output', '-o', type=str, dest='output')
    args = parser.parse_args()
    return args


def run_model_tuning(features, target, workers, model_type):
    """
    :param features:
    :param target:
    :return:
    """
    if model_type == 'rfcls':
        param_grid = {'bootstrap': [True], 'oob_score': [True],
                      'n_estimators': [750, 1250, 2000, 3000],
                      'min_samples_split': [2, 4, 16, 32]}
        model = RandomForestClassifier()
        score_method = 'f1_macro'
        label_one = target.sum()
        wt_one = 0.6 / label_one
        wt_zero = 0.4 / (target.size - label_one)
        smp_wt = np.array([wt_one if l == 1 else wt_zero for l in target], dtype=np.float32)
        fit_params = {'sample_weight': smp_wt}
        tuning = GridSearchCV(model, param_grid, scoring=score_method, pre_dispatch=workers,
                              cv=10, n_jobs=workers, refit=True, verbose=0, fit_params=fit_params)
    else:
        param_grid = {'bootstrap': [True], 'oob_score': [True],
                      'n_estimators': [500, 1000, 1500, 2000],
                      'min_samples_split': [2, 8, 32, 64]}
        model = RandomForestRegressor()
        score_method = 'r2'
        tuning = GridSearchCV(model, param_grid, scoring=score_method, pre_dispatch=workers,
                              cv=10, n_jobs=workers, refit=True, verbose=0)
    tuning = tuning.fit(features, target)
    return tuning


def load_datasets(fpath, indices):
    """
    :param fpath:
    :param indices:
    :return:
    """
    val_idx = indices['validation']
    training = []
    testing = []
    with pd.HDFStore(fpath, 'r') as hdf:
        for k in hdf.keys():
            train_idx = indices[k.strip('/')]
            test_idx = val_idx[k.strip('/')]
            data = hdf[k]
            train_sub = data.loc[train_idx[0]:train_idx[1], :]
            training.append(train_sub)
            test_sub = data.loc[test_idx[0]:, :]
            testing.append(test_sub)
    training = pd.concat(training, axis=0, ignore_index=True)
    testing = pd.concat(testing, axis=0, ignore_index=True)
    return training, testing


def train_models(args):
    """
    :return:
    """
    model_id = int(os.path.basename(args.index).split('.')[0].split('_')[-1])
    with open(args.index, 'rb') as cvindex:
        indices = pck.load(cvindex)
    training, validation = load_datasets(args.data, indices)
    sub_output = args.output.replace('.pck', '.h5').replace('rfcls', 'mlsplit')
    with pd.HDFStore(sub_output, 'w', complib='blosc', complevel=9) as hdf:
        hdf.put('training', training, format='table')
        hdf.put('validation', validation, format='table')
    feat_columns = sorted([c for c in training.columns if c.startswith('ft')])

    for mtype in ['rfcls', 'rfreg']:
        if mtype == 'rfcls':
            train_labels = np.array(training['output'] > 0, dtype=np.int8)
            valid_labels = np.array(validation['output'] > 0, dtype=np.int8)
        else:
            # this is to avoid tuning the regression models again,
            # their performance is insufficient, no point here...
            continue
            train_labels = np.log1p(training['output'])
            valid_labels = np.log1p(validation['output'])
        cv_result = run_model_tuning(training[feat_columns], train_labels, args.workers, mtype)
        best_model = cv_result.best_estimator_
        oob_score = float(cv_result.best_score_)
        label_one = valid_labels.sum()
        wt_one = 0.6 / label_one
        wt_zero = 0.4 / (valid_labels.size - label_one)
        smp_wt = np.array([wt_one if l == 1 else wt_zero for l in valid_labels], dtype=np.float32)
        valid_score = float(best_model.score(validation[feat_columns], valid_labels, sample_weight=smp_wt))
        model_output = args.output.replace('rfcls', mtype)
        with open(model_output, 'wb') as dump:
            pck.dump({'cv': cv_result, 'oob_score': oob_score, 'model': mtype,
                      'valid_score': valid_score, 'model_id': model_id}, dump)
    return


def load_models(folder):
    """
    :param folder:
    :return:
    """
    all_files = os.listdir(folder)
    model_files = [os.path.join(folder, fn) for fn in all_files if fn.endswith('.pck') and 'rfcls' in fn]
    return model_files


def make_predictions(params):
    """
    :param params:
    :return:
    """
    model_file, data_file, output, lock = params
    model_name = os.path.basename(model_file)
    model_name = 'model_' + model_name.split('_')[-1].replace('.pck', '')
    with open(model_file, 'rb') as dump:
        cv_dump = pck.load(dump)
        model = cv_dump['cv'].best_estimator_
        model_weight = cv_dump['valid_score']
    chromosomes = set()
    with pd.HDFStore(data_file, 'r') as hdf_in:
        for chrom in hdf_in.keys():
            chromosomes.add(chrom.strip('/'))
            chrom_data = hdf_in[chrom]
            feat_columns = sorted([c for c in chrom_data.columns if c.startswith('ft')])
            feat_matrix = chrom_data[feat_columns]
            pred_label = model.predict(feat_matrix)
            df_pred_prob = pd.DataFrame(model.predict_proba(feat_matrix),
                                        columns=model.classes_)
            pred_prob = df_pred_prob.lookup(df_pred_prob.index, pred_label)
            chrom_data['pred_label_' + model_name] = pred_label
            chrom_data['chrom'] = chrom.strip('/')
            chrom_data['pred_prob_' + model_name] = pred_prob
            out_cols = [c for c in chrom_data.columns if not c.startswith('ft')]
            with lock:
                with pd.HDFStore(output, 'a', complevel=9, complib='blosc') as hdf_out:
                    hdf_out.put(os.path.join(chrom, model_name), chrom_data[out_cols], format='table')
                    hdf_out.flush()
    return model_name, model_weight, chromosomes


def merge_predictions(params):
    """
    :param params:
    :return:
    """
    threshold = 11 * 0.6
    chrom, output_file, lock = params
    estimates = None
    with pd.HDFStore(output_file, 'r') as hdf:
        load_keys = [k for k in hdf.keys() if k.startswith('/' + chrom + '/')]
        for k in load_keys:
            next_data = hdf[k]
            if estimates is None:
                estimates = next_data
            else:
                estimates = estimates.merge(next_data, how='outer',
                                            suffixes=('', ''), on=None)
    label_columns = sorted([c for c in estimates.columns if 'pred_label_' in c])
    labels = estimates[label_columns]
    prob_columns = sorted([c for c in estimates.columns if 'pred_prob_' in c])
    probs = estimates[prob_columns]
    weighted_est = np.multiply(labels.values, probs.values)
    wt = weighted_est.sum(axis=1)
    estimates['total_weight'] = wt
    estimates['avg_weight'] = wt / labels.shape[1]
    select_idx = wt > threshold
    estimates = estimates.loc[select_idx, :]
    enh_collect = col.defaultdict(list)
    for row in estimates.itertuples():
        enh_collect[(row.name, row.chrom)].append((row.enh_name, row.enh_start, row.enh_end,
                                                   row.avg_weight, row.total_weight))
    assoc = []
    for k, v in enh_collect.items():
        v = sorted(v, key=lambda x: x[1])
        starts = ','.join([str(x[1]) for x in v])
        ends = ','.join([str(x[2]) for x in v])
        avg = ','.join([str(round(x[3], 2)) for x in v])
        total = ','.join([str(round(x[4], 2)) for x in v])
        assoc.append({'name': k[0], 'chrom': k[1], 'starts': starts,
                     'ends': ends, 'weights': avg, 'total_weights': total})
    out_columns = ['name', 'chrom', 'starts', 'ends', 'weights', 'total_weights']
    df = pd.DataFrame.from_records(assoc, index=range(len(assoc)),
                                   columns=out_columns)
    with lock:
        with pd.HDFStore(output_file, 'a', complib='blosc', complevel=9) as hdf:
            hdf.put('ega/' + chrom, df, format='table')
    return chrom


def apply_models(args):
    """
    :param args:
    :return:
    """
    model_files = load_models(args.models)
    chromosomes = set()
    models = set()
    with mp.Manager() as mng:
        lock = mng.Lock()
        params = [(mf, args.data, args.output, lock) for mf in model_files]
        with mp.Pool(args.workers) as pool:
            resit = pool.imap_unordered(make_predictions, params)
            for model, weight, chroms in resit:
                chromosomes = chromosomes.union(chroms)
                models.add((model, weight))
        params = [(c, args.output, lock) for c in chromosomes]
        with mp.Pool(args.workers) as pool:
            resit = pool.imap_unordered(merge_predictions, params)
            for chrom in resit:
                pass
    return


def extract_associations(args):
    """
    :param args:
    :return:
    """
    with pd.HDFStore(args.data, 'r') as hdf:
        headered = False
        with open(args.output, 'w') as out:
            load_keys = [k for k in hdf.keys() if k.startswith('/ega')]
            for k in load_keys:
                data = hdf[k]
                if not headered:
                    header = ['name', 'chrom', 'starts', 'ends', 'weights', 'total_weights']
                    _ = out.write('\t'.join(header) + '\n')
                    headered = True
                data.to_csv(out, sep='\t', index=False, header=False,
                            columns=header)
    return


if __name__ == '__main__':
    try:
        args = parse_command_line()
        if args.tune:
            train_models(args)
        elif args.apply:
            apply_models(args)
        elif args.extract:
            extract_associations(args)
        else:
            raise RuntimeError('No mode selected: {}'.format(args))
    except Exception as err:
        trb.print_exc()
        raise err
    else:
        sys.exit(0)
