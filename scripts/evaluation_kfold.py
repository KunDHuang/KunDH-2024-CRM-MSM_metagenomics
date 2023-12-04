#!/usr/bin/env python


"""
NAME: evaluation_kfold.py
DESCRIPTION: This script is to evaluate a model's prediction power using tenfold method.
DATE: 09.08.2022
AUTHOR: Kun D. Huang
"""

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from collections import namedtuple
from sklearn.calibration import CalibratedClassifierCV
from sklearn.feature_selection import SelectFromModel
import pandas as pd
import itertools 
import sys
import argparse
import textwrap
import subprocess
import random
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from sklearn.metrics import RocCurveDisplay
import matplotlib

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = 'Arial'

def read_args(args):
    # This function is to parse arguments

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                    description = textwrap.dedent('''\
                                    This script is to estimate ROC AUC based on metaphlan-style table with metadata being inserted.
                                    '''),
                                    epilog = textwrap.dedent('''\
                                    examples: evaluation_kfold.py --mpa_df <mpa_df.tsv> --md_rows 0,1,2,3,4 --target_row 3 --pos_feature <CRC> --neg_feature <Healthy> --fold_number 10 --repeat_time 20 --output ROC_AUC.svg  
                                    '''))

    parser.add_argument('--mpa_df',
                        nargs = '?',
                        help = 'Input a mpa-style table with metadata being inserted.',
                        type = str,
                        default = None)

    parser.add_argument('--md_rows',
                        nargs = '?',
                        help = 'Input row numbers for specifying metadata without considering header row, zero-based, comma delimited. for example, 0,1,2,3,4.',
                        default = None)

    parser.add_argument('--target_row',
                        nargs = '?',
                        help = 'Specify the row number for indicating target metadata to examine, zero-based without considering header row.',
                        type = int,
                        default = None)

    parser.add_argument('--pos_feature',
                        nargs = '?',
                        help = 'Specify the feature name to be labeled as posive, e.g. 1.',
                        type = str,
                        default = None)

    parser.add_argument('--neg_feature',
                        nargs = '?',
                        help = 'Specify the feature name to be labeld as negative, e.g. 0.',
                        type = str,
                        default = None)

    parser.add_argument('--fold_number',
                        nargs = '?',
                        help = 'Specify the fold number you want split the whole dataset.',
                        type = int,
                        default = None)

    parser.add_argument('--repeat_time',
                        nargs = '?',
                        help = 'Specify the repeat time you want to split the dataset.',
                        type = int,
                        default = None)

    parser.add_argument('--output',
                        nargs = '?',
                        help = 'Specify the output figure name.',
                        type = str,
                        default = None)

    parser.add_argument('--output_values',
                        nargs = '?',
                        help = 'Specify the output file name for storing ROC-AUC values.',
                        type = str,
                        default = None)

    parser.add_argument('--nproc',
                        nargs = '?',
                        help = 'Specify the number of processors you want to use. 4 by default.',
                        type = int,
                        default = 4)

    parser.add_argument('--transform',
                        nargs = '?',
                        help = 'Transform values in the matrix, [arcsin_sqrt] or [binary] or [None]. [None] by default',
                        type = str,
                        default = None)

    return vars(parser.parse_args())



def get_df_dropping_metadata(mpa4_df, row_number_list):
    # row_number_list: a list of row numbers in integer.
    # this function is to drop rows containing metadata.

    df_ = df_.drop(row_number_list)

    return df_

def get_target_metadata(mpa4_df, row_number):
    # row_number: the integer indicating the row which contains the metadata one wants to examine.
    # this function is to get a list of binary metadata for establishing ML model.
    
    features = mpa4_df.iloc[row_number].to_list()

    return features

def prepare_dataset(mpa4_style_md_df, pos_neg_dict, row_number_list, target_row, transform):
    # mpa4_style_md_df: the merged metaphlan4 table with metadata being inseted.
    # pos_neg_dict: the dictionary which maps examine value to 1 or 0.
    # This function is to prepare dataset for downstream analysis.

    df_no_md = mpa4_style_md_df.drop(row_number_list)
    sample_names = df_no_md.columns[1:].to_list()
    matrix = []
    for s in sample_names:
        values = [float(i) for i in df_no_md[s].to_list()]
        matrix.append(values)

    matrix = np.array(matrix)
    
    if transform == 'arcsin_sqrt':
        matrix = matrix/100
        print(matrix)
        matrix = np.arcsin(np.sqrt(matrix))

    elif transform == 'binary':
        matrix[np.where(matrix > 0)] = 1
    else:
        matrix = matrix

    features = get_target_metadata(mpa4_style_md_df, target_row)[1:]

    X = np.asarray(matrix)
    y = [pos_neg_dict[i] for i in features]
    y = np.asarray(y)

    return X, y


def roc_auc_curve(model, X, y, fold, repeat, output_name, output_values):
    # model: machine learning model to use.
    # X: the value matrix
    # y: the list of features, 1 and 0
    # fold: the number of fold to split the dataset.
    # repeat: the repeat number for splitting the dataset.
    # output_name: specify the output figure name.
    # outout_values: specify the output file name for storing estimated roc-auc values.

    cv = StratifiedKFold(n_splits=fold, shuffle=True)
    classifier = model

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    
    rocauc_opt = open(output_values, "w")

    rocauc_opt.write("repeat"+ "\t" + "fold" + "\t" + "roc_auc" + "\n")
    while repeat > 0:
        repeat -= 1
        for i, (train, test) in enumerate(cv.split(X, y)):
            classifier.fit(X[train], y[train])
            viz = RocCurveDisplay.from_estimator(
                classifier,
                X[test],
                y[test],
                name="ROC fold {}".format(i),
                alpha=0.3,
                lw=1,
            )
            interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)
            aucs.append(viz.roc_auc)
            rocauc_opt.write(str(repeat) + "\t" + str(i) + "\t" + str(viz.roc_auc) + "\n")

    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(
        mean_fpr,
        mean_tpr,
        color="b",
        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
        lw=2,
        alpha=0.8,
    )
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(
    mean_fpr,
    tprs_lower,
    tprs_upper,
    color="grey",
    alpha=0.2,
    label=r"$\pm$ 1 std. dev.",
    )
    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
    )
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.legend(loc="lower right")
    plt.savefig(output_name)
    rocauc_opt.close()



if __name__ == "__main__":

    pars = read_args(sys.argv)
    df = pd.read_csv(pars["mpa_df"], sep = "\t", index_col = False)
    row_number_list = [int(i) for i in pars["md_rows"].split(",")]
    pos_neg_dict = {pars["pos_feature"]:1, pars["neg_feature"]:0}
    X, y = prepare_dataset(df, pos_neg_dict, row_number_list, pars["target_row"], pars["transform"])
    model = RandomForestClassifier(n_estimators = 1000,
                                   criterion = 'entropy',
                                   min_samples_leaf = 1,
                                   max_features = 'sqrt',
                                   n_jobs = 4) # initiating a RF classifier
    roc_auc_curve(model, X, y, pars["fold_number"], pars["repeat_time"], pars["output"], pars["output_values"])
