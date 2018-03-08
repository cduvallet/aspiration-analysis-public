#!/usr/bin/env python
"""
This script makes many different types of classifiers and stores all of
its respective outputs.
"""

import argparse
import os
import sys
import multiprocessing

import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from scipy.stats import fisher_exact

# add util/ to path
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
import util


# def one_rep_rf(tmpotu, aspdict, random_state, n_cv, shuffle=False):
#     """
#     Do one replicate of the Random Forest workflow.
#     First, build a n_cv-fold cross-validated random forest (seeded by
#     random_state) to get the ROC curve and predictions.
#     Then, re-build a random forest based on all of the data to get
#     the feature importances and the oob score.
#     Return dataframes that will eventually be concatenated into tidy
#     dataframes.
#
#     Parameters
#     ----------
#     tmpotu : pandas DataFrame
#         OTUs in columns, samples in rows.
#     aspdict : dict
#         Dictionary mapping samples from tmpotu to 0/1 aspiration status
#     random_state : int
#         Seed for all random forests.
#     n_cv : int
#         Number of cross-validation folds.
#     shuffle : bool
#         whether to shuffle the Y labels before building classifier
#
#     Returns
#     -------
#     summaries : pandas DataFrame
#         cross-validated AUC, cross-validated fisher P, all data oob score
#     rocs : pandas DataFrame
#         'mean_tpr' and 'mean_fpr', interpolated ROC values from cross validation
#     predictions : pandas DataFrame
#         subject, true label, predicted label, probability. From cross-validation
#     features : pandas DataFrame
#         feature, importance. From non-cross-validated RF.
#     """
#
#     X = tmpotu.values
#     Y = [aspdict[i] for i in tmpotu.index]
#     if shuffle:
#         Y = np.random.permutation(Y)
#
#     rf = RandomForestClassifier(n_estimators=1000,
#                                 random_state=random_state)
#     ## Cross-validation results for summary stats (except oob), ROC,
#     ## and predictions
#     results = cv_and_roc(rf, X, Y, random_state=random_state, num_cv=n_cv)
#
#     # Update ROC info
#     rocs = pd.DataFrame(columns=['mean_tpr', 'mean_fpr'],
#                         data=np.array([results['mean_tpr'],
#                                        results['mean_fpr']]).T)
#     # Update predictions
#     predictions = pd.DataFrame(columns=['subject', 'true', 'predicted',
#                                         'probability'],
#                                data=np.array((tmpotu.index.tolist(),
#                                               results['y_true'],
#                                               results['y_predictions'],
#                                               results['y_prob'])).T)
#
#     ## Classifier trained on all data for oob_score and features
#     rf = RandomForestClassifier(n_estimators=1000,
#                                 random_state=random_state,
#                                 oob_score=True).fit(X, Y)
#     # Update features
#     features = pd.DataFrame(columns=['feature', 'importance'],
#                             data=np.array((tmpotu.columns.tolist(),
#                                            rf.feature_importances_)).T)
#
#     # Update classifier summary stats
#     summaries = pd.DataFrame({'auc': results['roc_auc'],
#                               'fisher_p': results['fisher_p'],
#                               'oob_score': rf.oob_score_,
#                               'sensitivity': results['sensitivity'],
#                               'specificity': results['specificity'],
#                               'ppv': results['ppv'],
#                               'npv': results['npv'],
#                               'N_asp': sum(np.array(Y) == 1),
#                               'N_nml': sum(np.array(Y) == 0)},
#                              index=[0])
#
#     return (summaries, rocs, predictions, features)

def make_combined_site_df(tidydf, sites, mbs_col):
    """
    Return a wide-form dataframe with data from all sites.

    Drops any rows with NaNs (i.e. subjects which are missing
    one of the sites). Also drops any subjects without mbs_col
    metadata.

    Parameters
    ----------
    tidydf : pandas DataFrame
        'subject_id', 'site', and 'otu_w_site' columns
    sites : list
        list of sites to keep
    """

    tmpotu = tidydf.query('site == @sites')\
            .dropna(subset=[mbs_col])\
            .pivot(index='subject_id', columns='otu_w_site',
                   values='abun')\
            .dropna(axis=0)

    return tmpotu

def multi_site_classifier((df, sites, iteration, random_state)):
    """
    Make a multi-site classifiers.

    First makes a datafarme with has subject IDs in rows
    and OTUs in columns. OTUs are labeled by which site
    they were in (e.g. k__Bacteria;...;d__denovo123-bal and
    k__Bacteria;...;d__denovo123-gastric_fluid are separate columns).

    Global variables used:
        meta (with column 'subject_id'), aspdict, mbs_col

    Parameters
    ----------
    df : pandas DataFrame
        subjects in index, OTUs (w site) in columns

    sites : list
        list of site(s) to consider

    iter : int [default : np.nan]
        iteration number, to track the ROC curve and predictions

    Returns
    -------
    summaries, rocs, predictions : pandas DataFrames

    """
    print(sites, df.shape, iteration, random_state)
    subjects = df.index.tolist()

    # Convert data to scikit-learn format
    # Using pandas.query matches order given in subjects, which matches
    # order in df (which becomes X)
    mbs_info = tidydf\
        .query('subject_id == @subjects')\
        [['subject_id', 'mbs_consolidated']]\
        .drop_duplicates()

    # Code up labels as 0 or 1
    Y = [aspdict[i] for i in mbs_info['mbs_consolidated']]
    X = df.values

    rf = RandomForestClassifier(n_estimators=1000, random_state=random_state)

    # Classify
    results = util.cv_and_roc(rf, X, Y, random_state=random_state)

    # Get ROC results
    mean_tpr = np.mean(results['tpr_list'], axis=0)
    rocs = pd.DataFrame(
        columns=['mean_tpr', 'mean_fpr'],
        data=np.array([mean_tpr, results['mean_fpr']]).T)
    rocs['iteration'] = iteration
    rocs['random_state'] = random_state
    rocs['site'] = '-'.join(sites)

    # Get predictions
    predictions = pd.DataFrame(
        columns=['subject', 'true', 'predicted', 'probability'],
        data=np.array((subjects, results['y_true'],
                       results['y_predictions'], results['y_prob'])).T)
    predictions['iteration'] = iteration
    predictions['random_state'] = random_state
    predictions['site'] = '-'.join(sites)

    # Get other single-valued summaries
    summaries = pd.DataFrame(
        {'auc': np.mean(results['auc_list']),
        'fisher_p': results['fisher_p'],
        'sensitivity': results['sensitivity'],
        'specificity': results['specificity'],
        'ppv': results['ppv'],
        'npv': results['npv'],
        'N_asp': sum(np.array(Y) == 1),
        'N_nml': sum(np.array(Y) == 0),
        'sites': '-'.join(sites),
        'iteration': iteration
        },
        index=[0])

    return (summaries, rocs, predictions)

def parallel_setup_and_classify(sites, SUMMS, ROCS, PREDS):
    """
    Set up the dataframe (once) with corresponding sites.
    Classify and update results.

    Global variables: tidydf, mbs_col, nreps,
    plus all the things multi_site_classifier uses
    """
    # Convert tidy dataframe to wide form dataframe
    df = make_combined_site_df(tidydf, sites, mbs_col)

    # Set up for parallel processing
    DFS = [df]*nreps
    SITES = [sites]*nreps
    ITERS = range(nreps)
    RANDSTATES = np.random.randint(0, 20000, size=nreps)

    p = multiprocessing.Pool()
    results = p.map(multi_site_classifier, zip(DFS, SITES, ITERS, RANDSTATES))
    p.close()
    p.join()

    # Update full results with these dataframes
    # Note: results is a list of tuples, where the first element of each tuple
    # is the summary dataframe, the second is the rocs dataframe, etc
    SUMMS += [res[0] for res in results]
    ROCS += [res[1] for res in results]
    PREDS += [res[2] for res in results]

    return (SUMMS, ROCS, PREDS)

p = argparse.ArgumentParser()
# Input files
p.add_argument('fnotu', help='Clean OTU table. Samples in rows, OTUs '
    + 'in columns.')
p.add_argument('fnmeta', help='Clean metadata. All samples (in rows) should'
    + ' be in the OTU table.')
# Output files
p.add_argument('fnsummaries', help='File where to write classifier summary '
    + 'metrics (AUC, fisher p).')
p.add_argument('fnrocs', help='File where to write classifier ROC data '
    + '(TPR, FPR).')
p.add_argument('fnpreds', help='File where to write classifier predictions '
    + '(subjects, predictions, true metadata).')
# Optional parameters
p.add_argument('--rfreps', help='Number of Random Forest replicates '
    + '[default: %(default)s]', default=100, type=int)
args = p.parse_args()

print('Reading files... '),
df = pd.read_csv(args.fnotu, sep='\t', index_col=0)
meta = pd.read_csv(args.fnmeta, sep='\t', index_col=0)
print('Finished.')

nreps = args.rfreps

# Define some global variables
nml = 'Normal'
asp = 'Aspiration/Penetration'
mbs_col = 'mbs_consolidated'
aspdict = {'Normal': 0, 'Aspiration/Penetration': 1}
sites = util.get_sites()

## Prepare the data for easy manipulations
# This creates a 3-column df with ['sample', 'otu', 'abun']
meta = meta.dropna(subset=[mbs_col])
tidydf = util.tidyfy_otu(df, meta, mbs_col)

## Initialize list of results
SUMMS = []
ROCS = []
PREDS = []

## Single-site classifiers
sites = ['bal']
SUMMS, ROCS, PREDS = parallel_setup_and_classify(sites, SUMMS, ROCS, PREDS)

sites = ['throat_swab']
SUMMS, ROCS, PREDS = parallel_setup_and_classify(sites, SUMMS, ROCS, PREDS)

sites = ['gastric_fluid']
SUMMS, ROCS, PREDS = parallel_setup_and_classify(sites, SUMMS, ROCS, PREDS)

## Two-site classifiers
sites = ['bal', 'throat_swab']
SUMMS, ROCS, PREDS = parallel_setup_and_classify(sites, SUMMS, ROCS, PREDS)

sites = ['bal', 'gastric_fluid']
SUMMS, ROCS, PREDS = parallel_setup_and_classify(sites, SUMMS, ROCS, PREDS)

sites = ['throat_swab', 'gastric_fluid']
SUMMS, ROCS, PREDS = parallel_setup_and_classify(sites, SUMMS, ROCS, PREDS)

## All three
sites = ['bal', 'throat_swab', 'gastric_fluid']
SUMMS, ROCS, PREDS = parallel_setup_and_classify(sites, SUMMS, ROCS, PREDS)

## Concatenate results dataframes and write
pd.concat(SUMMS).to_csv(args.fnsummaries, sep='\t', index=False)
pd.concat(ROCS).to_csv(args.fnrocs, sep='\t', index=False)
pd.concat(PREDS).to_csv(args.fnpreds, sep='\t', index=False)
#
# [full_results[0].append(lst[0]) for lst in results]
#
# ## Two-site classifiers
# for site1 in sites:
#     for site2 in sites[sites.index(site1)+1:]:
#         print(site1 + ' and ' + site2)
#         for n in range(nreps):
#             print(n)
#             random_state = np.random.randint(0, 20000)
#             tmpotu = tidyotu.query('(site == @site1) | (site == @site2)')\
#                         .dropna(subset=[mbscol])\
#                         .pivot(index='subject_id', columns='otu_w_site',
#                                values='abun')\
#                         .dropna()
#             # tmpres is: summaries, rocs, predictions, features
#             tmpres = one_rep_rf(tmpotu, aspdict, random_state, n_cv)
#
#             for df, lst in zip(tmpres, results):
#                 df['site'] = site1 + '_and_' + site2
#                 df['replicate'] = n
#                 df['random_state'] = random_state
#                 df['shuffled'] = False
#                 lst.append(df)
#
#             if shuffle:
#                 print(str(n) + '.1')
#                 # Re-do the same thing but shuffle the labels
#                 # tmpres is: summaries, rocs, predictions, features
#                 tmpres = one_rep_rf(tmpotu, aspdict, random_state, n_cv,
#                                     shuffle)
#
#                 for df, lst in zip(tmpres, results):
#                     df['site'] = site1 + '_and_' + site2
#                     df['replicate'] = n
#                     df['random_state'] = random_state
#                     df['shuffled'] = shuffle
#                     lst.append(df)
#
#
# ## Three-site classifiers
# print('All sites')
# for n in range(nreps):
#     print(n)
#     random_state = np.random.randint(0, 20000)
#     tmpotu = tidyotu.query('site == @sites')\
#                     .dropna(subset=[mbscol])\
#                     .pivot(index='subject_id', columns='otu_w_site',
#                            values='abun')\
#                     .dropna()
#     # tmpres is: summaries, rocs, predictions, features
#     tmpres = one_rep_rf(tmpotu, aspdict, random_state, n_cv)
#
#     for df, lst in zip(tmpres, results):
#         df['site'] = 'all_sites'
#         df['replicate'] = n
#         df['random_state'] = random_state
#         df['shuffled'] = False
#         lst.append(df)
#
#     if shuffle:
#         # Re-do the same thing but shuffle the labels
#         # tmpres is: summaries, rocs, predictions, features
#         print(str(n) + '.1')
#         tmpres = one_rep_rf(tmpotu, aspdict, random_state, n_cv, shuffle)
#
#         for df, lst in zip(tmpres, results):
#             df['site'] = 'all_sites'
#             df['replicate'] = n
#             df['random_state'] = random_state
#             df['shuffled'] = shuffle
#             lst.append(df)
#
# # results list has: summaries, rocs, predictions, features
#
# ## Convert all lists into tidy dataframes and write to files
# pd.concat(results[0]).to_csv(args.fnsummaries, sep='\t', index=False)
# pd.concat(results[1]).to_csv(args.fnrocs, sep='\t', index=False)
# # Update the predictions to include other patient metadata
# results[2] = pd.concat(results[2])
# results[2] = results[2].merge(meta[metacols],
#                               left_on='subject',
#                               right_on='subject_id')
# results[2].to_csv(args.fnpreds, sep='\t', index=False)
# pd.concat(results[3]).to_csv(args.fnfeatures, sep='\t', index=False)
