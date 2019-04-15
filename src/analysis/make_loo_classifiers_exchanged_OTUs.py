#!/usr/bin/env python
"""
This script makes many leave-one-out classifiers based on the exchanged
OTUs only.
"""

import argparse
import os
import sys
import copy

import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from scipy.stats import fisher_exact

from make_loo_classifiers import write_samples, make_combined_site_df, loo_classify

# add util/ to path
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
import util

def make_concordance_df(tidydf, mbs_col, site1, site2):
    """
    Make dataframe where the values are 0 and 1, indicating whether
    each OTU/column was concordant in both sites or not.
    """
    ## Get subjects with both sites
    tmp = tidydf\
        .query('(site == @site1) | (site == @site2)')\
        .dropna(subset=[mbs_col])\
        .groupby(['subject_id', 'otu'])\
        .size()
    tmp = tmp[tmp == 2].reset_index()

    subjs = tmp['subject_id'].unique().tolist()

    # Track the samples
    samples = tidydf\
        .query('subject_id == @subjs')\
        .query('(site == @site1) | (site == @site2)')\
        ['sample']\
        .unique().tolist()

    ## Make separate BAL and throat dataframes
    df1 = tidydf\
        .query('site == @site1')\
        .query('subject_id == @subjs')\
        .pivot(index='subject_id', columns='otu_w_site',
               values='abun')

    df2 = tidydf\
        .query('site == @site2')\
        .query('subject_id == @subjs')\
        .pivot(index='subject_id', columns='otu_w_site',
               values='abun')

    # Note: to adapt this code for classifiers which use more than just
    # the exchanged OTUs, you would need to remove any OTUs which are
    # zero in all patients in all sites. Right now, since we're specifically
    # looking at exchanged OTUs (which were defined a priori), we'll keep
    # all of them, even if no patients have any of them in any site.
    # Regardless, this shouldn't affect results except by introducing more
    # noise into the classifiers...

    ## Even though pandas should have done this already, double-make-sure
    ## that rows and columns match
    # Remove site label appended to OTUs
    df1.columns = [c.split('--')[0] for c in df1.columns]
    df2.columns = [c.split('--')[0] for c in df2.columns]
    # Re-order df1 so it matches df2
    df1 = df1.loc[df2.index, df2.columns]

    # Get concordance values: if both are greater than 0 OR both are 0,
    # returns 1. Otherwise 0
    concordance = ((df1 > 0) == (df2 > 0)).astype(int)

    return concordance, samples

def loo_classify_update_results(df, sites, mbs_info, cls_type,
                                random_state, PREDS, SUMMS):
    """
    Make classifier (using loo_classify()) and update lists
    of results with the classifier results and an indicator
    of which classifier was used.

    df : pandas DataFrame
    sites : list of str
    mbs_info : dict with {subject_id: 0 or 1}
    cls_type : str labeling classifier type (abundance or presence)
    random_state : int
    PREDS, SUMMS : lists of results to append these to
    """
    ## Make classifier
    resdf, roc_auc, fisher_p = loo_classify(df, mbs_info, random_state)

    ## Calculate N asp and N normal, for storing
    n_asp = sum([mbs_info[i] == 1 for i in df.index.tolist()])
    n_nml = sum([mbs_info[i] == 0 for i in df.index.tolist()])

    ## Store results
    resdf['site'] = '-'.join(sites)
    resdf['classifier_type'] = cls_type
    PREDS.append(resdf)

    SUMMS.append([
        '-'.join(sites), cls_type,
        n_asp, n_nml, df.shape[1],
        roc_auc, fisher_p])

    return PREDS, SUMMS

# def multi_site_classifier((df, sites, iteration, random_state, cls_label)):
#     """
#     Make a multi-site classifiers.
#
#     First makes a datafarme with has subject IDs in rows
#     and OTUs in columns. OTUs are labeled by which site
#     they were in (e.g. k__Bacteria;...;d__denovo123-bal and
#     k__Bacteria;...;d__denovo123-gastric_fluid are separate columns).
#
#     Global variables used:
#         meta (with column 'subject_id'), aspdict, mbs_col
#
#     Parameters
#     ----------
#     df : pandas DataFrame
#         subjects in index, OTUs (w site) in columns
#     sites : list
#         list of site(s) to consider
#     iter : int [default : np.nan]
#         iteration number, to track the ROC curve and predictions
#     cls_label : str
#         classifier label (e.g. 'presence' or 'abundance')
#
#     Returns
#     -------
#     summaries, rocs, predictions : pandas DataFrames
#
#     """
#     print(sites, cls_label, df.shape, iteration, random_state)
#     subjects = df.index.tolist()
#
#     # Convert data to scikit-learn format
#     # Using pandas.query matches order given in subjects, which matches
#     # order in df (which becomes X)
#     mbs_info = tidydf\
#         .query('subject_id == @subjects')\
#         [['subject_id', 'mbs_consolidated']]\
#         .drop_duplicates()
#
#     # Code up labels as 0 or 1
#     Y = [aspdict[i] for i in mbs_info['mbs_consolidated']]
#     X = df.values
#
#     rf = RandomForestClassifier(n_estimators=1000, random_state=random_state)
#
#     # Classify
#     results = util.cv_and_roc(rf, X, Y, random_state=random_state)
#
#     # Get ROC results
#     mean_tpr = np.mean(results['tpr_list'], axis=0)
#     rocs = pd.DataFrame(
#         columns=['mean_tpr', 'mean_fpr'],
#         data=np.array([mean_tpr, results['mean_fpr']]).T)
#     rocs['iteration'] = iteration
#     rocs['random_state'] = random_state
#     rocs['site'] = '-'.join(sites)
#     rocs['cls_type'] = cls_label
#     rocs['n_feats'] = X.shape[1]
#
#     # Get predictions
#     predictions = pd.DataFrame(
#         columns=['subject', 'true', 'predicted', 'probability'],
#         data=np.array((subjects, results['y_true'],
#                        results['y_predictions'], results['y_prob'])).T)
#     predictions['iteration'] = iteration
#     predictions['random_state'] = random_state
#     predictions['site'] = '-'.join(sites)
#     predictions['cls_type'] = cls_label
#     predictions['n_feats'] = X.shape[1]
#
#     # Get other single-valued summaries
#     summaries = pd.DataFrame(
#         {'auc': np.mean(results['auc_list']),
#         'fisher_p': results['fisher_p'],
#         'sensitivity': results['sensitivity'],
#         'specificity': results['specificity'],
#         'ppv': results['ppv'],
#         'npv': results['npv'],
#         'N_asp': sum(np.array(Y) == 1),
#         'N_nml': sum(np.array(Y) == 0),
#         'sites': '-'.join(sites),
#         'iteration': iteration,
#         'cls_type': cls_label,
#         'n_feats': X.shape[1]
#         },
#         index=[0])
#
#     return (summaries, rocs, predictions)
#
# def parallel_setup_and_classify(df, sites, cls_label,
#                                 SUMMS, ROCS, PREDS):
#     """
#     Set up the dataframe (once) with corresponding sites.
#     Classify and update results.
#
#     Parameters
#     ----------
#     df : pandas DataFrame
#     sites : list of strings
#     cls_label : str
#     SUMMS, ROCS, PREDS : lists of lists
#
#     Global variables: mbs_col, nreps,
#     plus all the things multi_site_classifier uses
#     """
#
#     # Set up for parallel processing
#     DFS = [df]*nreps
#     SITES = [sites]*nreps
#     ITERS = range(nreps)
#     RANDSTATES = np.random.randint(0, 20000, size=nreps)
#     CLS_LABELS = [cls_label]*nreps
#
#     p = multiprocessing.Pool()
#     results = p.map(multi_site_classifier,
#         zip(DFS, SITES, ITERS, RANDSTATES, CLS_LABELS))
#     p.close()
#     p.join()
#
#     # Update full results with these dataframes
#     # Note: results is a list of tuples, where the first element of each tuple
#     # is the summary dataframe, the second is the rocs dataframe, etc
#     SUMMS += [res[0] for res in results]
#     ROCS += [res[1] for res in results]
#     PREDS += [res[2] for res in results]
#
#     return (SUMMS, ROCS, PREDS)

p = argparse.ArgumentParser()
# Input files
p.add_argument('fnotu', help='Clean OTU table. Samples in rows, OTUs '
    + 'in columns.')
p.add_argument('fnmeta', help='Clean metadata. All samples (in rows) should'
    + ' be in the OTU table.')
p.add_argument('fnexchange', help='Path to file with exchanged OTUs. '
    + 'Classifiers are built based only on the exchanged OTUs. This file '
    + 'has the calculated prevalence values, and only includes the OTUs '
    + 'which have been determined to be exchanged. Important columns in '
    + 'this file are "site_comparison" and "otu".')
# Output files
p.add_argument('fnsummaries', help='File where to write classifier summary '
    + 'metrics (AUC, fisher p).')
#p.add_argument('fnrocs', help='File where to write classifier ROC data '
#    + '(TPR, FPR).')
p.add_argument('fnpreds', help='File where to write classifier predictions '
    + '(subjects, predictions, true metadata).')
args = p.parse_args()

print('Reading files... '),
df = pd.read_csv(args.fnotu, sep='\t', index_col=0)
meta = pd.read_csv(args.fnmeta, sep='\t', index_col=0)
print('Finished.')

# Define some global variables
nml = 'Normal'
asp = 'Aspiration/Penetration'
mbs_col = 'mbs_consolidated'
aspdict = {'Normal': 0, 'Aspiration/Penetration': 1}
sites = util.get_sites()
random_state = 12345

# Bad form but oh well: output file for patient lists
out_patients = 'data/patients/aspiration_classifiers_exchanged_OTUs'

## Prepare the data for easy manipulations
# This creates a 3-column df with ['sample', 'otu_w_site', 'abun']
meta = meta.dropna(subset=[mbs_col])
# All-caps tidydf doesn't change
TIDYDF = util.tidyfy_otu(df, meta, mbs_col)

## Make the aspiration status mapping dictionary
# Only need to do this once, since it has every
# subject's aspiration status
mbs_info = TIDYDF\
    [['subject_id', 'mbs_consolidated']]\
    .drop_duplicates()
# Convert to dictionary
mbs_info = dict(zip(mbs_info['subject_id'], mbs_info['mbs_consolidated']))
# Convert to 1 or 0
mbs_info = {k: aspdict[mbs_info[k]] for k in mbs_info}

## Since we're doing the exchanged OTUs only, prepare that dataframe too
exch = pd.read_csv(args.fnexchange, sep='\t')
# Remove all extra info, like prevalence
exch = exch[["otu", "site_comparison"]].drop_duplicates()

## Initialize list of results
SUMMS = []
#ROCS = []
PREDS = []

##########################################
## Classifiers based on BAL-throat OTUs ##
##########################################
keepotus = exch.query('site_comparison == "bal-throat_swab"')['otu'].tolist()
tidydf = TIDYDF.query("otu == @keepotus")

## Abundance-based classifiers
cls_label = 'abundance'
sites = ['bal']
print(cls_label, '-'.join(sites))
# Make wide dataframe with the right samples and OTUs
df, samples = make_combined_site_df(tidydf, sites, mbs_col)
# Track samples used in the classifier
write_samples(sites, out_patients, samples)
## Make classifier and update the PREDS and SUMMS lists
PREDS, SUMMS = loo_classify_update_results(
    df, sites, mbs_info, cls_label, random_state, PREDS, SUMMS)

sites = ['throat_swab']
print(cls_label, '-'.join(sites))
# Make wide dataframe with the right samples and OTUs
df, samples = make_combined_site_df(tidydf, sites, mbs_col)
# Track samples used in the classifier
write_samples(sites, out_patients, samples)
## Make classifier and update the PREDS and SUMMS lists
PREDS, SUMMS = loo_classify_update_results(
    df, sites, mbs_info, cls_label, random_state, PREDS, SUMMS)

sites = ['bal', 'throat_swab']
print(cls_label, '-'.join(sites))
# Make wide dataframe with the right samples and OTUs
df, samples = make_combined_site_df(tidydf, sites, mbs_col)
# Track samples used in the classifier
write_samples(sites, out_patients, samples)
## Make classifier and update the PREDS and SUMMS lists
PREDS, SUMMS = loo_classify_update_results(
    df, sites, mbs_info, cls_label, random_state, PREDS, SUMMS)

## Presence-based classifiers
# No need to re-track samples, since they're the same ones as above
tidydf.loc[:, 'abun'] = (tidydf['abun'] > 0).astype(int)
cls_label = 'presence'

sites = ['bal']
print(cls_label, '-'.join(sites))
df, _ = make_combined_site_df(tidydf, sites, mbs_col)
## Make classifier and update the PREDS and SUMMS lists
PREDS, SUMMS = loo_classify_update_results(
    df, sites, mbs_info, cls_label, random_state, PREDS, SUMMS)

sites = ['throat_swab']
print(cls_label, '-'.join(sites))
df, _ = make_combined_site_df(tidydf, sites, mbs_col)
## Make classifier and update the PREDS and SUMMS lists
PREDS, SUMMS = loo_classify_update_results(
    df, sites, mbs_info, cls_label, random_state, PREDS, SUMMS)

sites = ['bal', 'throat_swab']
print(cls_label, '-'.join(sites))
# With the two-site presence/absence classifier, look at the concordance
df, _ = make_concordance_df(tidydf, mbs_col, sites[0], sites[1])
## Make classifier and update the PREDS and SUMMS lists
PREDS, SUMMS = loo_classify_update_results(
    df, sites, mbs_info, cls_label, random_state, PREDS, SUMMS)


###########################################
## Classifiers based on BAL-gastric OTUs ##
###########################################
keepotus = exch.query('site_comparison == "bal-gastric_fluid"')['otu'].tolist()
tidydf = TIDYDF.query("otu == @keepotus")

# Abundance-based classifiers
cls_label = 'abundance'
sites = ['bal']
print(cls_label, '-'.join(sites))
# Make wide dataframe with the right samples and OTUs
df, samples = make_combined_site_df(tidydf, sites, mbs_col)
# Track samples used in the classifier
write_samples(sites, out_patients, samples)
## Make classifier and update the PREDS and SUMMS lists
PREDS, SUMMS = loo_classify_update_results(
    df, sites, mbs_info, cls_label, random_state, PREDS, SUMMS)

sites = ['gastric_fluid']
print(cls_label, '-'.join(sites))
# Make wide dataframe with the right samples and OTUs
df, samples = make_combined_site_df(tidydf, sites, mbs_col)
# Track samples used in the classifier
write_samples(sites, out_patients, samples)
## Make classifier and update the PREDS and SUMMS lists
PREDS, SUMMS = loo_classify_update_results(
    df, sites, mbs_info, cls_label, random_state, PREDS, SUMMS)

sites = ['bal', 'gastric_fluid']
print(cls_label, '-'.join(sites))
# Make wide dataframe with the right samples and OTUs
df, samples = make_combined_site_df(tidydf, sites, mbs_col)
# Track samples used in the classifier
write_samples(sites, out_patients, samples)
## Make classifier and update the PREDS and SUMMS lists
PREDS, SUMMS = loo_classify_update_results(
    df, sites, mbs_info, cls_label, random_state, PREDS, SUMMS)

# Presence-based classifiers
tidydf.loc[:, 'abun'] = (tidydf['abun'] > 0).astype(int)
cls_label = 'presence'

sites = ['bal']
print(cls_label, '-'.join(sites))
df, _ = make_combined_site_df(tidydf, sites, mbs_col)
## Make classifier and update the PREDS and SUMMS lists
PREDS, SUMMS = loo_classify_update_results(
    df, sites, mbs_info, cls_label, random_state, PREDS, SUMMS)

sites = ['gastric_fluid']
print(cls_label, '-'.join(sites))
df, _ = make_combined_site_df(tidydf, sites, mbs_col)
## Make classifier and update the PREDS and SUMMS lists
PREDS, SUMMS = loo_classify_update_results(
    df, sites, mbs_info, cls_label, random_state, PREDS, SUMMS)

sites = ['bal', 'gastric_fluid']
df, _ = make_concordance_df(tidydf, mbs_col, sites[0], sites[1])
PREDS, SUMMS = loo_classify_update_results(
    df, sites, mbs_info, cls_label, random_state, PREDS, SUMMS)

## Concatenate results dataframes and write
# SUMMS is a list of lists, turn into a dataframe
pd.DataFrame(SUMMS,
    columns=['site', 'classifier_type',
             'n_asp', 'n_nml', 'n_feats',
             'auc', 'fisher_p']
    ).to_csv(args.fnsummaries, sep='\t', index=False)
# PREDS is a list of dataframes, just concat
pd.concat(PREDS).to_csv(args.fnpreds, sep='\t', index=False)
