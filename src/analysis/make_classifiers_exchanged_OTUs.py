#!/usr/bin/env python
"""
This script makes many different types of classifiers and stores all of
its respective outputs.
"""

import argparse
import os
import sys
import multiprocessing
import copy

import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from scipy.stats import fisher_exact

# add util/ to path
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
import util

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

    # Keep only OTUs which are non-zero in these samples
    tmpotu = tmpotu.loc[:, tmpotu.sum() > 0]

    # Also return the samples used in this
    subjects = tmpotu.index.tolist()
    samples = tidydf\
        .query('subject_id == @subjects')\
        .query('site == @sites')\
        ['sample']\
        .unique().tolist()

    return tmpotu, samples

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

def multi_site_classifier((df, sites, iteration, random_state, cls_label)):
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
    cls_label : str
        classifier label (e.g. 'presence' or 'abundance')

    Returns
    -------
    summaries, rocs, predictions : pandas DataFrames

    """
    print(sites, cls_label, df.shape, iteration, random_state)
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
    rocs['cls_type'] = cls_label
    rocs['n_feats'] = X.shape[1]

    # Get predictions
    predictions = pd.DataFrame(
        columns=['subject', 'true', 'predicted', 'probability'],
        data=np.array((subjects, results['y_true'],
                       results['y_predictions'], results['y_prob'])).T)
    predictions['iteration'] = iteration
    predictions['random_state'] = random_state
    predictions['site'] = '-'.join(sites)
    predictions['cls_type'] = cls_label
    predictions['n_feats'] = X.shape[1]

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
        'iteration': iteration,
        'cls_type': cls_label,
        'n_feats': X.shape[1]
        },
        index=[0])

    return (summaries, rocs, predictions)

def parallel_setup_and_classify(df, sites, cls_label,
                                SUMMS, ROCS, PREDS):
    """
    Set up the dataframe (once) with corresponding sites.
    Classify and update results.

    Parameters
    ----------
    df : pandas DataFrame
    sites : list of strings
    cls_label : str
    SUMMS, ROCS, PREDS : lists of lists

    Global variables: mbs_col, nreps,
    plus all the things multi_site_classifier uses
    """

    # Set up for parallel processing
    DFS = [df]*nreps
    SITES = [sites]*nreps
    ITERS = range(nreps)
    RANDSTATES = np.random.randint(0, 20000, size=nreps)
    CLS_LABELS = [cls_label]*nreps

    p = multiprocessing.Pool()
    results = p.map(multi_site_classifier,
        zip(DFS, SITES, ITERS, RANDSTATES, CLS_LABELS))
    p.close()
    p.join()

    # Update full results with these dataframes
    # Note: results is a list of tuples, where the first element of each tuple
    # is the summary dataframe, the second is the rocs dataframe, etc
    SUMMS += [res[0] for res in results]
    ROCS += [res[1] for res in results]
    PREDS += [res[2] for res in results]

    return (SUMMS, ROCS, PREDS)

def write_samples(sites, fout_patients, samples):
    """
    Write the list of samples to a file.

    sites : list of str
    fout_patients : str, with stem of file name
    samples : list of str, sample IDs
    """
    fname = '.'.join([fout_patients, '_'.join(sites), 'samples', 'txt'])
    with open(fname, 'w') as f:
        f.write('\n'.join(samples))

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

# Bad form but oh well: output file for patient lists
out_patients = 'data/patients/aspiration_classifiers_exchanged_OTUs'

## Prepare the data for easy manipulations
# This creates a 3-column df with ['sample', 'otu_w_site', 'abun']
meta = meta.dropna(subset=[mbs_col])
# All-caps tidydf doesn't change
TIDYDF = util.tidyfy_otu(df, meta, mbs_col)

## Since we're doing the exchanged OTUs only, prepare that dataframe too
exch = pd.read_csv(args.fnexchange, sep='\t')
# Remove all extra info, like prevalence
exch = exch[["otu", "site_comparison"]].drop_duplicates()

## Initialize list of results
SUMMS = []
ROCS = []
PREDS = []

##########################################
## Classifiers based on BAL-throat OTUs ##
##########################################
keepotus = exch.query('site_comparison == "bal-throat_swab"')['otu'].tolist()
tidydf = TIDYDF.query("otu == @keepotus")

## Abundance-based classifiers
cls_label = 'abundance'
sites = ['bal']

# Make wide dataframe with the right samples and OTUs
df, samples = make_combined_site_df(tidydf, sites, mbs_col)
# Track samples used in the classifier
write_samples(sites, out_patients, samples)
SUMMS, ROCS, PREDS = parallel_setup_and_classify(
    df, sites, cls_label, SUMMS, ROCS, PREDS)

sites = ['throat_swab']
df, samples = make_combined_site_df(tidydf, sites, mbs_col)
write_samples(sites, out_patients, samples)
SUMMS, ROCS, PREDS = parallel_setup_and_classify(
    df, sites, cls_label, SUMMS, ROCS, PREDS)

sites = ['bal', 'throat_swab']
df, samples = make_combined_site_df(tidydf, sites, mbs_col)
write_samples(sites, out_patients, samples)
SUMMS, ROCS, PREDS = parallel_setup_and_classify(
    df, sites, cls_label, SUMMS, ROCS, PREDS)

## Presence-based classifiers
# No need to re-track samples, since they're the same ones as above
tidydf.loc[:, 'abun'] = (tidydf['abun'] > 0).astype(int)
cls_label = 'presence'

sites = ['bal']
df, _ = make_combined_site_df(tidydf, sites, mbs_col)
SUMMS, ROCS, PREDS = parallel_setup_and_classify(
    df, sites, cls_label, SUMMS, ROCS, PREDS)

sites = ['throat_swab']
df, _ = make_combined_site_df(tidydf, sites, mbs_col)
SUMMS, ROCS, PREDS = parallel_setup_and_classify(
    df, sites, cls_label, SUMMS, ROCS, PREDS)

sites = ['bal', 'throat_swab']
# With the two-site presence/absence classifier, look at the concordance
df, _ = make_concordance_df(tidydf, mbs_col, sites[0], sites[1])
SUMMS, ROCS, PREDS = parallel_setup_and_classify(
    df, sites, cls_label, SUMMS, ROCS, PREDS)


###########################################
## Classifiers based on BAL-gastric OTUs ##
###########################################
keepotus = exch.query('site_comparison == "bal-gastric_fluid"')['otu'].tolist()
tidydf = TIDYDF.query("otu == @keepotus")

# Abundance-based classifiers
cls_label = 'abundance'
sites = ['bal']
df, samples = make_combined_site_df(tidydf, sites, mbs_col)
write_samples(sites, out_patients, samples)
SUMMS, ROCS, PREDS = parallel_setup_and_classify(
    df, sites, cls_label, SUMMS, ROCS, PREDS)

sites = ['gastric_fluid']
df, samples = make_combined_site_df(tidydf, sites, mbs_col)
write_samples(sites, out_patients, samples)
SUMMS, ROCS, PREDS = parallel_setup_and_classify(
    df, sites, cls_label, SUMMS, ROCS, PREDS)

sites = ['bal', 'gastric_fluid']
df, samples = make_combined_site_df(tidydf, sites, mbs_col)
write_samples(sites, out_patients, samples)
SUMMS, ROCS, PREDS = parallel_setup_and_classify(
    df, sites, cls_label, SUMMS, ROCS, PREDS)

# Presence-based classifiers
tidydf.loc[:, 'abun'] = (tidydf['abun'] > 0).astype(int)
cls_label = 'presence'

sites = ['bal']
df, _ = make_combined_site_df(tidydf, sites, mbs_col)
SUMMS, ROCS, PREDS = parallel_setup_and_classify(
    df, sites, cls_label, SUMMS, ROCS, PREDS)

sites = ['gastric_fluid']
df, _ = make_combined_site_df(tidydf, sites, mbs_col)
SUMMS, ROCS, PREDS = parallel_setup_and_classify(
    df, sites, cls_label, SUMMS, ROCS, PREDS)

sites = ['bal', 'gastric_fluid']
df, _ = make_concordance_df(tidydf, mbs_col, sites[0], sites[1])
SUMMS, ROCS, PREDS = parallel_setup_and_classify(
    df, sites, cls_label, SUMMS, ROCS, PREDS)

## Concatenate results dataframes and write
pd.concat(SUMMS).to_csv(args.fnsummaries, sep='\t', index=False)
pd.concat(ROCS).to_csv(args.fnrocs, sep='\t', index=False)
pd.concat(PREDS).to_csv(args.fnpreds, sep='\t', index=False)
