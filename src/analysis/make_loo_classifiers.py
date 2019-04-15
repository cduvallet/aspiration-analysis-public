#!/usr/bin/env python
"""
This script makes classifiers based on the entire communities.
"""

import argparse
import os
import sys

import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, auc, roc_curve
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

def loo_classify(df, mbs_info, random_state):
    """
    Classify aspiration vs. non-aspiration based on the given wide
    dataframe, which has subjects in rows and features in columns.

    This performs leave-one-out classification, meaning thatL for each
    patient, a model is built on all patients except that one, and then
    used to predict the status of that patient. This prediction and
    probabilities are stored to build just one ROC curve for the
    classifier.

    Parameters
    ----------
    df : pandas DataFrame
        subjects in rows, samples in features
    mbs_info : dict
        dictionary with {subject_id: 0 or 1}, all subject in
        df.index should be keys in this dictionary
    random_state : int
        random state for classifier initialization

    Returns
    -------
    resdf : pandas DataFrame
        classification result from each subject
        columns=['subject_id', 'true_label', 'predicted_label', 'prob_class_1'])
    roc_auc : float
        area under the ROC curve, calculated from the leave-one-out
        probability predictions
    fisher_p : float
        fisher's exact p-value, calculated from the confusion matrix of
        the true labels and leave-one-out predicted labels
    """

    subjects = df.index.tolist()

    res = []
    for subj in subjects:

        # Remove the subject from the training dataframe
        train_X = df.drop(subj, axis=0)

        # Code up labels as 0 or 1.
        # aspdict contains {mbs_label: 0 or 1}; mbs_info is {subject: mbs_label}
        train_Y = [mbs_info[i] for i in train_X.index.tolist()]

        # Convert to numpy matrix
        train_X = train_X.values

        # Get test data (just that subject)
        test_X = df.loc[subj].values.reshape(1, -1)
        test_Y = mbs_info[subj]

        # Set up classifier
        rf = RandomForestClassifier(n_estimators=1000, random_state=random_state)
        rf = rf.fit(train_X, train_Y)

        # Probability of being class 1
        proba = rf.predict_proba(test_X)[:, 1]

        # This just gives "if proba > 0.5, classify as 1"
        pred = rf.predict(test_X)

        res.append([subj, test_Y, pred[0], proba[0]])

    # Convert all results to dataframe
    resdf = pd.DataFrame(data=res,
        columns=['subject_id', 'true_label', 'predicted_label', 'prob_class_1'])

    # Build confusion matrix and ROC curve
    confmat = confusion_matrix(
        resdf['true_label'], resdf['predicted_label'],
        labels=[0,1])
    _, fisher_p = fisher_exact(confmat)

    fpr, tpr, thresholds = roc_curve(
        resdf['true_label'], resdf['prob_class_1'])
    roc_auc = auc(fpr, tpr)

    return (resdf, roc_auc, fisher_p)

def loo_classify_update_results(df, sites, mbs_info, random_state, PREDS, SUMMS):
    """
    Make classifier (using loo_classify()) and update lists
    of results with the classifier results and an indicator
    of which classifier was used.
    """
    ## Make classifier
    resdf, roc_auc, fisher_p = loo_classify(df, mbs_info, random_state)

    ## Calculate N asp and N normal, for storing
    n_asp = sum([mbs_info[i] == 1 for i in df.index.tolist()])
    n_nml = sum([mbs_info[i] == 0 for i in df.index.tolist()])

    ## Store results
    resdf['site'] = '-'.join(sites)
    PREDS.append(resdf)

    SUMMS.append([
        '-'.join(sites), n_asp, n_nml, df.shape[1],
        roc_auc, fisher_p])

    return PREDS, SUMMS

def write_samples(sites, out_patients, samples):
    """
    Write the list of samples to a file.

    sites : list of str
    out_patients : str, with stem of file name
    samples : list of str, sample IDs
    """
    fname = '.'.join([out_patients, '_'.join(sites), 'samples', 'txt'])
    with open(fname, 'w') as f:
        f.write('\n'.join(samples))

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    # Input files
    p.add_argument('fnotu', help='Clean OTU table. Samples in rows, OTUs '
        + 'in columns.')
    p.add_argument('fnmeta', help='Clean metadata. All samples (in rows) should'
        + ' be in the OTU table.')
    # Output files
    p.add_argument('fnsummaries', help='File where to write classifier summary '
        + 'metrics (AUC, fisher p).')
    #p.add_argument('fnrocs', help='File where to write classifier ROC data '
    #    + '(TPR, FPR).')
    p.add_argument('fnpreds', help='File where to write classifier predictions '
        + '(subjects, true labels, predictions, probability of class 1).')

    ## Optional parameters
    #p.add_argument('--rfreps', help='Number of Random Forest replicates '
    #    + '[default: %(default)s]', default=100, type=int)
    args = p.parse_args()

    # Bad form but oh well: output file for patient lists
    out_patients = 'data/patients/aspiration_classifiers'

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

    ## Prepare the data for easy manipulations
    # This creates a 3-column df with ['sample', 'otu', 'abun']
    meta = meta.dropna(subset=[mbs_col])
    tidydf = util.tidyfy_otu(df, meta, mbs_col)

    ## Make the aspiration status mapping dictionary
    # Only need to do this once, since it has every
    # subject's aspiration status
    mbs_info = tidydf\
        [['subject_id', 'mbs_consolidated']]\
        .drop_duplicates()
    # Convert to dictionary
    mbs_info = dict(zip(mbs_info['subject_id'], mbs_info['mbs_consolidated']))
    # Convert to 1 or 0
    mbs_info = {k: aspdict[mbs_info[k]] for k in mbs_info}

    ## Initialize list of results, which will get turned into dataframes
    SUMMS = []
    #ROCS = []
    PREDS = []

    ## Single-site classifiers
    sites = ['bal']
    print('-'.join(sites))
    ## Make wide dataframe with the right samples and OTUs
    df, samples = make_combined_site_df(tidydf, sites, mbs_col)
    # Track samples used in the classifier
    write_samples(sites, out_patients, samples)
    ## Make classifier and update the PREDS and SUMMS lists
    PREDS, SUMMS = loo_classify_update_results(
        df, sites, mbs_info, random_state, PREDS, SUMMS)

    sites = ['throat_swab']
    print('-'.join(sites))
    ## Make wide dataframe with the right samples and OTUs
    df, samples = make_combined_site_df(tidydf, sites, mbs_col)
    # Track samples used in the classifier
    write_samples(sites, out_patients, samples)
    ## Make classifier and update the PREDS and SUMMS lists
    PREDS, SUMMS = loo_classify_update_results(
        df, sites, mbs_info, random_state, PREDS, SUMMS)

    sites = ['gastric_fluid']
    print('-'.join(sites))
    ## Make wide dataframe with the right samples and OTUs
    df, samples = make_combined_site_df(tidydf, sites, mbs_col)
    # Track samples used in the classifier
    write_samples(sites, out_patients, samples)
    ## Make classifier and update the PREDS and SUMMS lists
    PREDS, SUMMS = loo_classify_update_results(
        df, sites, mbs_info, random_state, PREDS, SUMMS)

    ## Two-site classifiers
    sites = ['bal', 'throat_swab']
    print('-'.join(sites))
    ## Make wide dataframe with the right samples and OTUs
    df, samples = make_combined_site_df(tidydf, sites, mbs_col)
    # Track samples used in the classifier
    write_samples(sites, out_patients, samples)
    ## Make classifier and update the PREDS and SUMMS lists
    PREDS, SUMMS = loo_classify_update_results(
        df, sites, mbs_info, random_state, PREDS, SUMMS)

    sites = ['bal', 'gastric_fluid']
    print('-'.join(sites))
    ## Make wide dataframe with the right samples and OTUs
    df, samples = make_combined_site_df(tidydf, sites, mbs_col)
    # Track samples used in the classifier
    write_samples(sites, out_patients, samples)
    ## Make classifier and update the PREDS and SUMMS lists
    PREDS, SUMMS = loo_classify_update_results(
        df, sites, mbs_info, random_state, PREDS, SUMMS)

    sites = ['throat_swab', 'gastric_fluid']
    print('-'.join(sites))
    ## Make wide dataframe with the right samples and OTUs
    df, samples = make_combined_site_df(tidydf, sites, mbs_col)
    # Track samples used in the classifier
    write_samples(sites, out_patients, samples)
    ## Make classifier and update the PREDS and SUMMS lists
    PREDS, SUMMS = loo_classify_update_results(
        df, sites, mbs_info, random_state, PREDS, SUMMS)

    ## All three
    sites = ['bal', 'throat_swab', 'gastric_fluid']
    print('-'.join(sites))
    ## Make wide dataframe with the right samples and OTUs
    df, samples = make_combined_site_df(tidydf, sites, mbs_col)
    # Track samples used in the classifier
    write_samples(sites, out_patients, samples)
    ## Make classifier and update the PREDS and SUMMS lists
    PREDS, SUMMS = loo_classify_update_results(
        df, sites, mbs_info, random_state, PREDS, SUMMS)

    ## Concatenate results dataframes and write
    # SUMMS is a list of lists, turn into a dataframe
    pd.DataFrame(SUMMS,
        columns=['site', 'n_asp', 'n_nml', 'n_feats', 'auc', 'fisher_p']
        ).to_csv(args.fnsummaries, sep='\t', index=False)
    # PREDS is a list of dataframes, just concat
    pd.concat(PREDS).to_csv(args.fnpreds, sep='\t', index=False)
