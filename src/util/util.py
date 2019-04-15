#!/usr/bin/env python
import copy
import pandas as pd
import numpy as np
# Stats functions
from scipy.stats.mstats import kruskalwallis
from scipy.stats import ranksums, mannwhitneyu, fisher_exact
from scipy import interp
# FDR correction
from statsmodels.sandbox.stats.multicomp import multipletests
# ML functions
from sklearn.cross_validation import StratifiedKFold
from sklearn.metrics import confusion_matrix, auc, roc_curve
# Plotting
import matplotlib.pyplot as plt


def raw2abun(df):
    return df.divide(df.sum(axis=1), axis=0)

def collapse_taxonomic_contents_df(OTU_table, taxonomic_level,
                                   keep_unanno=False):
    """
    Collapses OTU table to given taxonomic level by string-matching.

    Adapated from Thomas Gurry's code at
    https://github.com/thomasgurry/amplicon_sequencing_pipeline

    Parameters
    ----------
    OTU_table : pandas dataframe
        OTUs in columns, samples in rows.
        Taxonomic levels in OTU strings should be semicolon-delimited,
        starting with kingdom level.
        Unannotated taxonomic levels should end with '__' (e.g. ''...;g__Roseburia;s__')
    taxonomic_level : str
        kingdom, phylum, class, order, family, genus, or species
    keep_unanno : bool [default: False]
        Whether to discard unannotated OTUs or keep them (e.g. g__)
        NOTE: I forgot, this isn't actually legit and shouldn't be used.
        You don't want to collapse unannotated genera bc then you have one
        thing that represents many different genera (but only has one abundance
        value), so it's apples and oranges...

    Returns
    -------
    newdf : pandas dataframe
        OTUs in columns, samples in rows.
        OTUs are collapsed to the given taxonomic level.
        Matching values (for annotated taxa) are summed for each sample.
        Values corresponding to unannotated taxa are discarded.
    """

    OTU_IDs = list(OTU_table.columns)

    # Collapse to the right level
    if(taxonomic_level == "kingdom"):
        OTU_taxa = [OTU_ID.split(';')[0] for OTU_ID in OTU_IDs]
    if(taxonomic_level == "phylum"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:2]) for OTU_ID in OTU_IDs]
    if(taxonomic_level == "class"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:3]) for OTU_ID in OTU_IDs]
    if(taxonomic_level == "order"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:4]) for OTU_ID in OTU_IDs]
    if(taxonomic_level == "family"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:5]) for OTU_ID in OTU_IDs]
    if(taxonomic_level == "genus"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:6]) for OTU_ID in OTU_IDs]
    if(taxonomic_level == "species"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:7]) for OTU_ID in OTU_IDs]

    # Get indices in original dataframe corresponding to each unique taxon
    taxa_indices = {}
    for i in range(len(OTU_taxa)):
        # Append index if the taxa is annotated
        if not OTU_taxa[i].endswith("__"):
            # Initialize list if not already in dict
            if OTU_taxa[i] not in taxa_indices:
                taxa_indices[OTU_taxa[i]] = []
            taxa_indices[OTU_taxa[i]].append(i)
        # If we should keep unannotated taxa, also update that
        if keep_unanno and OTU_taxa[i].endswith("__"):
            # Initialize list if not already in dict
            if OTU_taxa[i] not in taxa_indices:
                taxa_indices[OTU_taxa[i]] = []
            taxa_indices[OTU_taxa[i]].append(i)

    # Make new empty df with the same samples as original df
    # and taxa names
    newdf = pd.DataFrame(index=OTU_table.index,
                         columns=taxa_indices.keys(),
                         data=0)

    # Get sample contents for each taxa of the chosen level and put into newdf
    for key in taxa_indices:
        indices = taxa_indices[key]
        newcol = OTU_table.iloc[:, indices].sum(axis=1)
        newdf[key] = copy.copy(newcol)

    return newdf

def compare_otus_teststat(df, Xsmpls, Ysmpls, method='kruskal-wallis', multi_comp=None):
    """
    Compares columns between Xsmpls and Ysmpls, with statistical method=method.
    Returns dataframe with both the qvals ('p') and test statistic ('test-stat')

    parameters
    ----------
    df             dataframe, samples are in rows and OTUs in columns
    X,Ysmpls       list of samples to compare
    method         statistical method to use for comparison
    multi_comp     str, type of multiple comparison test to do.
                   Currently accepts 'fdr' or None

    outputs
    -------
    results        dataframe with OTUs in rows and 'p' and 'test-stat' in columns

    """
    if method == 'kruskal-wallis':
        pfun = kruskalwallis
    elif method == 'wilcoxon' or method == 'ranksums':
        pfun = ranksums
    elif method == 'mann-whitney':
        pfun = mannwhitneyu
        # Note: prob wanna add some kwargs here to say whether 2sided or not

    results = pd.DataFrame(index=df.columns, columns=['test-stat', 'p'])
    for o in df.columns:
        try:
            h, p = pfun(df.loc[Xsmpls, o].tolist(), df.loc[Ysmpls, o].tolist())
        except:
            p = 1
            h = 0
        results.loc[o, 'p'] = p
        results.loc[o, 'test-stat'] = h

    if multi_comp == 'fdr':
        _, results['q'], _, _ = multipletests(results['p'], method='fdr_bh')

    return results

def get_sites():
    return ['bal', 'gastric_fluid', 'throat_swab']

def get_metacols():
    """
    Returns some metadata columns of interest. I think this should be
    OK to use as a stock function, but we'll see...
    """
    metacols =  ['subject_id', 'site',
                 'mbs_consolidated', 'ppi_consolidated',
                 'number of full colum events/total events',
                 'Total number of reflux episodes (acid+non-acid)',
                 'percent proximal total',
                 'percent distal total',
                 'Was Bile CA detected?', 'Was Bile DCA detected?',
                 'Was Bile LCA detected?', 'Was Bile TCA detected?']
    return metacols

def add_metadata_to_beta_list(beta_lst, meta, metacols):
    """
    Add metadata values to the samples in beta_lst.

    Parameters
    ----------
    beta_lst : list
        list of lists, where each sublist has:
        [sample1, sample2, beta]
    meta : pandas DataFrame
        dataframe with metadata columns and 'subject_id', 'site'
        columns
    metacols : list
        additional columns to add for within-patient comparisons

    Returns
    -------
    df_beta : pandas DataFrame
        dataframe with the original list values ('sample1', 'sample2',
        'beta'), 'site1', 'site2', 'site_comparison', 'patient_comp',
        'subject', and metacols columns
    """

    # And get the other metadata values...
    smpl2subj = dict(zip(meta.index, meta['subject_id']))
    smpl2site = dict(zip(meta.index, meta['site']))
    smpl2batch = dict(zip(meta.index, meta['batch']))

    beta_all = []
    for lst in beta_lst:
        s1 = lst[0]
        s2 = lst[1]
        subj1 = smpl2subj[s1]
        subj2 = smpl2subj[s2]

        site1 = smpl2site[s1]
        site2 = smpl2site[s2]
        site_comp = site1 + '-' + site2

        batch1 = smpl2batch[s1]
        batch2 = smpl2batch[s2]

        if subj1 == subj2:
            patient_comp = 'within'
            subject = subj1
            metadata = list(meta.loc[s1, metacols].values)
        else:
            patient_comp = 'between'
            subject = np.nan
            metadata = len(metacols)*[np.nan]

        beta_all.append(
            lst
            + [site1, site2, site_comp, patient_comp, subject, batch1, batch2]
            + metadata)

    df_beta = pd.DataFrame(
        data=beta_all,
        columns=['sample1', 'sample2', 'beta', 'site1', 'site2',
                 'site_comparison', 'patient_comp', 'subject',
                 'batch1', 'batch2'] +
                 metacols)

    return df_beta

## Classification
def cv_and_roc(rf, X, Y, num_cv=5, random_state=None, shuffle=False):
    """
    Perform cross validated training and testing and return the aggregate
    interpolated ROC curve and confusion matrices.

    Parameters
    ----------
    rf : any sklearn classifier object
    X : array-like or sparse matrix, shape = [n_samples, n_features]
        The input samples to be split into train and test folds and
        cross-validated.
    Y : list or array
        array-like, shape = [n_samples] or [n_samples, n_outputs]
        The target values (class labels in classification).
    num_cv : int (default: 5)
        number of cross-validation folds
    random_state : int (default 12345)
        random state seed for StratifiedKFold

    Returns
    -------
    d : dict with the following keys:
        'roc_auc': area under the ROC curve
        'conf_mat': confusion matrix array, looks like:
                              pred
                             0   1
                    true  0  -   -
                          1  -   -
        'fisher_p': fisher exact pvalue of conf_mat above
        'y_probs': probability of being class 1
        'y_trues': true labels
        'mean_fpr', 'mean_tpr': interpolated values used to build ROC curve

    """
    if isinstance(Y, list):
        Y = np.asarray(Y)
    cv = StratifiedKFold(Y, num_cv, shuffle=True, random_state=random_state)
    tpr_list = []
    auc_list = []
    mean_fpr = np.linspace(0, 1, 100)
    conf_mat = np.asarray([[0,0],[0,0]])
    y_probs = np.empty_like(Y, dtype=float)
    y_predictions = np.empty_like(Y)

    for train_index, test_index in cv:
        X_train, X_test = X[train_index], X[test_index]
        Y_train, Y_test = Y[train_index], Y[test_index]
        probs = rf.fit(X_train, Y_train).predict_proba(X_test)[:,1]
        preds = rf.predict(X_test)

        # Store probability and predicted values for X_test
        y_probs[test_index] = probs
        y_predictions[test_index] = preds

        # Compute ROC curve and area under the curve
        fpr, tpr, thresholds = roc_curve(Y_test, probs)
        roc_auc = auc(fpr, tpr)
        auc_list.append(roc_auc)
        tpr_list.append(interp(mean_fpr, fpr, tpr))
        tpr_list[-1][0] = 0.0

        # Compute confusion matrix
        conf_mat += confusion_matrix(Y_test, preds, labels=[0,1])

    # Get fisher p value, sensitivity/specificity, positive/negative pred value
    _, fisher_p = fisher_exact(conf_mat)
    # conf_mat: rows are the true labels, columns are the predicted labels
    conf_mat = pd.DataFrame(conf_mat)
    conf_mat.columns = ['pred 0', 'pred 1']
    conf_mat.index = ['true 0', 'true 1']
    ppv = conf_mat.loc['true 1', 'pred 1'] / float(conf_mat['pred 1'].sum())
    npv = conf_mat.loc['true 0', 'pred 0'] / float(conf_mat['pred 0'].sum())
    spec = \
        conf_mat.loc['true 0', 'pred 0'] / float(conf_mat.loc['true 0'].sum())
    sens = \
        conf_mat.loc['true 1', 'pred 1'] / float(conf_mat.loc['true 1'].sum())

    return {i: j for i, j in
            zip(('auc_list', 'conf_mat', 'mean_fpr', 'tpr_list',
                'fisher_p', 'y_prob', 'y_predictions', 'y_true',
                'sensitivity', 'specificity', 'ppv', 'npv'),
               (auc_list, conf_mat, mean_fpr, tpr_list, fisher_p,
               y_probs, y_predictions, Y,
               sens, spec, ppv, npv))}

def plot_auc_from_list(tpr_list, mean_fpr, color='blue', ax=None):
    """
    tpr_list : list
        list of true positive rates (len(tpr_list) = number of CV folds)
    mean_fpr : list
        np.linspace(0, 1) with the same number of elements as tpr_list
    color : str
        color of the mean ROC line to plot

    Adapted from Isaac's code at
    https://github.com/irockafe/revo_healthcare/blob/master/src/visualization/viz.py
    """
    # get mean tpr and fpr
    mean_tpr = np.mean(tpr_list, axis=0)
    # make sure it ends up at 1.0
    mean_tpr[-1] = 1.0

    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(mean_fpr, mean_tpr, color=color)
    # Plot luck line
    ax.plot([0, 1], [0, 1], linestyle='--', color='k', alpha=0.5)

    # plot 1-std
    std_tpr = np.std(tpr_list, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color=color,
                    alpha=0.2)

    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')

    return ax

def tidyfy_otu(df, meta, mbscol, cols=None):
    """
    Turn OTU table into a tidy dataframe with subject, sample ID, site, OTU,
    OTU labeled with the site, abundance, and aspiration status.

    This script also removes some unwanted samples (e.g. second time point,
    lung transplant samples)

    To make our OTU-based classifiers (after tidyfying), we will:
    1. Query tidyotu for the site(s) of interest and drop any samples
       without aspiration metadata.
    2. Pivot queried table to have subjects in rows and otu-w-site in columns
    3. Assemble the output (Y) from the metadata, using a subject to
       aspiration status dict

    Parameters
    ----------
    df : pandas DataFrame
        OTUs in columns and samples in rows
    meta : pandas DataFrame
        metadata in columns and samples in rows
    mbscol : str
        column in meta with aspiration status
    cols : str
        additional metadata columns to keep. Should be a list, even if it's
        only one column.
    """

    print('Tidying data... ')
    tidyotu = pd.melt(df.reset_index(), id_vars='index', value_name='abun',
                      var_name='otu')
    tidyotu = tidyotu.rename(columns={'index': 'sample'})

    # Remove any of the 2nd time point samples
    exclude = ['2', 'F', 'sick', 'F2T']
    for s in exclude:
        tidyotu = tidyotu[~tidyotu['sample'].str.endswith(s)]
    # And remove any lung transplant samples
    tidyotu = tidyotu[~tidyotu['sample'].str.startswith('05')]

    # Add the subject ID, site, and aspiration metadata associated with each sample
    keepcols = [mbscol, 'site', 'subject_id']
    if cols is not None:
        keepcols += cols
    tidyotu = tidyotu.merge(meta[keepcols], left_on='sample', right_index=True)
    # Label OTUs with their associated site
    tidyotu['otu_w_site'] = tidyotu.apply(
        lambda row: row['otu'] + '--' + row['site'],
        axis=1)
    print('Finished.')
    return tidyotu

def convert_to_latex(row):
    """
    Convert a row from a pandas DataFrame (i.e. a Series)
    into different cells in a Latex table. Also replaces
    underscores with escaped underscore.
    """
    return ' & '.join(['\_'.join(i.split('_')) for i in row.astype(str)]) + ' \\\ '
