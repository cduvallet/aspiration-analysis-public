#!/usr/bin/env python
"""
Reads in metadata and OTU table and performs univariate test comparing
aspirators/penetrators with non-aspirators.
"""

import argparse
import pandas as pd
import os, sys

# add util/ to path
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
import util

def univariate_by_site(aspmeta, df, method, meta_col, case_lbl, ctrl_lbl):
    """
    Iterate through all sites in aspmeta['site'] and do univariate
    test on 'Aspiration/Penetration' vs. 'Normal' samples.

    Parameters
    ----------
    aspmeta : pandas dataframe
        samples in rows, 'site' and 'mbs_consolidated' columns (at least)
    df : pandas dataframe
        samples in rows (as in aspmeta), OTUs/genera in columns.
        Stats comparison is performed directly on values (i.e. no collapsing
        or relative abundance-ing is done)
    method : str
        method to use to call qvalues. see util.compare_otus_teststat for
        available options
    meta_col, case_lbl, ctrl_lbl : str
        metadata column and associate labels to compare between

    Returns
    -------
    allqdf : pandas dataframe
        tidy dataframe with ['otu', 'site', 'p', 'q', 'n_asp', 'n_nml'] columns,
        where 'n_asp' and 'n_nml' are the number of samples in each
        category, respectively.
    """
    allq = []
    for site, sitedf in aspmeta.groupby('site'):
        asp_smpls = sitedf[
            sitedf[meta_col] == case_lbl
            ].index.tolist()
        nml_smpls = sitedf[
            sitedf[meta_col] == ctrl_lbl
            ].index.tolist()

        # Remove OTUs which are absent in these samples
        subdf = df.loc[asp_smpls + nml_smpls]
        subdf = subdf.loc[:, subdf.sum() != 0]

        results = util.compare_otus_teststat(subdf, asp_smpls, nml_smpls,
                                             method=method, multi_comp='fdr')
        results.index.name = 'otu'
        results = results.reset_index()
        results['site'] = site
        results['n_asp'] = len(asp_smpls)
        results['n_nml'] = len(nml_smpls)
        results['total_taxa'] = subdf.shape[1]
        allq.append(results)

    return pd.concat(allq, ignore_index=True)

p = argparse.ArgumentParser()
p.add_argument('fnotu', help='clean OTU table. samples in rows, OTUs in columns')
p.add_argument('fnmeta', help='clean metadata. All samples (in rows) should'
    + ' be in the OTU table')
p.add_argument('fout', help='file to write univariate p and q-values to')
p.add_argument('--method', help='method to call pvalues', default='wilcoxon')
p.add_argument('--metacol', help='metadata column', default='mbs_consolidated')
p.add_argument('--caselabel', help='label for cases', default='Aspiration/Penetration', type=str)
p.add_argument('--ctrllabel', help='label for controls', default='Normal', type=str)
args = p.parse_args()

df = pd.read_csv(args.fnotu, sep='\t', index_col=0)
meta = pd.read_csv(args.fnmeta, sep='\t', index_col=0)

meta_col = args.metacol
case_lbl = args.caselabel
ctrl_lbl = args.ctrllabel

## For each site, compare aspirators and non-aspirators
aspmeta = meta.dropna(subset=[meta_col])
# Convert the non-aspiration columns into strings
if aspmeta[meta_col].dtype == 'int' or aspmeta[meta_col].dtype == 'float64':
    aspmeta[meta_col] = aspmeta[meta_col].astype(int).astype(str)
qvals_otu = univariate_by_site(aspmeta, df, args.method, meta_col, case_lbl, ctrl_lbl)
qvals_otu['level'] = 'otu'

## What if I do univariate at genus level?
genusdf = util.collapse_taxonomic_contents_df(df, 'genus').loc[aspmeta.index]
## Note: there are like 10 samples with fewer than 50% of reads annotated to genus level
# genusdf.loc[aspmeta.index].sum(axis=1).sort_values()
## 10 of these seem dominated by k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__;g__;s__;d__denovo34,
## which BLASTS to Bacillus thermoamylovorans
qvals_genus = univariate_by_site(aspmeta, genusdf, args.method, meta_col, case_lbl, ctrl_lbl)
qvals_genus['level'] = 'genus'

allqdf = pd.concat((qvals_otu, qvals_genus), ignore_index=True)
#allqdf = qvals_otu
allqdf.to_csv(args.fout, sep='\t', index=False)
