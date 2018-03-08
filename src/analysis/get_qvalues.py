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

def univariate_by_site(aspmeta, df):
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

    Returns
    -------
    allqdf : pandas dataframe
        tidy dataframe with ['otu', 'site', 'p', 'q', 'n_asp', 'n_nml'] columns,
        where 'n_asp' and 'n_nml' are the number of samples in each
        category, respectively.
    """
    allq = []
    for site, sitedf in aspmeta.groupby('site'):
        asp_smpls = sitedf.query('mbs_consolidated == "Aspiration/Penetration"').index
        nml_smpls = sitedf.query('mbs_consolidated == "Normal"').index
        results = util.compare_otus_teststat(df, asp_smpls, nml_smpls,
                                             method='wilcoxon', multi_comp='fdr')
        results.index.name = 'otu'
        results = results.reset_index()
        results['site'] = site
        results['n_asp'] = len(asp_smpls)
        results['n_nml'] = len(nml_smpls)
        allq.append(results)

    return pd.concat(allq, ignore_index=True)

p = argparse.ArgumentParser()
p.add_argument('fnotu', help='clean OTU table. samples in rows, OTUs in columns')
p.add_argument('fnmeta', help='clean metadata. All samples (in rows) should'
    + ' be in the OTU table')
p.add_argument('fout', help='file to write univariate p and q-values to')
args = p.parse_args()

df = pd.read_csv(args.fnotu, sep='\t', index_col=0)
meta = pd.read_csv(args.fnmeta, sep='\t', index_col=0)

## For each site, compare aspirators and non-aspirators
aspmeta = meta.dropna(subset=['mbs_consolidated'])
qvals_otu = univariate_by_site(aspmeta, df)
qvals_otu['level'] = 'otu'

## What if I do univariate at genus level?
genusdf = util.collapse_taxonomic_contents_df(df, 'genus').loc[aspmeta.index]
## Note: there are like 10 samples with fewer than 50% of reads annotated to genus level
# genusdf.loc[aspmeta.index].sum(axis=1).sort_values()
## 10 of these seem dominated by k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__;g__;s__;d__denovo34,
## which BLASTS to Bacillus thermoamylovorans
qvals_genus = univariate_by_site(aspmeta, genusdf)
qvals_genus['level'] = 'genus'

allqdf = pd.concat((qvals_otu, qvals_genus), ignore_index=True)
allqdf.to_csv(args.fout, sep='\t', index=False)
