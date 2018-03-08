#!/usr/bin/env python
"""
This script calculates the prevalence of bacteria in two sites.

Given:
- tidy exchange file
- r and q thresholds (otherwise none)
- otu and metadata

Return:
- prevalence of shared bugs (those that meet criteria), for each metadata category
- columns = otu, meta_col, meta_val, prevalence
"""
import argparse

import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests

import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
import util

def site_comp_metadata(meta, sites, metacols):
    """
    Make tidy metadata dataframe containing samples from the same
    subject, for subjects which have two sites sampled.

    Parameters
    ----------
    meta : pandas dataframe
        metadata dataframe with column 'site'
    sites : list
        list of sites in meta['site'] to consider
    metacols : list
        columns with additional metadata to also include

    Returns
    -------
    tidymeta : pandas dataframe
        Tidy data with columns ['sample1', 'sample2', 'site_comparison',
        metacols]. Contains data for samples which come from the same subject
        and are from different sites. Excludes samples which end in
        'F', '2', 'sick' or start with '05'.
    """
    lst_tidymeta = []

    for site1 in sites:
        for site2 in sites[sites.index(site1)+1:]:
            # Drop any patients who don't have OTU present in both sites
            # Note: I need to groupby subject_id and site first because
            # some patients have multiple samples per site, so I need to
            # collapse those...
            patients = meta\
                        .query('(site == @site1) | (site == @site2)')\
                        .groupby(['subject_id', 'site'])\
                        .size()\
                        .reset_index()\
                        .groupby('subject_id')\
                        .size()
            patients = patients[patients == 2].index.tolist()

            for p in patients:
                s1 = meta.query('(subject_id == @p) & (site == @site1)')\
                    .index.values
                s2 = meta.query('(subject_id == @p) & (site == @site2)')\
                    .index.values

                #if len(s1) > 1 or len(s2) > 1:
                #    s1 = [i for i in s1
                #          if not i.endswith('2')
                #          and not i.endswith('F')
                #          and not i.endswith('sick')
                #          and not i.endswith('F2T')
                #          and not i.endswith('F2')
                #          and not i.startswith('05')]
                #    s2 = [i for i in s2
                #          if not i.endswith('2')
                #          and not i.endswith('F')
                #          and not i.endswith('sick')
                #          and not i.endswith('F2T')
                #          and not i.endswith('F2')
                #          and not i.startswith('05')]
                if len(s1) != 1 or len(s2) != 1:
                    continue

                # s1 and s2 are index objects...
                s1 = s1[0]
                s2 = s2[0]
                lst_tidymeta.append(
                    [p, s1, s2, site1 + '-' + site2] +
                    [meta.loc[s1, i] for i in metacols])

    tidymeta = pd.DataFrame(data=lst_tidymeta,
                            columns=['subject', 'sample1', 'sample2',
                                     'site_comparison'] + metacols)
    return tidymeta

def calculate_exchange_preva(col, s1smpls, s2smpls):
    """
    Calculate the exchange/sharedness value for abundances in col
    between s1smpls and s2smpls.

    Parameters
    ----------
    col : pandas Series
        Values to correlate and calculate 'exchange' for, e.g. relative
        abundances of one OTU. Index should have at least s1smpls and
        s2smpls. Values should be 0 if the OTU is not present and greater
        than zero if it is.

    s1smpls, s2smpls :  lists
        Samples to consider. Should be the same length and paired (i.e. the
        first sample in s1smpls comes from the same patient as the first
        sample in s2smpls.)

    Returns
    -------
    sharedness : float
        Percent of (s1, s2) pairs where both samples are non-zero.
        In other words, percentage of x, y points which are off the axes.
    """
    preva = sum([i and j for i, j in zip(
                    (col.loc[s1smpls] > 0).values,
                    (col.loc[s2smpls] > 0).values)])\
            /float(len(s1smpls))
    return preva

p = argparse.ArgumentParser()
p.add_argument('exchange', help='file with exchange (correlations)')
p.add_argument('otu', help='otu table, OTUs in columns samples in rows. '
    + 'Can be either relative abundance or reads - only thing that gets '
    + 'counted is presence/absence.')
p.add_argument('meta', help='metadata file')
p.add_argument('fout', help='out file')
p.add_argument('--nthresh', help='number of patients that need to have the '
    + 'otu present in both sites for it to be considered for being exchanged.',
    default=5, type=int)
p.add_argument('--qthresh', help='q-value threshold to use to consider '
    + 'which OTUs get prevalence calculated for.', default=1.0, type=float)

# Parse args
args = p.parse_args()
exchange = pd.read_csv(args.exchange, sep='\t')
otu = pd.read_csv(args.otu, sep='\t', index_col=0)
meta = pd.read_csv(args.meta, sep='\t', index_col=0)
nthresh = args.nthresh
qthresh = args.qthresh

# Get subject-wise metadata and samples for all site comparisons
sites = util.get_sites()
metacols =  ['subject_id', 'site',
             'mbs_consolidated', 'ppi_consolidated',
             'Was Bile CA detected?', 'Was Bile DCA detected?',
             'Was Bile LCA detected?', 'Was Bile TCA detected?']
# Remove unwanted samples from the metadata
meta = meta[~meta['sample_id.1'].str.endswith('F')]
meta = meta[~meta['sample_id.1'].str.endswith('2')]
meta = meta[~meta['sample_id.1'].str.endswith('sick')]
meta = meta[~meta['sample_id.1'].str.endswith('F2T')]
meta = meta[~meta['sample_id.1'].str.endswith('F2')]
meta = meta[~meta['sample_id.1'].str.startswith('05')]

# Tidyfy metadata
tidymeta = site_comp_metadata(meta, sites, metacols)
metacols = metacols[2:]

# Query OTUs that match r and qthresh criteria, using partial corr results
exchange = exchange.query('n_partial >= @nthresh')
_, exchange['q_partial'], _, _ = multipletests(
    exchange['p_partial'], method='fdr_bh')
exchange = exchange.query('(q_partial < @qthresh)')
exchange['site_comparison'] = exchange['site1'] + '-' + exchange['site2']

allres = []
for g, subexchange in exchange.groupby('site_comparison'):
    print(g)
    sitemeta = tidymeta.query('site_comparison == @g')

    # These are the samples with data for both sites
    s1smpls = sitemeta['sample1']
    s2smpls = sitemeta['sample2']

    # These are the exchanged OTUs
    otus = subexchange['otu']

    # This returns an unnamed Series with index = otus
    res = otu[otus]\
        .apply(calculate_exchange_preva, args=(s1smpls, s2smpls)).T\
        .reset_index()

    res.columns = ['otu', 'prevalence_exchange']
    res['meta_var'] = 'all_patients'
    res['meta_val'] = 'all_patients'
    res['site_comparison'] = g
    res['n_patients'] = len(s1smpls)
    allres.append(res)

    ## Now calculate exchange for each subgroup of patients
    for mcol in metacols:
        print(mcol)
        for mval, sitemeta_sub in sitemeta.dropna(subset=[mcol]).groupby(mcol):
            # These are the samples with data for both sites
            s1smpls = sitemeta_sub['sample1']
            s2smpls = sitemeta_sub['sample2']

            # This returns an unnamed Series with index = otus
            res = otu[otus]\
                .apply(calculate_exchange_preva, args=(s1smpls, s2smpls)).T\
                .reset_index()
            res.columns = ['otu', 'prevalence_exchange']
            res['meta_var'] = mcol
            res['meta_val'] = mval
            res['site_comparison'] = g
            res['n_patients'] = len(s1smpls)
            allres.append(res)

# allres is a list of dataframes with the same column names
resdf = pd.concat(allres)
resdf.to_csv(args.fout, sep='\t', index=False)
