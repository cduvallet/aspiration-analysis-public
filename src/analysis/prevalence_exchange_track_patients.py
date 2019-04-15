#!/usr/bin/env python
"""
This script is essentialy a bare-bones copy of the script that
calculates the prevalence of exchanged OTUs in aspirators and
non-aspirators.

Its goal is to track the samples I used in calculating these
prevalences. I'll only track the aspiration patients (because
these are the only ones I report in the paper), but you could
easily extend it to the other metadata columns.
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

p = argparse.ArgumentParser()
p.add_argument('exchange', help='file with exchange (correlations)')
p.add_argument('otu', help='otu table, OTUs in columns samples in rows. '
    + 'Can be either relative abundance or reads - only thing that gets '
    + 'counted is presence/absence.')
p.add_argument('meta', help='metadata file')
p.add_argument('fout', help='stem of out file, with the samples used in '
    + 'each calculation')
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

# Unlike original code, we're only going to look at the aspiration patients
#metacols = metacols[2:]
metacols = ['mbs_consolidated']

# I really don't need this code to track the patients, but whatever
# Query OTUs that match r and qthresh criteria, using partial corr results
exchange = exchange.query('n_partial >= @nthresh')
_, exchange['q_partial'], _, _ = multipletests(
    exchange['p_partial'], method='fdr_bh')
exchange = exchange.query('(q_partial < @qthresh)')
exchange['site_comparison'] = exchange['site1'] + '-' + exchange['site2']

for g, subexchange in exchange.groupby('site_comparison'):
    all_samples = []
    print(g)
    # The code below looks at all patients with the two
    # sites. Ignore it for these purposes.
    #
    sitemeta = tidymeta.query('site_comparison == @g')
    #
    # # These are the samples with data for both sites
    # s1smpls = sitemeta['sample1']
    # s2smpls = sitemeta['sample2']
    #
    # # These are the exchanged OTUs
    # otus = subexchange['otu']
    #
    # # This returns an unnamed Series with index = otus
    # res = otu[otus]\
    #     .apply(calculate_exchange_preva, args=(s1smpls, s2smpls)).T\
    #     .reset_index()
    #
    # res.columns = ['otu', 'prevalence_exchange']
    # res['meta_var'] = 'all_patients'
    # res['meta_val'] = 'all_patients'
    # res['site_comparison'] = g
    # res['n_patients'] = len(s1smpls)
    # allres.append(res)

    ## Now calculate exchange for each subgroup of patients
    for mcol in metacols:
        print(mcol)
        for mval, sitemeta_sub in sitemeta.dropna(subset=[mcol]).groupby(mcol):
            # These are the samples with data for both sites
            s1smpls = sitemeta_sub['sample1'].tolist()
            s2smpls = sitemeta_sub['sample2'].tolist()

            all_samples += s1smpls
            all_samples += s2smpls

            # Similar code as above in the original, which I removed
            # here to reduce clutter

        # Write the samples used in that site comparison and metadata
        # column to a file
        fname = '.'.join([args.fout, g, mcol, 'samples', 'txt'])
        with open(fname, 'w') as f:
            f.write('\n'.join(all_samples))
