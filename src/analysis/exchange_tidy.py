#!/usr/bin/env python
"""
This script calculates the correlation and partial correlations of all
OTUs across sites.

For each OTU, this script calculates the correlation between the abundance
in two sites across patients (e.g. corr(siteA_i, siteB_i) for i patients)

It writes a file with the following columns:
otu, site1, site2, site_control, r, p, n_patients

Where site1 and site2 are the two sites being correlated, site_control
is the site controlled for (in the case of partial correlation) or NaN,
and n_patients is the number of patients with non-zero abundance in
both site1 and site2.
"""
import argparse
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
import util

def corr_formula(r_xy, r_xz, r_zy):
    return (r_xy - r_xz*r_zy) / np.sqrt((1 - r_xz**2)*(1 - r_zy**2))

def my_partial_corr(x, y, z):
    """
    Calculate partial correlation of x and y conditioned on z, using
    Spearman correlations and the formula from Wikipedia:

    r_xy_z = (r_xy - r_xz * r_zy) / sqrt((1 - r_xz^2)(1 - r_zy^2))
    """
    r_xy, _ = spearmanr(x, y)
    r_xz, _ = spearmanr(x, z)
    r_zy, _ = spearmanr(z, y)
    return corr_formula(r_xy, r_xz, r_zy)

def partial_corr_p(x, y, z, r_true):
    """
    Calculate the pvalue of the partial correlation of x and y conditioned
    on z.

    Shuffle values of x and y [TODO: do I need to shuffle z also?] to get
    the null.

    Return the one-tailed p-value given by r_null > r_true (probability of
    finding a larger correlation than what was observed).
    """
    r_nulls = []
    iters = 2000
    for i in range(iters):
        nullx = np.random.permutation(x)
        nully = np.random.permutation(y)
        r_nulls.append(my_partial_corr(nullx, nully, z))
    return sum(r_nulls > r_true)/float(iters)

def corr_site12(otudf, patients, site1, site2, shuffle=False):
    """
    Calculate spearman correlation of abundances in site1 and site2.
    If shuffle = True, shuffles the abundances within sites (i.e.
    scrambles the patient labels)
    """

    twositesdf = otudf.dropna().query('subject_id == @patients')
    # Make sure that samples from the same patient are grouped
    # I think the .query call does this automatically, but just in case...
    twositesdf = twositesdf.sort_values(by='subject_id')

    # Calculate correlation between non-zero abun in site1 and site2
    if shuffle:
        r, p = spearmanr(
                twositesdf.query('site == @site1')['abun'].sample(frac=1),
                twositesdf.query('site == @site2')['abun'].sample(frac=1))
    else:
        r, p = spearmanr(
                twositesdf.query('site == @site1')['abun'],
                twositesdf.query('site == @site2')['abun'])
    return r, p

def corr_site12_partial3(otudf, patients, site1, site2, site3, shuffle=False):
    """
    Get partial correlation of sites 1 and 2 in patients, partialled on
    site3.

    If shuffle=True, shuffles the abundances within sites (i.e. scrambles
    patient labels)
    """

    # Convert to wide dataframe with sites in columns,
    # subjects in rows
    threesitesdf = otudf\
        .query('subject_id == @patients')\
        .pivot(index='subject_id',
               columns='site', values='abun')

    # Make sure that samples from the same patient are grouped
    # I think the .query call does this automatically, but just in case...
    threesitesdf = threesitesdf.reset_index().sort_values(by='subject_id')

    # If all values in site3 are nan, replace with some noise
    # Otherwise fill NaNs in site3 with the minimum (i.e. the zero value)
    minval = threesitesdf.min().min()
    if threesitesdf[site3].isnull().all():
        threesitesdf[site3] = \
            np.log10(10**minval +
                     np.random.normal(0, 1e-6, threesitesdf[site3].shape))
    else:
        threesitesdf[site3] = threesitesdf[site3].fillna(minval)

    if shuffle:
        # pd.apply(axis=0) applies function to columns independently
        # np.random.permutation permutes values, not in-place
        threesitesdf = threesitesdf.apply(np.random.permutation, axis=0)

    # Partial correlation btw sites 1 and 2, controlled on site3
    r_partial = my_partial_corr(
                        threesitesdf[site1],
                        threesitesdf[site2],
                        threesitesdf[site3])

    # Permutation-based partial correlation p value
    p_partial = partial_corr_p(
                        threesitesdf[site1],
                        threesitesdf[site2],
                        threesitesdf[site3],
                        r_partial)

    return r_partial, p_partial

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('fn_otu', help='OTU table file path. OTUs in columns, '
                                  + 'samples in rows. Relative abundances.')
    p.add_argument('fn_meta', help='metadata file path.')
    p.add_argument('fn_exchange', help='output exchange file path')
    p.add_argument('--shuffle', action='store_true')
    args = p.parse_args()

    # Set random state, so that calls to permutations, etc
    # return the same permutations across different runs of
    # this script
    np.random.seed(12345)

    ## Read in OTU table and metadata
    df = pd.read_csv(args.fn_otu, sep='\t', index_col=0)
    meta = pd.read_csv(args.fn_meta, sep='\t', index_col=0)

    # Convert to log and turn 0's into NaNs
    df = np.log10(df).replace(-np.inf, np.nan)
    df.index.name = 'sample'

    ## Melt OTU table and add 'site' and 'subject_id' columns
    dflong = pd.melt(df.reset_index(), id_vars=['sample'], var_name='otu',
                     value_name='abun')
    dflong = dflong.merge(
                meta[['site', 'subject_id']],
                left_on='sample', right_index=True)
    # Remove any samples corresponding to second time point
    # Fundo samples end in F2 (corresponding to patient X-F1)
    # gastric/throat time points end in GI/GF or TI/TF
    exclude = ['2', 'F', 'sick', 'F2T']
    for s in exclude:
        dflong = dflong[~dflong['sample'].str.endswith(s)]
    # And remove any lung transplant samples
    dflong = dflong[~dflong['sample'].str.startswith('05')]

    ## Correlations
    sites = util.get_sites()
    res = []
    for o, otudf in dflong.groupby('otu'):
        print(o)
        # For each pairwise site combo:
        for site1 in sites:
            for site2 in sites[sites.index(site1)+1:]:
                # Drop any patients who don't have OTU present in both sites
                # Note: I need to groupby subject_id and site first because
                # some patients have multiple samples per site, so I need to
                # collapse those...
                patients = otudf\
                            .dropna()\
                            .query('(site == @site1) | (site == @site2)')\
                            .groupby(['subject_id', 'site'])\
                            .size()\
                            .reset_index()\
                            .groupby('subject_id')\
                            .size()
                patients = patients[patients == 2].index.tolist()

                if len(patients) > 4:
                    # Calculate correlation between non-zero abun in site1 and site2
                    r, p = corr_site12(
                        otudf, patients, site1, site2, shuffle=args.shuffle)

                    ## Partial correlation conditioned on site 3
                    ## Algorithm: of the patients with non-zero abundance in both
                    ## sites1 and 2, keep only patients who had site3 sequenced
                    ## (regardless of whether OTU was non-zero in that site).
                    ## If OTU was non-zero in site3, fill that NaN value with the
                    ## min value of both sites1 and 2. Then, find the partial
                    ## correlation between sites1 and 2 conditioned on site3, using
                    ## the formula from Wikipedia with spearman correlations.

                    # Of the pairwise patients, drop any who don't have
                    # site3 sequenced
                    site3 = [s for s in sites if s != site1 and s != site2][0]
                    three_patients = otudf\
                        .query('subject_id == @patients')\
                        .query('site == @site3')\
                        ['subject_id']\
                        .values\
                        .tolist()
                    if len(three_patients) > 4:
                        r_partial, p_partial = corr_site12_partial3(
                            otudf, patients, site1, site2, site3,
                            shuffle=args.shuffle)
                    else:
                        r_partial = np.nan
                        p_partial = np.nan

                    # Track results
                    res.append([o, site1, site2, site3,
                                r, p, len(patients),
                                r_partial, p_partial, len(three_patients),
                                args.shuffle])

    resdf = pd.DataFrame(data=res,
                         columns=['otu', 'site1', 'site2', 'site3',
                                  'r_site12', 'p_site12', 'n_site12',
                                  'r_partial', 'p_partial', 'n_partial',
                                  'shuffled'])
    resdf.to_csv(args.fn_exchange, sep='\t', index=False)
