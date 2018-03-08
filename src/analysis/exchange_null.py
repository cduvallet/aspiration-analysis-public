#!/usr/bin/env python
"""
This script calculate the null exchange by shuffling sample labels before
calculating correlations.

This script uses functions defined in exchange_tidy, so it needs to be in
the same directory as it.
"""
import pandas as pd
import numpy as np
import argparse
import multiprocessing

import exchange_tidy as exch

import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
import util


def get_corr_results(patients, otudf, site1, site2, site3, nthresh):
    """
    Get correlation results.

    Parameters
    ----------
    patients : list
        list of patients to consider, i.e. who have both sites sequenced
    otudf : pandas DataFrame
        dataframe with 'subject_id', 'site', and 'abun' columns
    site1, site2, site3 : str
        correlation btw sites1 and 2 is calculated, partialled on site3
    nthresh : int
        threshold number of patients which need the sites sequenced in
        order for the correlation to even be calculated
    """

    # Calculate correlation between non-zero abun in site1 and site2
    r, p = exch.corr_site12(
        otudf, patients, site1, site2)

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
    three_patients = otudf\
        .query('subject_id == @patients')\
        .query('site == @site3')\
        ['subject_id']\
        .values\
        .tolist()
    if len(three_patients) >= nthresh:
        r_partial, p_partial = exch.corr_site12_partial3(
            otudf, patients, site1, site2, site3)
    else:
        r_partial = np.nan
        p_partial = np.nan

    # Track results
    res = [site1, site2, site3,
           r, p, len(patients),
           r_partial, p_partial, len(three_patients)]

    return res

def parallel_corr((o, patients, otudf, site1, site2, site3, nthresh, n)):
    """See get_corr_results for parameter definitions."""
    print(o, site1, site2, n)
    # Get the correlation and partial correlation
    oneres = get_corr_results(
        patients, otudf, site1, site2, site3, args.nthresh)
    oneres = [o] + oneres + [n]

    return oneres


p = argparse.ArgumentParser()
p.add_argument('fn_otu', help='OTU table file path. OTUs in columns, '
                              + 'samples in rows. Relative abundances.')
p.add_argument('fn_meta', help='metadata file path.')
p.add_argument('fn_exchange', help='output exchange file path with null corrs')
p.add_argument('--nshuffle', help='number of iterations [default: %(default)s]',
    default=10, type=int)
p.add_argument('--nthresh', help='threshold for number of patients OTU must '
    + 'be present in for consideration. Set this to a number higher than the '
    + 'default when you know you"ll be excluding bugs from downstream '
    + 'analyses to speed up this code. [default: %(default)s]', default=4,
    type=int)
args = p.parse_args()

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

## Set up the inputs to parallel code
fxn_data = []
print('Setting up input data...')
for n in range(args.nshuffle):
    print(n)
    # Shuffle patient IDs within site-OTU combination
    dflong['subject_id'] = dflong\
        .groupby(['site', 'otu'])['subject_id']\
        .transform(np.random.permutation)
    # This is very slow. I could probably make it faster somehow?
    for o, otudf in dflong.groupby('otu'):
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
                site3 = [s for s in sites if s != site1 and s != site2][0]
                if len(patients) >= args.nthresh:
                    # Code is going to need a list of
                    # (o, patients, otudf, site1, site2, site3, nthresh, n)
                    fxn_data.append(
                        [o, patients, otudf, site1, site2, site3,
                         args.nthresh, n])

print('Setting up input data... Done.')

## Run the parallel code
p = multiprocessing.Pool()
res = p.map(parallel_corr, fxn_data)
p.close()
p.join()


resdf = pd.DataFrame(data=res,
                     columns=['otu', 'site1', 'site2', 'site3',
                              'r_site12', 'p_site12', 'n_site12',
                              'r_partial', 'p_partial', 'n_partial',
                              'shuffle_iter'])
resdf.to_csv(args.fn_exchange, sep='\t', index=False)
