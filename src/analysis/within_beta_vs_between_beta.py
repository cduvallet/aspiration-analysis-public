#!/usr/bin/env python
"""
This script calculates the average between-person JSD for the same
site vs. the within-person JSD for different sites.
"""

import pandas as pd
import numpy as np

import argparse

p = argparse.ArgumentParser()
p.add_argument('fjsd', help='File with tidy JSD. Has site_comparison, '
    + 'subject, sample1, sample2, patient_comp, and beta columns')
p.add_argument('fmeta', help='metadata file path.')
p.add_argument('fout', help='output file path')
args = p.parse_args()

fout = args.fout

## Read in data
fjsdlong = args.fjsd
fnmeta = args.fmeta

jsd = pd.read_csv(fjsdlong, sep='\t')
meta = pd.read_csv(fnmeta, sep='\t', index_col=0)

# Keep just within-batch comparisons
jsd2014 = jsd.query('(batch1 == 2014) & (batch2 == 2014)')
jsd2016 = jsd.query('(batch1 == 2016) & (batch2 == 2016)')
jsd = pd.concat((jsd2014, jsd2016))

allcomps = ['bal-throat_swab', 'bal-gastric_fluid', 'gastric_fluid-throat_swab']

reslst = []

for subj, ptjsd in jsd.query('site_comparison == @allcomps').groupby('subject'):
    samples = list(set(ptjsd['sample1'].tolist() + ptjsd['sample2'].tolist()))
    for s in samples:
        site = meta.loc[s, 'site']

        # This has the between patients and same site comparison,
        # between the patient we're looking at and all others
        site_comp = site + '-' + site
        sitejsd = jsd\
            .query('patient_comp == "between"')\
            .query('site_comparison == @site_comp')\
            .query('(sample1 == @s) | (sample2 == @s)')

        # Track mean and median beta diversity between that patient
        # and all others (same site comparison)
        mean_btw = np.mean(sitejsd['beta'])
        median_btw = np.median(sitejsd['beta'])

        # Get the site comparisons that patient has that include that site
        # note: do this because I don't necessarily know the comparison string
        # wrt the sample we're iterating through now
        relevant_comps = [i for i in ptjsd['site_comparison'] if site in i]

        # Get the within-patient beta div for each of those comparisons
        for comp in relevant_comps:

            # get the site1-site2 within patient jsd
            beta = ptjsd.query('site_comparison == @comp')['beta'].values[0]

            # how many site-site jsd's are lower than this within-patient one?
            n_lower = sum(sitejsd['beta'] < beta)

            # Prepare data to store
            row = [subj, site, s, comp,
                   sitejsd.shape[0], n_lower,
                   beta, mean_btw, median_btw]

            reslst.append(row)

# Convert into dataframe
## Note: "n_lower" means "number of inter-patient comparisons
## (site1 vs site1) which were more similar than that patient's
## respective intra-patient comparison (site1 vs. site2)"
resdf = pd.DataFrame(
    reslst,
    columns=['subject', 'site', 'sample', 'site_comparison',
             'n_btw_patient', 'n_lower',
             'within_beta', 'mean_inter_beta', 'median_inter_beta'])
resdf['percent_lower'] = resdf['n_lower']/resdf['n_btw_patient'].astype(float)

# Write to file
resdf.to_csv(fout, sep='\t', index=False)
