#!/usr/bin/env python
"""
This script takes in a file with the between-site correlations and some
threshold values for correlation and q-values to identify and label
which OTUs are exchanged across which sites.
"""
import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
import argparse


p = argparse.ArgumentParser()
p.add_argument('exchangein', help='file with between-site correlations. '
    + 'Should have "n_partial", "r_partial", "p_partial", "otu", "site1", '
    + 'and "site2" columns.')
p.add_argument('exchangeout', help='file to write labeled exchanged '
    + 'OTUs to. This file will have each site_comparison as a column, '
    + 'OTUs in rows, and 1 for exchanged status (NaN if not exchanged).')
p.add_argument('--nthresh', help='minimum number of patients that an OTU '
    + 'must be non-zero in in both sites to be considered for exchange.',
    type=int)
p.add_argument('--qthresh', help='qvalue threshold to be considered '
    + 'exchanged. Multiple-test correction is performed after selecting '
    + 'OTUs which pass the nthresh patient threshold.', type=float)
args = p.parse_args()

df = pd.read_csv(args.exchangein, sep='\t')

df['site_comparison'] = df['site1'] + '-' + df['site2']

# Select subset of OTUs that meet N patient criteria
# i.e. these OTUs were non-zero in both sites of at least n_thresh patients
df = df.query('n_partial >= @args.nthresh')

# Correct for multiple tests among these OTUs
# Note: this isn't stratified by site comparison or anything. Is that bad?
# Should I be correcting across sites separately? I don't think so...
_, df['q_partial'], _, _ = multipletests(df['p_partial'], method='fdr_bh')

# Get the exchanged OTUs
# Note: we're no longer using a correlation threshold, bc all r values with
# q < qthresh should automatically be greater than zero.
exchange = df\
    .query('n_partial >= @args.nthresh')\
    .query('q_partial < @args.qthresh')

#print(exchange.groupby('site_comparison').size())

# Convert the tidy dataframe to longform with site_comparison in columns
# and 1 for exchanged OTUs, NaN for not exchanged OTUs
exchange = exchange\
    .pivot(columns='site_comparison', index='otu', values='r_partial')\
    .notnull().astype(int)\
    .replace(0, np.nan)

#print(exchange.sum())

# Write
exchange.to_csv(args.exchangeout, sep='\t')
