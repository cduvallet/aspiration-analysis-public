#!/usr/bin/env python
"""
This script calculates alpha diversity for all samples.
"""

"""
This script calculates the alpha diversity for each non-collapsed
OTU table.
"""
import numpy as np
import pandas as pd
import argparse

import skbio.diversity.alpha as alph

import os, sys
# Add src/util to path and import modules from files there
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
from util import get_metacols

def make_alpha_df(df, metric):
    """
    Make the alpha diversity tidy dataframe.

    Parameters
    ----------
    df: pandas DataFrame
        samples in rows, OTUs in columns
    meta : pandas DataFrame
        samples in rows, 'DiseaseState' column
    metric : str
        alpha diversity metric
    metacols : list
        list of metadata columns to add to alpha dataframe

    Returns
    -------
    alpha : pandas DataFrame
        tidy dataframe with ['sample', 'alpha', 'alpha_metric'] columns
    """
    alpha = alpha_diversity(df, metric).reset_index()
    alpha.columns = ['sample', 'alpha']
    alpha['alpha_metric'] = metric

    return alpha

def alpha_diversity(df, metric='shannon'):
    """
    Calculate shannon diversity of all samples in df.

    Parameters
    ----------
    df : pandas DataFarme
        dataframe with samples in rows, OTUs in columns
    metric : str
        'shannon', 'chao1', 'simpson'

    Returns
    -------
    alpha : pandas Series
        pandas series with samples in rows
    """

    alphafun = alph.shannon
    if metric == 'chao1':
        alphafun = alph.chao1
    elif metric == 'simpson':
        alphafun = alph.simpson
    elif metric != 'shannon':
        print('Unknown alpha diversity metric. Doing Shannon Index.')

    return df.apply(alphafun, axis=1)

p = argparse.ArgumentParser()
p.add_argument('counts_in', help='path to OTU table with raw counts, OTUs in '
    + 'columns, samples in rows.')
p.add_argument('meta', help='path to metadata file')
p.add_argument('alphas_out', help='out file with all alpha diversities for '
    + 'all samples')
args = p.parse_args()

## Read data
df = pd.read_csv(args.counts_in, sep='\t', index_col=0)
meta = pd.read_csv(args.meta, sep='\t', index_col=0)

# Metadata columns pre-defined in util module
metacols = get_metacols()

alphas = []
for metric in ['shannon', 'chao1', 'simpson']:
    alpha = make_alpha_df(df, metric)
    alphas.append(alpha)

alphasdf = pd.concat(alphas, ignore_index=True)

# Add sample metadata
alphasdf = pd.merge(meta[metacols], alphasdf,
                    left_index=True, right_on='sample',
                    how='right')

alphasdf.to_csv(args.alphas_out, sep='\t', index=False)
