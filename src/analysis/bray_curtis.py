#!/usr/bin/env python
"""
This script calculates the bray-curtis distance between all sample pairs.
"""

import argparse
import pandas as pd
import numpy as np
from skbio.diversity import beta_diversity

import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
from util import add_metadata_to_beta_list

p = argparse.ArgumentParser()
p.add_argument('otu', help='path to OTU table [counts]')
p.add_argument('meta_in', help='path to metadata file, samples in rows')
p.add_argument('beta', help='path to write beta diversity to')
args = p.parse_args()

df = pd.read_csv(args.otu, sep='\t', index_col=0)

samples = list(df.index)

# Don't include samples from second time point or lung transplants
exclude = ['2', 'F', 'sick', 'F2T']
for s in exclude:
    samples = [i for i in samples if not i.endswith(s)]
samples = [i for i in samples if not i.startswith('05')]

df = df.loc[samples]

beta_div = beta_diversity('braycurtis', df.values, ids=df.index)

# Convert to redundant, wide format
widebeta = pd.DataFrame(beta_div.redundant_form(),
    index=df.index, columns=df.index)

# Save the wide-form dataframe
fname = args.beta.split('.txt')[0] + '.wide.txt'
widebeta.to_csv(fname, sep='\t')

## Convert to long-form, tidy format
widebeta.index.name = 'sample1'
# Remove redundant comparisons before converting to tidy.
# k=1 removes the diagonal, so we don't keep the self-self comparisons
betadf = widebeta.where(np.triu(np.ones(widebeta.shape, dtype=bool), k=1))

betadf = pd.melt(betadf.reset_index(),
    id_vars='sample1', var_name='sample2', value_name='beta')
betadf = betadf.dropna()

# Add metadata to each comparison
meta = pd.read_csv(args.meta_in, sep='\t', index_col=0)
aspcols = ['Results of MBS closest to enrollment date',
           'Results of worst MBS',
           'mbs_consolidated']
ppicols = ['ppi_consolidated', 'On PPI currently?', 'PPI Status',
           'Patient taking PPI', 'Patient taking PPI?', 'ACIDSUP']
refluxcols = ['total duration of acid reflux',
              'Reflux - total number of episodes',
              'Number of acid reflux episodes',
              'Total number of reflux episodes (acid+non-acid)',
              'SI - Reflux', 'SSI - Reflux', 'SAP - Reflux',
              'Number of non-acid reflux episodes',
              'percent distal nonacid', 'percent proximal total',
              'percent distal acid', 'percent proximal acid',
              'percent proximal nonacid', 'percent distal total',
              'number of full colum events/total events',
              'Number of full column episodes']
bilecols = ['Was Bile CA detected?',
            'Was Bile DCA detected?',
            'Was Bile LCA detected?',
            'Was Bile TCA detected?']
metacols = aspcols + ppicols + refluxcols + bilecols

betadf = add_metadata_to_beta_list(betadf.values.tolist(), meta, metacols)

betadf.to_csv(args.beta, sep='\t', index=False)
