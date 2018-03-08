#!/usr/bin/env python
"""
This script calculates the weighted unifrac distance between all samples
in an OTU table.
"""
import argparse
import pandas as pd
import numpy as np
from skbio import TreeNode
from skbio.diversity import beta_diversity

import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
from util import add_metadata_to_beta_list


p = argparse.ArgumentParser()
p.add_argument('otu', help='path to OTU table [counts]')
p.add_argument('tree', help='path to tree with "denovo#" IDs [newick]')
p.add_argument('meta_in', help='path to metadata file, samples in rows')
p.add_argument('beta', help='path to write beta diversity to')
p.add_argument('--weighted', help='flag to calculate weighted unifrac',
    action='store_true')
args = p.parse_args()

df = pd.read_csv(args.otu, sep='\t', index_col=0)
tree = TreeNode.read(args.tree).root_at_midpoint()

# Get OTU IDs, just the "denovo#" part (to match the tree)
otu_ids = [i.split(';')[-1][3:] for i in df.columns]

# weighted or unweighted?
if args.weighted:
    metric = 'weighted_unifrac'
    kwargs = {'normalized': True}
else:
    metric = 'unweighted_unifrac'
    kwargs = {}

# Add tree and otu_ids to kwargs
kwargs.update({'tree': tree, 'otu_ids': otu_ids})

beta_div = beta_diversity(metric, df.values, ids=df.index,  **kwargs)

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

betadf = pd.melt(widebeta.reset_index(),
    id_vars='sample1', var_name='sample2', value_name='beta')

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
