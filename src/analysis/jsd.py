#!/usr/bin/env python
"""
This script takes in an OTU table and calculates JSD between all samples.

It writes a file with the following columns:
    sample1, sample2: which samples JSD is calculated between
    jsd : Jensen-Shannon distance (square root of divergence)
    patient_comp : 'within' or 'between', indicating if comparison is
                   between two patients or within the same subject
    site_comp : 'bal_throat', 'bal_gastric', 'throat_gastric', 'bal_bal',
               'gastric_gastric', 'throat_throat'
As well as respective aspiration, ppi, and reflux metadata columns (only
filled in for 'within' patient comparisons)
"""
import argparse
import pandas as pd
import numpy as np
import warnings
import sys

import os
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
from util import add_metadata_to_beta_list

def jsd(x,y):
    """
    Calculate Jensen-Shannon distance (square root of Jensen-Shannon
    divergence).

    Code from Thomas Gurry, who go it from somewhere else I think.

    Parameters
    ----------
    x, y : array-like
        vectors of relative abundance for OTUs, in each sample being
        compared
    """
    x = np.array(x)
    y = np.array(y)
    d1 = x*np.log2(2*x/(x+y))
    d2 = y*np.log2(2*y/(x+y))
    d1[np.isnan(d1)] = 0
    d2[np.isnan(d2)] = 0
    d = 0.5*np.sum(d1+d2)
    return np.sqrt(d)

p = argparse.ArgumentParser()
p.add_argument('otu_in', help='path to OTU table with relative abundances, '
                              + 'samples in rows, OTUs in columns.')
p.add_argument('meta_in', help='path to metadata file, samples in rows')
p.add_argument('jsd_out', help='path to write JSD values to')
args = p.parse_args()

## Read in files
df = pd.read_csv(args.otu_in, sep='\t', index_col=0).astype(float)
meta = pd.read_csv(args.meta_in, sep='\t', index_col=0)

# Set up metadata labels of interest
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
# Tranpose OTU table so that operations are performed column-wise
df = df.T

# Calculate JSD for all pairwise samples
samples = list(df.columns)

# Don't include samples from second time point or lung transplants
exclude = ['2', 'F', 'sick', 'F2T']
for s in exclude:
    samples = [i for i in samples if not i.endswith(s)]
samples = [i for i in samples if not i.startswith('05')]

# Update df with the right samples
df = df[samples]

jsd_lst = []
print("Calculating JSD...")
for i in range(0, len(samples)):
    if i%20 == 0:
        sys.stdout.write('{} '.format(i))
    for j in range(i+1, len(samples)):
        # Division by zero in np.log2() yields RunTimeWarning - ignore it
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            jsdval = jsd(df.iloc[:, i], df.iloc[:, j])
        jsd_lst.append([samples[i], samples[j], jsdval])

dfjsd = add_metadata_to_beta_list(jsd_lst, meta, metacols)

dfjsd.to_csv(args.jsd_out, sep='\t', index=False)

## Note: something about how I made my tidy dataframe is giving
## me headaches when I try to pivot it to a wide form dataframe.
## Since it's feasible, let's just make the wide dataframe here
## manually. Womp, tidy data: 1, Claire: 0.
samples = set(
    dfjsd['sample1'].unique().tolist() + dfjsd['sample2'].unique().tolist())
widejsd = pd.DataFrame(index=samples, columns=samples)
for g, subdf in dfjsd.groupby('sample1'):
    widejsd.loc[subdf['sample2'], g] = subdf['beta'].values
    widejsd.loc[g, subdf['sample2']] = subdf['beta'].values
widejsd = widejsd.fillna(0)
fname = args.jsd_out.split('.txt')[0] + '.wide.txt'
widejsd.to_csv(fname, sep='\t')

## Tidify data (with respect to reflux and bile metadata)
# Need to replace the 1's in 'number of full colum events/total events'
# with nan's
dfjsd['number of full colum events/total events'] = \
    dfjsd['number of full colum events/total events'].replace('1.0', np.nan)

## Turn dfjsd into tidy-format with reflux_type and bile_type aggregated
# Melt the bile acid metadata into variable, value pairs
tidybilejsd = pd.melt(dfjsd, value_vars=bilecols, var_name='bile_type',
                  value_name='bile_value',
                  id_vars=[i for i in dfjsd.columns if i not in bilecols])
# Remove any entries without bile measurement (or for between-patient comps)
tidybilejsd = tidybilejsd.dropna(subset=['bile_value'])
# Just keep the bile acid name
tidybilejsd['bile_type'] = tidybilejsd['bile_type'].apply(
    lambda x: x.split('Was Bile')[1].split('detected?')[0].strip())

# Melt the reflux metadata into variable, value pairs
tidyrefluxjsd = pd.melt(dfjsd, value_vars=refluxcols, var_name='reflux_type',
                  value_name='reflux_value',
                  id_vars=[i for i in dfjsd.columns if i not in refluxcols])
# Remove nan's
tidyrefluxjsd['reflux_value'] = \
    tidyrefluxjsd['reflux_value'].replace('n/a', np.nan)
tidyrefluxjsd = tidyrefluxjsd.dropna(subset=['reflux_value'])

# Merge them together
tidyjsd = pd.merge(tidybilejsd, tidyrefluxjsd,
                   on=[i for i in tidybilejsd.columns if i in
                       tidyrefluxjsd.columns],
                   how='outer')
tidyjsd = tidyjsd[[i for i in tidyjsd.columns if i
                   not in bilecols + refluxcols]]

## Write files
tidyjsd.to_csv(args.jsd_out.rsplit('.', 1)[0] + '.tidy_reflux_bile.txt',
               sep='\t', index=False)
print("Calculating JSD... Finished.")
