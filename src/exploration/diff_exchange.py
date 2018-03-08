#!/usr/bin/env python
"""
This script prints the differential exchange between aspirators and
non-aspirators for each site.

Note that only OTUs with the correct min_patients, r_thresh, and q_thresh
(set in exchange.py) have values for differential exchange.
"""

import pandas as pd
import argparse
import os

p = argparse.ArgumentParser()
p.add_argument('fnexchange')
p.add_argument('outdir')

try:
    args = p.parse_args()
except:
    fnexchange = '~/github/aspiration-analysis/data/analysis/exchange.txt'
    outdir = '~/github/aspiration-analysis/src/exploration'
    args = p.parse_args([fnexchange])

os.chdir(args.outdir)

df = pd.read_csv(args.fnexchange, sep='\t')

df = df.dropna(subset=['diff_exchange'])

for g, subdf in df.groupby('site_comparison'):
    # Make sitedf, a dataframe with aspirators in columns and OTUs in rows
    sitedf = subdf.pivot(index='otu', columns='patient_type', values='exchange')
    subdf.index = subdf['otu']
    # Join on OTUs. subdf has duplicate rows bc each OTU has one row for
    # aspirators and one row for non-aspirators
    sitedf = sitedf.join(subdf['diff_exchange']).drop_duplicates()
    sitedf = sitedf.sort_values(by='diff_exchange', ascending=False)
    sitedf.to_csv('differential_exchange.{}.txt'.format(g), sep='\t')
    sitedf.to_excel('differential_exchange.{}.xls'.format(g))
