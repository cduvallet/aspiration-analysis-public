#!/usr/bin/env python
"""
This script reads in the raw OTU table and wrangled metadata,
and returns an OTU table with:
- relative abundances
- OTUs with more than 100 reads
- samples with more than 100 reads
- samples with both metadata and 16S data
"""
import argparse
import pandas as pd

import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
from util import raw2abun

def remove_shallow_smpls(df, n_reads):
    """
    Removes samples with fewer than n_reads from dataframe df.

    Parameters
    -----------
    df : pandas dataframe
        samples are in rows, OTUs are in columns
    n_reads : int
        minimum number of reads per sample for sample to be kept
    """

    total_reads = df.sum(axis=1)
    shallow_smpls = [smpl for smpl in total_reads.index \
                     if total_reads.loc[smpl] <= n_reads]
    df = df.drop(shallow_smpls)

    return df

def remove_shallow_otus(df, perc_samples=None, n_reads=None):
    """
    Removes OTUs which are present in fewer than 100*perc_samples
    percent of samples OR which have fewer than n_reads reads.

    Parameters
    ----------
    df : pandas dataframe
        Samples are in rows. OTUs are in columns.
    perc_samples : float
        min percent of samples that an OTU must be present in to not
        be thrown out.
    n_reads : int
        min number of reads an OTU must have in df to not be thrown
        out.

    Either perc_samples or n_reads must be specified. If both are specified,
    the perc_samples filtering is done first and then OTUs with fewer than
    n_reads total are thrown out.

    """
    if perc_samples is not None:
        presencemap = lambda x: 1 if x else 0
        otus_perc_present = df.applymap(presencemap).sum() / df.shape[0]
        keepotus = list(
            otus_perc_present[otus_perc_present > perc_samples].index)
        df = df[keepotus]

    if n_reads is not None:
        # Removes any OTUs with fewer than n_reads from the raw and abun dfs
        # samples are in rows and OTUs are in columns
        total_reads = df.sum(axis=0)
        shallow_col_indices = [i for i in range(len(total_reads.index)) \
                               if total_reads.iloc[i] < n_reads]
        shallow_otus = df.columns[shallow_col_indices]
        df = df.drop(shallow_otus, axis=1)

    return df

def read_rosen_data(fnotu, fnmeta, sample_reads):

    df = pd.read_csv(fnotu, sep='\t', index_col=0).T
    df.index = [i if i != '04-087-1g' else '04-087-1G' for i in df.index]

    meta = pd.read_csv(fnmeta, sep='\t', index_col=0)
    meta.index = [i if i != '04-087-1g' else '04-087-1G' for i in meta.index]

    # Add read depth to metadata
    meta['total_reads'] = df.sum(axis=1)

    # Remove OTUs with fewer than 100 reads over ALL samples
    df = remove_shallow_otus(df, n_reads=100)
    # Remove sample with fewer than 100 reads
    df = remove_shallow_smpls(df, n_reads=sample_reads)
    # Update metadata to reflect new OTU table
    meta = meta.loc[df.index]
    abundf = raw2abun(df)

    return df, abundf, meta

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('fnotu', help='path to raw OTU table')
    p.add_argument('fnmeta', help='path to wrangled metadata')
    p.add_argument('fnotu_out', help='path to write clean OTU table '
        + ' to (counts)')
    p.add_argument('fnabun_out', help='path to write clean OTU table '
        + ' to (relative abundance)')
    p.add_argument('fnmeta_out', help='path to write clean metadata to')
    # Optional arguments
    p.add_argument('--samplereads', help='minimum number of reads per sample.'
        + ' [default: %(default)s]', default=100, type=int)
    args = p.parse_args()

    df, abundf, meta = read_rosen_data(
        args.fnotu, args.fnmeta, args.samplereads)
    df.to_csv(args.fnotu_out, sep='\t')
    abundf.to_csv(args.fnabun_out, sep='\t')
    meta.to_csv(args.fnmeta_out, sep='\t')
