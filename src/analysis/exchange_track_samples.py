#!/usr/bin/env python
"""
This script is essentialy a bare-bones copy of the script that
calculates the correlation and partial correlations of all
OTUs across sites.

Its goal is to track the samples I used in calculating the OTU
correlations (AKA their exchanged-ness).
"""
import argparse
import pandas as pd
import numpy as np

import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
import util



if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('fn_otu', help='OTU table file path. OTUs in columns, '
                                  + 'samples in rows. Relative abundances.')
    p.add_argument('fn_meta', help='metadata file path.')
    p.add_argument('fn_patients', help='stem of output file (i.e. the part'
        + ' of the file name before ".site.samples.txt"')
    # Add this argument so that we only track the samples we end up using
    # in the final paper
    p.add_argument('--npatients', help='number of patients used as the cutoff',
        default=4, type=int)
    p.add_argument('--shuffle', action='store_true')
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

    sites = util.get_sites()

    ## Initialize sample tracking
    all_samples = {}
    for site1 in sites:
        for site2 in sites[sites.index(site1)+1:]:
            all_samples[site1 + '-' + site2] = []

    for o, otudf in dflong.groupby('otu'):
        #print(o)
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

                if len(patients) >= args.npatients:
                    # Here, the real code calculates correlation between non-zero abun in site1 and site2

                    ## Next, it does the partial correlation conditioned on site 3
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

                    # If there are enough patients for this calculation to
                    # be considered:
                    if len(three_patients) > args.npatients:
                        # Track the patients and samples used
                        new_samples = otudf.query('subject_id == @three_patients')['sample'].tolist()
                        # Add to previous list
                        samples_so_far = all_samples[site1 + '-' + site2]
                        samples_so_far += new_samples
                        all_samples[site1 + '-' + site2] = list(set(samples_so_far))

    ## Write all the results to files
    for k in all_samples:
        fname = '.'.join([args.fn_patients, k, 'samples', 'txt'])
        with open(fname, 'w') as f:
            f.write('\n'.join(all_samples[k]))
