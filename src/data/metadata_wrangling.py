#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 15:07:06 2017

@author: claire
"""
import os
import argparse
import pandas as pd
import numpy as np
from collections import Counter

def check_duplicate_samples_df(df, fname='', index_col=None):
    """
    Check if one df has duplicate samples.
    Returns boolean nodups, set to True if duplicate samples were found
    """
    if index_col is not None:
        smplcounter = Counter(df[index_col])
    else:
        smplcounter = Counter(df.index)
    dupsmpl = [i for i in smplcounter if smplcounter[i] > 1]
    if len(dupsmpl) > 0:
        print('File {} has duplicate samples'.format(fname) + '\n\t' + ', '.join(dupsmpl))
        nodups = False
    else:
        nodups = True
    return nodups, dupsmpl

def check_duplicate_samples_dict(metadata):
    """
    Check if any df's in metadata dict have duplicate samples

    metadata is a dict with {key: {'df': pandas dataframe, 'sample_col': str}}
    """
    nodups = True
    for m in metadata:
        df = metadata[m]['df']
        nodups, _ = check_duplicate_samples_df(df, m, metadata[m]['sample_col'])
    if nodups:
        print('No files have duplicate samples!')
    return None

def remove_duplicate_rows_from_k23aim4form3(metadata):
    """
    Fix and remove the duplicate samples in K23 Aim 4 Form 3
    """

    m = 'k23aim4form3'
    df = metadata[m]['df']

    ## Find the columns with values that don't match

    #s = '04-256-8'
    #tmp = df[ df[metadata[m]['sample_col']] == s].dropna(how='all', axis=1)
    #print('The following subset of columns have different metadata for sample {}'.format(s))
    #print(tmp.loc[:, tmp.apply(lambda col: True if col.iloc[0] != col.iloc[1] else False)])
    # Keep the more recent of the two rows
    df = df.drop(199)

    #s = '04-029-8'
    #tmp = df[ df[metadata[m]['sample_col']] == s].replace('n/a', np.nan).dropna(how='all', axis=1)
    #mismatches = tmp.loc[:, tmp.apply(lambda col: True if col.iloc[0] != col.iloc[1] else False)]
    #print('{} columns have different metadata for sample {}'.format(mismatches.shape[1], s))
    #mismatches.applymap(lambda x: 1 if isinstance(x, float) and np.isnan(x) else 0).sum(axis=1)
    # Just looking at the data, seems like the 29th row has 3 nan's, corresponding to the MBS results
    # Row 28 has those MBS results
    mbscols = ['Date of MBS closest to enrollment date',
               'Months between enrollment and MBS closest to enrollment date',
               'Results of MBS closest to enrollment date']
    df.loc[29, mbscols] = df.loc[28, mbscols]
    df = df.drop(28)

    return df

def concatenate_and_merge_rows(metadata):
    """
    Rename index according to sample_col,
    concatenate all the metadata files,
    and remove any duplicate rows by merging their non-overlapping columns.
    Returns a dataframe with subjects in the index.

    Does not:
        - handle situations where same subject has conflicting metadata for a column
        - same subject has more than 2 entries
    """
    # Relabel index of each df
    for m in metadata:
        df = metadata[m]['df']
        df.index = df[metadata[m]['sample_col']]
        metadata[m]['df'] = df

    bigdf = pd.concat((metadata[m]['df'] for m in metadata))
    _, dupsmpls = check_duplicate_samples_df(bigdf)
    ## Note: no samples have duplicate columns, and this takes a long time. Skip it
    ## See how many duplicate columns each sample has
    #for smpl in dupsmpls:
    #    dupcols = (bigdf.loc[smpl].apply(lambda col: len(col.dropna().unique())) > 1).sum()
    #    if dupcols > 0:
    #        print('Sample {} has {} duplicate cols'.format(smpl, dupcols))

    # Looks like none have conflicting columns! Just combine the rows using df.combine_first()
    print('Combining duplicate samples in big metadata (bc none of their columns are conflicting)...')
    for s in dupsmpls:
        subdf = bigdf.loc[s]
        if subdf.shape[0] > 2:
            print('Warning! Sample {} has {} entries! Only first 2 rows will be combined!'.format(s, subdf.shape[0]))
        subdf = subdf.iloc[0].combine_first(subdf.iloc[1])
        bigdf = bigdf.drop(s)
        bigdf.loc[s] = subdf

    return bigdf

def find_closest_match_subject(smpls2subj_map, meta):
    """
    Print the closest matching subject in metadata for samples in
    smpls2subj dict with 'unk' subject.
    """
    otu_no_meta = [i for i in smpls2subj_map if smpls2subj_map[i] not in meta.index]
    for s in set(otu_no_meta):
        print('{}: {}'.format(s, smpls2subj_map[s]))
        print('\t' + ', '.join([i for i in meta.index if smpls2subj_map[s][:6] in i]))
    return None

def read_ppi_metadata(datadir, fnppi):
    """
    Read in the ppiMetadata.csv file. Fix weird formatting and
    rename the one duplicate sample so that it matches its sample ID in
    OTU table
    """
    ppidf = pd.read_csv(os.path.realpath(os.path.join(datadir, fnppi)))
    ppidf['STUDYID'] = ppidf['STUDYID'].apply(lambda x: x.strip())
    ppidf['SAMPLEID'] = ppidf['SAMPLEID'].apply(lambda x: x.strip())
    # Edit the duplicate 05-013-1G sample
    ppidf.loc[173, 'SAMPLEID'] = 'dup05-013-1G'
    ppidf.loc[172, 'SAMPLEID'] = 'dup05-013-1B'
    return ppidf

def sample2subject(sample):
    """
    Return the subject corresponding to the given sample.

    sample is a string corresponding to the sample ID in the OTU table.

    The normal parsing is: <##>-<##>-<#S>, where S indicates the sample type
    (BAL, gastric, throat swab).

    There are also fundoplication samples, which can end in F1, F2, TI, and TF.
    These have one digit in their middle number (i.e. <##>-<#>-<F1>).
    """

    # First pass: sample ID belongs to subjects of the form <##>-<##>-<one letter/number>
    subj = '-'.join(sample.split('-')[0:2]) + '-' + sample.split('-')[2][0]

    # If sample ID middle number is only one digit, it's a fundo sample
    # and should be turned into <##>-<#>-F1 to match subject ID in metadata
    if len(subj.split('-')[1]) == 1:
        subj = subj.rsplit('-',1)[0] + '-F1'
        # These fundoplication samples are also missing the leading 0 in
        # the OTU table (but have it in the metadata)
        if len(subj.split('-')[0]) == 2:
            subj = '0' + subj

    ## If sample ends in TF or TI, it is also a fundoplication subject.
    ## Note: This part is commented out bc these samples also all have
    ## only one digit in their middle ID
    ## TODO: double check this with Rachel!!
    #if subj.split('-')[2] == "TI" or subj.split('-')[2] == "TF":
    #    subj = '-'.join(subj.split('-')[0:2]) + '-F1'

    # Special cases
    if sample in ['03-076-8-B-SF', '03-076-8-B-TF']:
        return '03-076-8'

    return subj

def sample2batch(sample, files2014, files2016):
    """
    Return which batch of data the sample came from.

    fn2014batch and fn2016batch are paths to files with lists of fastq
    files of the form sampleID.merged.fastq
    """
    if sample in files2014 and sample not in files2016:
        batch = '2014'
    elif sample in files2014 and sample in files2016:
        batch = 'both'
    elif sample in files2016 and sample not in files2014:
        batch = '2016'
    else:
        batch = 'unk'
    return batch

def sample2site(sample):
    """
    Return the site corresponding to give sample
    """
    # Define abbreviations to site dict
    d = {'B': 'bal', 'T': 'throat_swab', 'R': 'rectal_swab',
         'G': 'gastric_fluid', 'g': 'gastric_fluid',
         'S': 'stool'}
    # Standard case is that the samples are labeled <subject_ID><site>, where
    # <site> is either B, T, or G
    if sample[-1] in ['B', 'T', 'G', 'g', 'R']:
        return d[sample[-1]]

    # Another case is that the sample ID ends with 'TF', 'TI', 'SF', or 'SI'
    # 'TI' and 'TF' are initial and final throat swabs before and after PPI
    # 'SI' and 'SF' are initial and final stool samples before and after PPI
    # These are R01 samples
    elif sample[-2:] in ['TI', 'TF', 'SF', 'SI']:
        return d[sample[-2]]

    # The fundoplication samples are throat swabs, and end with 'F1' or 'F2'
    elif sample[-2:] in ['F1', 'F2']:
        return d['T']

    # And finally, there are some samples that end in 'Xsick' or 'XI' or 'XF',
    # where X is 'B', 'G', 'R', or 'T'. I'm not sure what study these samples
    # belong to, but I think we can safely assume they come from the respective
    # sites
    elif sample[-2:] in ['RI', 'GI', 'RF', 'GF']:
        return d[sample[-2]]
    elif sample.endswith('sick'):
        return d[sample[-5]]

    else:
        return 'unk'

def sample2baseline(sample):
    """
    Return whether the sample is a baseline sample or additional time point
    for that patient.

    For example, the K23 samples just have one sample per subject (per body
    site). On the other hand R01 samples have pre and post PPI therapy,
    with stool/rectal swab and throat samples for each patient. If we ever
    want to use these in an analysis where we need to have only one sample
    per patient, we need to keep only the baseline samples (so each patient
    doesn't have biological duplicates.)

    Returns 1 if baseline sample, 2 if follow-up sample.
    Doesn't handle more than 2 time points so far...
    """

    # The K23 samples are easy: they are all labeled <##>-<###>-<#S>, where S
    # is the body site. They also should only have one sample per patient
    splitsample = sample.split('-')
    if (len(splitsample) ==3 and len(splitsample[0]) == 2
            and len(splitsample[1]) == 3 and len(splitsample[2]) == 2):
        return '1'

    # The fundoplication samples are of the form: '56-3-TI' or '030-5-F1'.
    # There are also a few fundoplication samples which are like '024-1-F1T'
    # They all have a middle ID that is one digit
    if len(splitsample[1]) == 1:
        # Look at the second character in the third delimited ID, e.g.
        # 1 or 2 (for the F1/F2 samples) or I or F (for the TI/TF samples)
        if splitsample[2][1] == '1' or splitsample[2][1] == 'I':
            return '1'
        elif splitsample[2][1] == '2' or splitsample[2][1] == 'F':
            return '2'

    # The R01 samples have pre- and post-PPI therapy
    # The R01 samples are of the form '01-263-4RI', and all have three-character
    # third delimiters. Note that some of the fundo samples have this
    # too: '026-8-F2T'
    if len(splitsample[2]) == 3 and not splitsample[2].startswith('F'):
        if splitsample[2][-1] == 'I':
            return '1'
        elif splitsample[2][-1] == 'F':
            return '2'

    return 'unk'

def investigate_acidsup_ppicol_discrepancies(meta, ppidf, ppicols):
    """
    Look at the rows which have values in meta[ppicols] which conflict
    with value in ppidf['ACIDSUP'], where ppidf was read from the
    ppiMetadata.csv file and meta is the concatenation of all the Excel files.
    """
    # Look at rows which have ACIDSUP defined and at least one of the ppicols defined
    tmp = pd.merge(meta, ppidf, left_on='subject_id', right_on='STUDYID',
                   how='outer') \
                  [ppicols + ['ACIDSUP', 'STUDYID', 'SAMPLEID']] \
                  .dropna(subset=['ACIDSUP']) \
                  .dropna(subset=ppicols, how='all')

    # Find the rows where a ppicol value clashes with the ACIDSUP value
    ppimap = {'yes': 'On', 'no': 'Off',
              'off': 'Off', 'on': 'On',
              'On': 'On', 'Off': 'Off',
              'Yes': 'On', 'No': 'Off',
              'n/a': np.nan, np.nan: np.nan}

    def no_conflict(row, ppicols, ppimap):
        acidsup = row['ACIDSUP']
        others = row[ppicols].apply(lambda x: ppimap[x])
        if len(set(others.dropna().values)) == 1 and acidsup in others.values:
            return True
        else:
            return False
    tmp = tmp.loc[tmp.apply(lambda row: no_conflict(row, ppicols, ppimap), axis=1) == 0]
    return tmp

def consolidate_ppi(row, ppicols):
    """
    Return the consolidated value across all ppicols in row.

    If non-nan values do not conflict, return that value.
    If non-nan values do conflict, return 'conflicting'
    """
    ppimap = {i: 'on' for i in ['on', 'On', 'yes', 'Yes']}
    ppimap.update({i: 'off' for i in ['off', 'Off', 'No', 'no']})
    ppimap['n/a'] = 'n/a'

    if row[ppicols].isnull().sum() == len(ppicols):
        return np.nan

    tmp = row[ppicols].dropna().apply(lambda x: ppimap[x])

    if tmp.shape[0] == 1:
        return tmp.values[0]

    elif tmp.eq(tmp.iloc[0]).sum() == tmp.shape[0]:
        # If all values in the row are equal, return that value
        return tmp.iloc[0]
    else:
        return 'conflicting'

def consolidate_mbs(row, worstmbs, recentmbs):
    """
    Consolidate conflicting/missing MBS results

    If both a most recent and a worst MBS result are given, return the most
    recent.
    If only a worst or most recent MBS result is given, return that.

    Values of Aspiration or Penetration are converted to "Aspiration/Penetration"
    """
    def is_not_null(val):
        if isinstance(val, float) and not np.isfinite(val):
            return False
        else:
            return True

    aspmap = {'Normal': 'Normal', 'Aspiration': 'Aspiration/Penetration',
              'Penetration': 'Aspiration/Penetration', 'Other': 'Other'}

    # If recent is defined, return it. If it isn't, return worst.
    # If neither are defined, return nan.
    if is_not_null(row[recentmbs]):
        return aspmap[row[recentmbs]]
    elif is_not_null(row[worstmbs]):
        return aspmap[row[worstmbs]]
    else:
        return np.nan


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('datadir', help='directory with ppiMetadata.csv file and '
                   + '`june2016-emails` folder that contains the Excel files.')
    p.add_argument('fnotu', help='path to raw OTU table with samples in columns '
                   + 'and OTUs in rows. Only samples in this OTU table are '
                   + 'wrangled.')
    p.add_argument('fn2014batch', help='file with the *.merged.fastq files in '
                   + 'the 2014 batch of data., used to label samples with '
                   + 'which batch they came from.')
    p.add_argument('fn2016batch', help='file with the *.merged.fastq files in '
                   + 'the 2016 batch of data, used to label samples with '
                   + 'which batch they came from.')
    p.add_argument('fn_metaout', help='path to output metadata file')
    args = p.parse_args()

#    datadir = '../../data/raw/metadata'
#    fnotu = '/Users/claire/github/disease_dataset_project/data/reflux_rosen/rosen_mincount10_maxee2_trim200_results/RDP/rosen_mincount10_maxee2_trim200.otu_table.99.denovo.rdp_assigned'
#    fn2014batch = '/Users/claire/github/disease_dataset_project/data/reflux_rosen/rosen_mincount10_maxee2_trim200_results/2014_batch.txt'
#    fn2016batch = '/Users/claire/github/disease_dataset_project/data/reflux_rosen/rosen_mincount10_maxee2_trim200_results/2016_batch.txt'
#    fn_metaout = '../../data/clean/all_metadata.txt'
    datadir = args.datadir
    fnotu = args.fnotu
    fn2014batch = args.fn2014batch
    fn2016batch = args.fn2016batch
    fn_metaout = args.fn_metaout


    ## TODO: add this later!
    # This will be its own separate business
    metabfile = ['CHHB-01-15VW CLIENT DATA TABLE.XLSX']

    #### PPI STUDY METADATA FILE ####
    ## Set up metadata file info
    fnppi = 'ppiMetadata.csv'
    ppidf = read_ppi_metadata(datadir, fnppi)

    #### EXCEL FILES ####
    print('Reading in the excel files...')
    # Label the metadata files
    names = ['fundo_baseline', 'k23aim4form3', 'r01form2', 'r01form3aims12', 'k23aim2form3',
             'k23aim5form3', 'r01form3aim3']
    excelfiles = ['june2016-emails/Fundo Baseline_6.2.16.xlsx',
                  'june2016-emails/K23 Aim 4 Form 3 Procedure Results_6.2.16.xlsx',
                  'june2016-emails/R01 Form 2 Baseline Aim 3.xlsx',
                  'june2016-emails/R01 Form 3 Aims 1 and  2.xlsx',
                  'june2016-emails/K23 Aim 2 and 3 Form 3 Procedure Results_6.2.16.xlsx',
                  'june2016-emails/K23 Aim 5 Form 3 Procedure Results_6.2.16.xlsx',
                  'june2016-emails/R01 Form 3 Aim 3.xlsx']
    excelfiles = [os.path.realpath(os.path.join(datadir, i)) for i in excelfiles]
    # Column with sample ID
    samplecols = ['Subject ID', 'Subject ID number', 'A1. Subject ID number:',
                  'Subject ID', 'Subject ID number', 'Subject ID number',
                  'Subject ID number']

    # Make dict to hold all the metadata info
    metadata = {n: {'excel_file': e, 'sample_col': s} for n, e, s in zip(names, excelfiles, samplecols)}

    ## Read in the excel files
    for m in metadata:
        metadata[m]['df'] = pd.read_excel(metadata[m]['excel_file'], na_values=' ') \
                .dropna(how='all', axis=1).dropna(how='all', axis=0) \
                .dropna(subset=[metadata[m]['sample_col']]) \
                .drop_duplicates()

    ## Clean up duplicate samples
    check_duplicate_samples_dict(metadata)
    # Fix the duplicate samples in the one file with 2 duplicate samples
    metadata['k23aim4form3']['df'] = remove_duplicate_rows_from_k23aim4form3(metadata)
    check_duplicate_samples_dict(metadata)

    ## Concat all the metadatas
    meta = concatenate_and_merge_rows(metadata)
    check_duplicate_samples_df(meta)

    # And add the subject IDs (which are in index) to own column
    meta['subject_id'] = meta.index

    ### Next up: grab the subjects with 16S data, add sample-wise metadata
    ### (i.e. data batch, sample site)

    #### OTU TABLE - INFO FROM SAMPLE IDS ####
    print('Getting sample IDs and metadata from sample IDs in OTU table...')
    ## Read in OTU table
    otu = pd.read_csv(fnotu, sep='\t', index_col=0).T

    # Map samples in OTU table to subject IDs in metadata
    smpls2subj_map = {o: sample2subject(o) for o in otu.index}
    otuinfo = pd.DataFrame(smpls2subj_map.items(), columns=['sample_id', 'subject_id'])

    # If you feel like it, look at the nearest subjects in meta for the samples without matches
    if False:
        find_closest_match_subject(smpls2subj_map, meta)

    # Make dicts to map sample to batch and body site
    files2014 = [i.split('.merged.fastq')[0] for i in open(fn2014batch, 'r').readlines()]
    files2016 = [i.split('.merged.fastq')[0] for i in open(fn2016batch, 'r').readlines()]

    otuinfo['batch'] = otuinfo['sample_id'].apply(lambda x: sample2batch(x, files2014, files2016))
    otuinfo['site'] = otuinfo['sample_id'].apply(sample2site)
    otuinfo['sample_number'] = otuinfo['sample_id'].apply(sample2baseline)


    #### MERGE ALL METADATA ####
    print('Merging all metadata sources...')
    # Right merge the subject-wise massive metadata (meta) on the subjects in the
    # OTU table (otuinfo). Do right merge so that we keep all samples with 16S.
    # Note: There are 260 unique subjects in the OTU table.
    #           otuinfo['subject_id'].unique().shape
    #       And 713 unique subjects in the massive metadata
    #           meta['subject_id'].unique().shape
    #       But only 241 subjects with both 16S data and metadata:
    #           allmeta = pd.merge(meta, otuinfo, left_on='subject_id',
    #                              right_on='subject_id', how='inner')\
    #                              ['subject_id'].unique().shape
    if 'sample_id.1' in meta.columns:
        print('meta, line 476, has sample_id.1')
    if 'sample_id.1' in otuinfo.columns:
        print('otuinfo, line 476, has sample_id.1')
    allmeta = pd.merge(meta, otuinfo, left_on='subject_id', right_on='subject_id',
                       how='right')
    if 'sample_id.1' in allmeta.columns:
        print('sample_id.1 in allmeta after merging meta with otuinfo')

    # Merge this overall metadata with the PPI metadata, keeping only samples
    # from the already merged df (i.e. only samples with 16S data)
    allmeta = pd.merge(allmeta, ppidf, left_on='sample_id', right_on='SAMPLEID',
                       how='left')
    if 'sample_id.1' in allmeta.columns:
        print('sample_id.1 in allmeta after merging with ppidf')
    #### CLEAN PPI AND MBS METADATA ####
    print('Consolidating PPI and MBS metadata...')
    ## PPI columns
    # There are 4 columns with PPI information in the Excel files.
    # Looks like 'On PPI currently?' and 'PPI Status' belong to the fundoplication
    # samples, and contain conflicting information in a few cases (i.e. 'PPI Status'
    # is 'on' but 'On PPI currently?' is 'No'). We'll need to ask Rachel for how
    # to consider these.
    # The other two columns, 'Patient taking PPI' and 'Patient taking PPI?' have
    # non-overlapping information (i.e. if a sample has a value in one column,
    # it does not have it in the other).
    # Finally, the ppiMetadata.csv file has an ACIDSUP column. We'll need to see
    # if this contains conflicting information as well.
    ppicols = ['On PPI currently?', 'PPI Status', 'Patient taking PPI',
               'Patient taking PPI?', 'ACIDSUP']
    ## MBS columns
    worstmbs = 'Results of worst MBS'
    recentmbs = 'Results of MBS closest to enrollment date'

    if False:
        # This function returns a dataframe with the rows with discrepancies
        investigate_acidsup_ppicol_discrepancies(meta, ppidf, ppicols)

    allmeta['ppi_consolidated'] = allmeta.apply(lambda row: consolidate_ppi(row, ppicols), axis=1)
    allmeta['mbs_consolidated'] = allmeta.apply(lambda row: consolidate_mbs(row, worstmbs, recentmbs), axis=1)

    allmeta.index = allmeta['sample_id']
    # Some of the column labels have newlines in them
    allmeta.columns = [i.replace('\n', '') for i in allmeta.columns]
    # Also at least one column label has a weird left-to-right mark (u'\u200e') in it...
    allmeta.columns = [i.encode('utf-8') for i in allmeta.columns]

    # There's one patient with a value for recentmbs but not worstmbs
    s = '02-300-2'
    for site in ['T', 'G']:
        allmeta.loc[s+site, worstmbs] = allmeta.loc[s+site, recentmbs]

    # And also one with 'Other' for recent mbs and 'Aspiration' for worst
    s = '05-176-4T'
    allmeta.loc[s, 'mbs_consolidated'] = 'Aspiration/Penetration'

    # Remove any empty columns
    allmeta = allmeta.dropna(how='all', axis=1)
    allmeta.to_csv(fn_metaout, sep='\t')
