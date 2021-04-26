'''
Created on 1 Apr 2021

@author: Erik Alm (ECDC)

A variant is defined as:
 All sequences with the same lineage and the same set of significant mutations (defined by the scores in mutscores.csv and the SCORE_CUTOFF, see config)

Scoring rules for variants are:
 1) Sum of all mutation scores, and defined in the file mutscores.csv, positions not listed use a default score (see config)
 2) Log10 of number of total sequences for the variant
 3) 2.5 * The fraction of sequences of the variant that are newly submitted + 2.5 * The fraction of samples for the variant that are newly collected
The three scores are summed, there are weights available in the config

Indata (GISAID daily metadata dump) is expected in the INDIR as defined in the config. Outputs will be stored in csv format in OUTDIR.

'''

import tarfile
import pandas as pd
from glob import glob
import os
import math
from collections import Counter
import json

with open('config.json') as json_file:
    cfg = json.load(json_file)

INDIR = cfg['INDIR'] # Directory where GISAID EpiCoV metadata and definition files are located 
OUTDIR = cfg['OUTDIR'] # Directory for output

LOCATION_FIELD = cfg['LOCATION_FIELD'] # GISAID EpiCoV field for location
STRAIN_FIELD = cfg['STRAIN_FIELD'] # GISAID EpiCoV field for virus name
MUTATIONS_FIELD = cfg['MUTATIONS_FIELD'] # GISAID EpiCoV field for mutations
LINEAGE_FIELD = cfg['LINEAGE_FIELD'] # GISAID EpiCoV field for Pango lineage
COL_DATE_FIELD = cfg['COL_DATE_FIELD'] # GISAID EpiCoV field for collection date
SUB_DATE_FIELD = cfg['SUB_DATE_FIELD'] # GISAID EpiCoV field for submission date
DAYS_NEW_UPLOAD = cfg['DAYS_NEW_UPLOAD'] # Days since upload to count as a newly uploaded sequences
DAYS_NEW_COLLECTION = cfg['DAYS_NEW_COLLECTION'] # Days since collection to count as a newly collected sample
DAYS_NEW_VARIANT = cfg['DAYS_NEW_VARIANT'] # Days from first detection that counts as an early detection
MIN_SEQUENCES = cfg['MIN_SEQUENCES'] # Minimum number of sequences to display a variant, using 1 is very noisy
SCORE_CUTOFF = cfg['SCORE_CUTOFF'] # The minimum mutation score for a mutation to be included in the mutation label
DEFAULT_MUTATION_SCORE = cfg['DEFAULT_MUTATION_SCORE'] # The mutation score for unlisted genomic regions in the score definitions

WEIGHT_MUTATION = cfg['WEIGHT_MUTATION']
WEIGHT_NUMBER = cfg['WEIGHT_NUMBER']
WEIGHT_TREND = cfg['WEIGHT_TREND']

METADATA_FILE_PATTERN = cfg['METADATA_FILE_PATTERN']


''' Method for reading GISAID EpiCoV xz-compressed metadata file (published daily by GISAID) '''
def read_gisaid(infile):
    
    ''' Extract metadata file and read to dataframe '''
    with tarfile.open(infile, "r:xz") as tar:
        csv_path = tar.getnames()[1]
        df = pd.read_csv(tar.extractfile(csv_path), sep="\t", index_col=STRAIN_FIELD, dtype={LINEAGE_FIELD:str})
        df = df[~df.index.duplicated(keep='first')]
    
    
    
    ''' Split location field into parts, and change country to title case '''
    df[['region', 'country', 'division']] = df[LOCATION_FIELD].str.split('\s?/\s?', 2, expand=True)
    df['country'] = df['country'].str.title()
    df.loc[df['country']=='Usa', 'country'] = 'USA'
    
    
    ''' Convert dates to datetime and check whether each entry is newly collected and/or uploaded'''
    df['dt_submitted'] = pd.to_datetime(df[SUB_DATE_FIELD], errors='coerce', format='%Y-%m-%d')
    df['dt_collected'] = pd.to_datetime(df[COL_DATE_FIELD], errors='coerce', format='%Y-%m-%d')
    now = pd.to_datetime('now')
    df['newly uploaded'] = df['dt_submitted'].between(now - pd.Timedelta(DAYS_NEW_UPLOAD, 'd'), now)
    df['newly collected'] = df['dt_collected'].between(now - pd.Timedelta(DAYS_NEW_COLLECTION, 'd'), now)
    
    print(df)
    return df

''' Create a readable list of countries, sorted so that the ones with highest number of sequences come first '''
def agg_countrylist(country_list):
    counts = Counter(country_list).items()
    sorted_counts = sorted(counts, key=lambda x: x[1], reverse=True)
    result = ', '.join(c + '(' + str(v) + ')' for c, v in sorted_counts if c != '')
    return result


''' Method for extracting all the mutation information from a dataframe containing GISAID EPiCoV metadata '''
def extract_mutations(df, df_scores):
    
    ''' Extract mutation ref-base, position, and alt-base from metadata using regex '''
    pattern = r'(?P<mutation>[A-Za-z0-9]+_[A-Za-z0-9]+)'
    df_mutations = df[MUTATIONS_FIELD].str.extractall(pattern)
    df_mutations[['gene', 'sub']] = df_mutations['mutation'].str.split('_', 1, expand=True)
    pattern2 = r'(?P<ref>[A-Za-z]+)(?P<pos>[0-9]+)(?P<alt>[A-Za-z]+)'
    df_mutations[['ref', 'pos', 'alt']] = df_mutations['mutation'].str.extract(pattern2, expand=True)
    
    ''' Output any failed parsings, there should be none '''
    print('List of mutations with failed parsing:')
    print(df_mutations[df_mutations['pos'].isnull()])
    
    ''' Ensure positions are integer '''
    df_mutations['pos'] = df_mutations['pos'].astype(int)
    
    ''' Merge mutations with score table, set unlisted genomic positions to default score '''
    df_mutations = df_mutations.reset_index().merge(df_scores, how='left', on=['gene', 'pos']).set_index(STRAIN_FIELD)
    df_mutations.loc[df_mutations['score'].isnull(), 'score'] = DEFAULT_MUTATION_SCORE
    
    ''' Only include mutations scoring higher than a cut-off in the mutation profile used to label the variant '''
    df_mutations['mutation label'] = df_mutations['mutation']
    df_mutations.loc[df_mutations['score']<SCORE_CUTOFF, 'mutation label'] = ''
    df_mutations['mutation label'] = df_mutations['mutation label'].astype(str)
    
    return df_mutations

''' Method for calculating labels from a dataframe of pre-assigned labels and raw mutation and lineage labels dataframe '''
def calc_labels(df_variants, df_assigned):
    labels = {}
    status = {}
    comment = {}
    
    ''' Split mutation label into a list of mutations '''
    df_assigned['mutation_list'] = df_assigned['mutation label'].apply(lambda x: str(x).split(';'))
    df_variants['mutation_list'] = df_variants['mutation label'].apply(lambda x: str(x).split(';'))
    
    ''' Iterate over detected variants '''
    for index1, row in df_variants.iterrows():
        label = ''
        labels[index1] = label
        status[index1] = ''
        comment[index1] = ''
        min_mismatch = 100
        mismatch = 0
        
        ''' Iterate over pre-defined list of known variants of interest and concern '''
        for index2, row_ref in df_assigned.iterrows():
            
            ''' If lineage does not match, do not assign this as the label '''
            if row[LINEAGE_FIELD] != row_ref[LINEAGE_FIELD]:
                continue
            
            ''' Calculate which mutations that are missing from the pre-defined list variant '''
            missing_muts = []
            for m in row_ref['mutation_list']:
                if m not in row['mutation_list']:
                    mismatch += 1
                    missing_muts.append(m)
            
            ''' Calculate which additional mutations the variant has '''
            extra_muts = []
            for m in row['mutation_list']:
                if m!='nan' and m!='' and m not in row_ref['mutation_list']:
                    mismatch += 1
                    extra_muts.append(m)
            
            if mismatch < min_mismatch:
                min_mismatch = mismatch
            else:
                continue
            
                
            ''' If all characteristic mutations are found, assign VOC/VOI label '''
            status[index1] = ''
            if not missing_muts:
                status[index1] = row_ref['Status']
            elif row_ref['definition_flexible'] == 'no':
                ''' If the variant definition is not listed as flexible, do not assign label if characteristic mutations are missing'''
                continue
            
            ''' Assign this lineage as a label '''
            label = row_ref['Label']
            
            ''' List missing mutations (if any) after label '''
            for m in missing_muts:
                label = label + '-' + m.split('_')[1]
            
            ''' List additional mutations (if any) after label '''
            for m in extra_muts:
                label = label + '+' + m.split('_')[1]
            labels[index1] = label
           
            ''' If exact match, also assign the comment '''
            
            if not missing_muts and not extra_muts:
                comment[index1] = row_ref['Comment']
            else:
                comment[index1] = ''    
            
    
    ''' Apply the labels, statuses and comments to all variants '''
    df_variants['Label'] = df_variants.index.map(labels)
    df_variants['Status'] = df_variants.index.map(status)
    df_variants['Comment'] = df_variants.index.map(comment)
    
    return df_variants


if __name__ == '__main__':
    ''' Find and read latest GISAID EpiCoV snapshot file '''
    files = glob(os.path.join(INDIR, METADATA_FILE_PATTERN))
    files.sort()
    infile = files[-1]
    snapshot_date = '-'.join(os.path.basename(infile).split('_')[2:5]).split('.')[0]
    print('Latest snapshot found: ' + snapshot_date)
    df = read_gisaid(infile)
    
    ''' Read file containing already assigned variants '''
    df_assigned = pd.read_csv('assigned_variants.csv', sep=',', dtype={'mutation label':str})
    df_assigned.set_index("Label")
    print('Loaded pre-defiend variants:')
    print(df_assigned)
    
    ''' Read file containing scores for genomic positions '''
    df_scores = pd.read_csv('mutscores.csv', sep=",", dtype={'score':'float32'})
    print('Loaded mutation scores:')
    print(df_scores)
    

    ''' Extract and score mutations from GISAID EpiCoV snapshot  '''
    df_mutations = extract_mutations(df, df_scores)
    
    ''' Aggregate mutation data back to strains, creating a mutation label based on the mutations
        Adding up scores for each individual mutation
        Join with GISAID EpiCoV metadata '''
    df_viruses = df_mutations.groupby(by=STRAIN_FIELD, as_index=True).agg(
        {'mutation label': lambda x: ';'.join(a for a in sorted([str(i) for i in x]) if a and a != 'nan'),
         'mutation': lambda x: ','.join(a for a in sorted([str(i) for i in x]) if a and a != 'nan'),
         'score':'sum'})
    df_viruses = df_viruses.reset_index().merge(df, how='left', on=[STRAIN_FIELD]).set_index(STRAIN_FIELD)
    
    ''' Ensure that columns that will be aggregated over are string datatype '''
    df_viruses[LINEAGE_FIELD] = df_viruses[LINEAGE_FIELD].astype(str)
    df_viruses['country'] = df_viruses['country'].astype(str)
    
    ''' Get earliest collection date for each variant '''
    df_earliest_dates = df_viruses[[LINEAGE_FIELD, 'mutation label', 'dt_collected']].reset_index().groupby(
                                by=['mutation label', LINEAGE_FIELD]).agg({'dt_collected' : 'min'})
    
    ''' Label viruses that are newly collected or uploaded '''
    df_viruses['new countries uploaded'] = df_viruses.apply(lambda x : x['country'] if x['newly uploaded'] else '', axis=1)
    df_viruses['new countries collected'] = df_viruses.apply(lambda x : x['country'] if x['newly collected'] else '', axis=1)
    
    ''' Label viruses that were detected early, as defined in the config '''
    df_viruses['earliest countries'] = df_viruses.apply(
        lambda x : x['country'] if x['dt_collected'] <=\
            df_earliest_dates.loc[(x['mutation label'], x[LINEAGE_FIELD]), 'dt_collected'] + pd.Timedelta(DAYS_NEW_VARIANT, 'd') else '', axis=1)
    
    ''' Aggregate strains into variants based on the Pango lineage and mutation label '''
    df_groups = df_viruses[['score', LINEAGE_FIELD, 'mutation label', 'country', 'mutation', 'newly collected', 'newly uploaded',
                            'new countries uploaded', 'new countries collected', 'dt_collected', 'earliest countries']].reset_index().groupby(
                                by=['mutation label', LINEAGE_FIELD])
    
    df_variants = df_groups.agg({STRAIN_FIELD : 'count',
                                'score' : 'median',
                                'country' : lambda x : agg_countrylist(x),
                                'dt_collected' : 'min',
                                'earliest countries' : lambda x : agg_countrylist(x),
                                'new countries uploaded' : lambda x : agg_countrylist(x),
                                'new countries collected' : lambda x : agg_countrylist(x),
                                'mutation' : lambda x: x.value_counts().index[0],
                                'newly collected' : 'sum',
                                'newly uploaded' : 'sum',}).reset_index()
    
    ''' Rename some aggregated columns to more suitable labels '''    
    df_variants = df_variants.rename(columns={'Virus name': 'count','mutation': 'most common mutation profile',
                                      'score':'mutation profile score', 'country':'countries', 'dt_collected':'earliest collection date'})
    
    ''' Filter variants by the minimum number of sequences required '''
    df_variants = df_variants[df_variants['count']>=MIN_SEQUENCES]
    
    ''' Calculate an overall score based on mutation score number of sequences and the trend '''
    df_variants['number score'] = df_variants.apply(lambda x: math.log10(x['count']), axis=1)
    df_variants['trend score'] = df_variants.apply(lambda x: 2.5 * x['newly collected'] / x['count'] + 2.5 * x['newly uploaded'] / x['count'], axis=1)
    df_variants['overall score'] = df_variants.apply(
        lambda x: WEIGHT_MUTATION * x['mutation profile score'] + WEIGHT_NUMBER * x['number score'] + WEIGHT_TREND * x['trend score'], axis=1)
    
    ''' Sort variants by overall score '''
    df_variants = df_variants.sort_values(by='overall score', ascending=False)
    
    ''' Calculate labels for the variants '''
    df_variants = calc_labels(df_variants, df_assigned) 
    
    print(df_variants.columns)
    ''' Make sure the most important columns are displayed first '''
    prefix_cols = ['Label', 'Status', 'Comment', 'mutation label', 'Pango lineage', 'earliest collection date', 'count', 'mutation profile score', 'number score', 'trend score', 'overall score']
    cols = prefix_cols + [col for col in df_variants if col not in prefix_cols]
    df_variants = df_variants[cols]

    ''' Output variant table '''
    print('Variants:')
    print(df_variants)
    df_variants.to_csv(os.path.join(OUTDIR, 'variants_' + snapshot_date + '.csv'), sep=',')
    print('DONE')
    
