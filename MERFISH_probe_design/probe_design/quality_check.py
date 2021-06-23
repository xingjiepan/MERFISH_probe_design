#!/usr/bin/env python3


import numpy as np
import pandas as pd


def check_and_standardize_transcriptome(transcriptome:pd.core.frame.DataFrame,
        remove_non_standard_columns:bool=False):
    '''Check the quality of the transcriptome and standardize it.
    Return a standardized transcriptome.'''
    standard_transcriptome = transcriptome

    # Check the existance of standarad columns
    standard_columns = ['transcript_id', 'sequence', 'gene_id', 'gene_short_name', 'FPKM']
    for sc in standard_columns:
        if not sc in transcriptome.columns:
            print(f'\033[91m ERROR: missing the standard column {sc}!')

    # Remove non-standard columns 
    if remove_non_standard_columns:
        nscs = [c for c in transcriptome.columns if c not in standard_columns]
        standard_transcriptome = transcriptome.drop(columns=nscs)

    # Check that the transcript ids are unique 
    t_ids, counts = np.unique(np.array(standard_transcriptome['transcript_id']), return_counts=True)
    for i in range(len(t_ids)):
        if counts[i] > 1:
            print(f'\033[91m ERROR: the transcript {t_ids[i]} have {counts[i]} entries!')
    
    return standard_transcriptome

def generate_transcript_level_report(probe_dict:dict, transcriptome:pd.core.frame.DataFrame):
    '''Generate a data frame of transcript level metrics.
    '''
    metrics = {'gene_id':[], 'gene_short_name':[], 'transcript_id':[], 'FPKM':[],
            'length':[], 'N_probes':[]}

    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            transcript_df = transcriptome[transcriptome['transcript_id'] == tk]
            
            metrics['gene_id'].append(transcript_df.iloc[0]['gene_id'])
            metrics['gene_short_name'].append(transcript_df.iloc[0]['gene_short_name'])
            metrics['transcript_id'].append(transcript_df.iloc[0]['transcript_id'])
            metrics['FPKM'].append(transcript_df.iloc[0]['FPKM'])
            metrics['length'].append(len(transcript_df.iloc[0]['sequence']))
            metrics['N_probes'].append(probe_dict[gk][tk].shape[0])

    return pd.DataFrame(metrics)
