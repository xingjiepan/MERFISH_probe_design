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

def barcode_str_to_array(bc_str:str):
    return np.array([int(c) for c in bc_str])

def barcode_array_to_str(bc_array:np.ndarray):
    return ''.join(['1' if i > 0 else '0' for i in bc_array])

def coverage_string(probe_bit_counts:np.ndarray):
    return ':'.join([str(i) for i in probe_bit_counts if i > 0])

def max_N_non_overlapping_probes(shifts:list, target_length:int):
    ss = sorted(shifts)
    N_accepted = 0
    tail = -1

    for s in ss:
        if s > tail:
            N_accepted += 1
            tail = s + target_length - 1

    return N_accepted

def generate_transcript_level_report(probe_dict:dict, transcriptome:pd.core.frame.DataFrame):
    '''Generate a data frame of transcript level metrics.
    '''
    metrics = {'gene_id':[], 'gene_short_name':[], 'transcript_id':[], 'FPKM':[],
            'length':[], 'barcode':[], 'N_probes':[], 'probe_bit_coverage':[],
            'max_N_non_overlapping_probes':[]}

    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            transcript_df = transcriptome[transcriptome['transcript_id'] == tk]

            # Add the basic metrics
            metrics['gene_id'].append(transcript_df.iloc[0]['gene_id'])
            metrics['gene_short_name'].append(transcript_df.iloc[0]['gene_short_name'])
            metrics['transcript_id'].append(transcript_df.iloc[0]['transcript_id'])
            metrics['FPKM'].append(transcript_df.iloc[0]['FPKM'])
            metrics['length'].append(len(transcript_df.iloc[0]['sequence']))
            metrics['N_probes'].append(probe_dict[gk][tk].shape[0])

            # Calculate barcode related metrics
            probe_barcodes = [barcode_str_to_array(bc_str) for bc_str in probe_dict[gk][tk]['probe_barcode']]
            probe_bit_counts = np.sum(probe_barcodes, axis=0)
            metrics['barcode'].append(barcode_array_to_str(probe_bit_counts))
            metrics['probe_bit_coverage'].append(coverage_string(probe_bit_counts))

            # Calculate the overlapping matric
            target_length = len(probe_dict[gk][tk].iloc[0]['target_sequence']) 
            metrics['max_N_non_overlapping_probes'].append(max_N_non_overlapping_probes(
                list(probe_dict[gk][tk]['shift']), target_length))

    return pd.DataFrame(metrics)
