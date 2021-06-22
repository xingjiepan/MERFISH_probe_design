#!/usr/bin/env python3
'''Functions for constructing and manipulating the probe dictionary.
'''

import pandas as pd
from Bio.Seq import reverse_complement


def init_probe_dict(target_gene_ids:list, transcriptome:pd.core.frame.DataFrame, 
                        gene_id_key:str, K:int):
    '''Initialize the probe dictionary.
    Arguments:
        target_gene_ids: The ids of target genes.
        transcriptome: The transcriptome data frame.
        gene_id_key: The column names for the gene ids. Options are gene_short_name and gene_id.
        K: The size of each target region.
    
    Return:
        probe_dict: A dictionary of dictionary of data frames. The probes of a
                        transcript are stored in the data frame probe_dict[gene_id][transcript_id].
    ''' 
    # Get a sub-transcriptome of the target genes
    sub_transcriptome = transcriptome[transcriptome[gene_id_key].isin(target_gene_ids)]
    print(f'Found {sub_transcriptome.shape[0]} transcripts for {len(target_gene_ids)} target genes.')
    
    # Split all the target genes into overlapping K-mers
    probe_dict = {}
    for index, row in sub_transcriptome.iterrows():
        seq = row['sequence']
       
        gene_id = row[gene_id_key]
        transcript_id = row['transcript_id']
        
        # Add the entry to the dictionary if necessary
        if gene_id not in probe_dict:
            probe_dict[gene_id] = {}
        
        # Add the probes
        shifts = []
        sequences = []
        for shift in range(len(seq) - K + 1):
            shifts.append(shift)
            sequences.append(seq[shift:shift + K])
            
        probe_dict[gene_id][transcript_id] = pd.DataFrame({'gene_id':[gene_id] * len(shifts), 
            'transcript_id':[transcript_id] * len(shifts), 'shift':shifts, 'target_sequence':sequences})

    return probe_dict

def print_probe_dict(probe_dict:dict):
    '''Print the number of probes for each genes and transcripts.'''
    print('Gene\tTranscript\tN_probes')
    for gk in probe_dict.keys():
        print(gk)
        
        for tk in probe_dict[gk].keys():
            print(f'\t{tk}\t{probe_dict[gk][tk].shape[0]}')

def probe_dict_to_df(probe_dict:dict):
    '''Convert the probe dictionary into a single data frame.'''
    df = None
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            if df is None:
                df = probe_dict[gk][tk].copy()
            else:
                df = df.append(probe_dict[gk][tk], ignore_index=True)
    return df

def select_transcripts_by_ids(probe_dict:dict, transcript_ids:set):
    '''Select transcripts with specified ids.
    Return a new dictionary of probes.
    '''
    new_probe_dict = {}

    for gk in probe_dict.keys():
        new_probe_dict[gk] = {}
        
        for tk in probe_dict[gk].keys():
            if tk in transcript_ids:
                new_probe_dict[gk][tk] = probe_dict[gk][tk].copy()

    return new_probe_dict

def select_transcripts_by_num_probes(probe_dict:dict):
    '''For each gene, select the transcript that has most probes.
    Return a new dictionary of probes.
    '''
    new_probe_dict = {}

    for gk in probe_dict.keys():
        new_probe_dict[gk] = {}
        
        counts = [(tk, probe_dict[gk][tk].shape[0]) for tk in probe_dict[gk].keys()]
        tk_max = max(counts, key=lambda x:x[1])[0]

        new_probe_dict[gk][tk_max] = probe_dict[gk][tk_max].copy()

    return new_probe_dict

def get_rc_sequences(probe_dict:dict, input_column:str, output_column:str):
    '''Get the reverse complementary sequences of a column of sequences.'''
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            rc_seqs = [reverse_complement(seq) for seq in probe_dict[gk][tk][input_column]]
            probe_dict[gk][tk][output_column] = pd.Series(rc_seqs, index=probe_dict[gk][tk].index)


