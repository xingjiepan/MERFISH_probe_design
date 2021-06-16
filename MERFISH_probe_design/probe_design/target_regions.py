#!/usr/bin/env python3
'''Functions for selecting target regions.
'''

import pandas as pd


def init_target_regions(target_gene_ids:list, transcriptome:pd.core.frame.DataFrame, 
                        gene_id_key:str, K:int):
    '''Initialize the target regions.
    Arguments:
        target_gene_ids: The ids of target genes.
        transcriptome: The transcriptome data frame.
        get_id_key: The column names for the gene ids. Options are gene_short_name and gene_id.
        K: The size of each target region.
    
    Return:
        target_regions: A dictionary of dictionary of data frames. The target sequences of a
                        transcript are stored in the data frame target_regions[gene_id][transcript_id].
    ''' 
    # Get a sub-transcriptome of the target genes
    sub_transcriptome = transcriptome[transcriptome[gene_id_key].isin(target_gene_ids)]
    print(f'Found {sub_transcriptome.shape[0]} transcripts for {len(target_gene_ids)} target genes.')
    
    # Split all the target genes into overlapping K-mers
    target_regions = {}
    for index, row in sub_transcriptome.iterrows():
        seq = row['sequence']
       
        gene_id = row[gene_id_key]
        transcript_id = row['transcript_id']
        
        # Add the entry to the dictionary if necessary
        if gene_id not in target_regions:
            target_regions[gene_id] = {}
        
        # Add the target regions
        shifts = []
        sequences = []
        for shift in range(len(seq) - K + 1):
            shifts.append(shift)
            sequences.append(seq[shift:shift + K])
            
        target_regions[gene_id][transcript_id] = pd.DataFrame({'shift':shifts, 'sequence':sequences})

    return target_regions

def print_target_regions(target_regions:dict):
    '''Print the number target regions for each genes and transcripts.'''
    print('Gene\tTranscript\tN_target_regions')
    for gk in target_regions.keys():
        print(gk)
        
        for tk in target_regions[gk].keys():
            print(f'\t{tk}\t{target_regions[gk][tk].shape[0]}')
    

