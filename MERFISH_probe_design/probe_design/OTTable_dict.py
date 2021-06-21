#!/usr/bin/env python3
'''Python dictionary based implementation of the off-target table.
'''

import re
import pickle
import numpy as np
import pandas as pd


class OTTable (dict):
    '''A python dictionary based implementation of the off-target table.'''
    def __missing__ (self, key):
        return 0

    def add_seq(self, seq, weight):
        '''Add a sequence to the OTTable with a given weight.'''
        self[seq] = self[seq] + weight

    def save_pkl(self, file_name):
        '''Save the OTTable as a pickle file.'''
        with open(file_name, 'wb') as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
        print(f'Wrote the OTTable to {file_name}.')
    
    @staticmethod
    def load_pkl(file_name):
        '''Load a OTTable from a pickle file.'''
        with open(file_name, 'rb') as f:
            ottable = pickle.load(f)
        print(f'Load the OTTable from {file_name}.')
        return ottable


def get_OTTable_for_sequences(sequences:list, K:int, weights:list=[], verbose:bool=False):
    '''Get an OTTable for a list of sequences.
    Arguments:
        sequences: A list of sequences.
        K: The size of K-mers to extract from the sequences.
        weights: The weights of sequences.
        verbose: Print the status for every 1,000 sequences. 
    '''
    # Use uniform weights if no weights are given
    if len(weights) == 0:
        weights = [1] * len(sequences)

    assert(len(weights) == len(sequences))

    table = OTTable()

    for i in range(len(sequences)):
        seq = sequences[i]
        w = weights[i]

        # Find all K-mers and add to the OTTable
        for j in range(len(seq) - K + 1): 
            table.add_seq(seq[j:j+K], w)

        if verbose and (i + 1) % 1000 == 0:
            print(f'Processed {i + 1}/{len(sequences)} sequences.')

    return table

def get_OTTable_for_rtRNAs(ncRNAs:pd.core.frame.DataFrame, K:int):
    '''Get an OTTable for the rRNAs and tRNAs given a data frame of non-coding RNAs.
    Arguments:
        ncRNAs: A data frame of non-coding RNAs.
        K: The size of K-mers for the OTTable.
    '''
    # Extract all rRNA and tRNA sequences
    biotypes_to_keep = ['rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA']

    rt_sequences = []
    for index, row in ncRNAs.iterrows():
        match = re.search(r'gene_biotype:(\S+) ', row['description'])
        if match is not None and match.group(1) in biotypes_to_keep:
            rt_sequences.append(row['sequence'])
    
    print(f'Found {len(rt_sequences)} rRNAs/tRNAs from {ncRNAs.shape[0]} non-coding RNAs.')

    # Create a OTTable using the sequences of rRNAs/tRNAs

    return get_OTTable_for_sequences(rt_sequences, K)

def get_OTTable_for_transcriptome(transcriptome:pd.core.frame.DataFrame, K:int, FPKM_threshold:float=0):
    '''Get an OTTable for the transcriptome.
    Arguments:
        transcriptome: A data frame of a transcriptome.
        K: The size of K-mers for the OTTable.
        FPKM_threshold: A transcript is included only if its FPKM is greater than this threshold.
    '''
    transcriptome_kept = transcriptome[transcriptome['FPKM'] > FPKM_threshold]
    print(f'Construct a OTTable using {transcriptome_kept.shape[0]}/{transcriptome.shape[0]} transcripts with FPKM > {FPKM_threshold}.')
   
    return get_OTTable_for_sequences(list(transcriptome_kept['sequence']), K, 
            weights=list(transcriptome_kept['FPKM']), verbose=True)

def get_gene_OTTables(transcriptome:pd.core.frame.DataFrame, target_gene_ids:list, gene_id_key:str, 
        K:int, FPKM_threshold:float=0):
    '''Get a dictionary of gene-level OTTables.
    Arguments:
        transcriptome: A data frame of a transcriptome.
        target_gene_ids: The ids of target genes.
        gene_id_key: The column names for the gene ids. Options are gene_short_name and gene_id.
        K: The size of K-mers for the OTTable.
        FPKM_threshold: A transcript is included only if its FPKM is greater than this threshold.
    '''
    gene_ottable_dict = {}
    
    for gene_id in target_gene_ids:
        print(f'Generate OTTable for gene {gene_id}.')
        gene_transcriptome = transcriptome[transcriptome[gene_id_key] == gene_id]
        gene_ottable_dict[gene_id] = get_OTTable_for_transcriptome(gene_transcriptome, K, FPKM_threshold=FPKM_threshold)

    return gene_ottable_dict
        
def calc_OTs(probe_dict:dict, ottable:OTTable, seq_key:str, ot_key:str, K:int):
    '''Calculate off-targets for sequences.
    Arguments:
        probe_dict: The probe dictionary.
        ottable: The OTTable for off-target calculation.
        seq_key: The column key for the sequences.
        ot_key: The column key to save the off-targets.
        K: The size of K-mers for the OTTable.
    '''
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            ot_counts = []

            for seq in probe_dict[gk][tk][seq_key]:
                ot_count = 0

                for i in range(len(seq) - K + 1):
                    ot_count += ottable[seq[i:i+K]]

                ot_counts.append(ot_count)

            probe_dict[gk][tk][ot_key] = pd.Series(ot_counts, index=probe_dict[gk][tk].index)
                
def calc_specificity(probe_dict:dict, ottable:OTTable, gene_ottable_dict:dict, transcript_fpkms:dict, 
        seq_key:str, speci_key:str, isospeci_key:str, K:int):
    '''Calculate specificities and isoform specificities for sequences.
    Definition:
        specificity of K-mer: the FPKM-weighted counts of the K-mer in all transcripts of a gene
                              divided by the FKPM-weighted counts of the K-mer in all transcripts.
        isospecificity of K-mer: the FPKM of the transcript divided by the FPKM-weighted counts of
                                the K-mer in all transcripts of a gene.
        specificity of sequence: the average of specificities of K-mers of the sequence.
        isospecificity of sequence: the average of isospecificities of K-mers of the sequence.

    Arguments:
        probe_dict: The probe dictionary.
        ottable: The OTTable for specificity calculation.
        gene_ottable_dict: A dictionary of OTTables for all genes of interest.
        transcript_fpkms: A dictionary that maps the names of transcripts to their FPKMs.
        seq_key: The column key for the sequences.
        speci_key: The column key to save the specificities.
        isospeci_key: The column key to save the isospecificities.
        K: The size of K-mers for the OTTable.
    '''
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            # Set the specificities to be zeros if the transcript do not express
            if 0 == transcript_fpkms[tk]:
                specificities = [0] * probe_dict[gk][tk].shape[0]
                isospecificities = [0] * probe_dict[gk][tk].shape[0]

            # Calculate specificities for probes that target expressed transcripts
            else:
                specificities = []
                isospecificities = []
                for seq in probe_dict[gk][tk][seq_key]:
                    
                    # Calculate the specificities for each K-mer
                    speci_K_mer = []
                    isospeci_K_mer = []
                    for i in range(len(seq) - K + 1):
                        speci_K_mer.append(gene_ottable_dict[gk][seq[i:i+K]] / ottable[seq[i:i+K]])
                        isospeci_K_mer.append(transcript_fpkms[tk] / gene_ottable_dict[gk][seq[i:i+K]])

                    # The specificity of the sequence is the average of its K-mer specificities
                    specificities.append(np.mean(speci_K_mer))
                    isospecificities.append(np.mean(isospeci_K_mer))

            # Update the data frame
            probe_dict[gk][tk][speci_key] = pd.Series(specificities, index=probe_dict[gk][tk].index)
            probe_dict[gk][tk][isospeci_key] = pd.Series(isospecificities, index=probe_dict[gk][tk].index)



