#!/usr/bin/env python3
'''Python dictionary based implementation of the off-target table.
'''

import re
import pickle
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


def get_OTTable_for_sequences(sequences:list, K:int, weights:list=[]):
    '''Get an OTTable for a list of sequences.
    Arguments:
        sequences: A list of sequences.
        K: The size of K-mers to extract from the sequences.
        weights: The weights of sequences.
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

    return table

def get_OTTable_for_rtRNAs(ncRNAs:pd.core.frame.DataFrame, K:int):
    '''Get an OTTable for the rRNAs and tRNAs given a data frame of non-coding RNAs.
    Arguments:
        ncRNAs: A data frame of non-coding RNAs.
        K: THe size of K-mers for the OTTable.
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
