#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_values_from_probe_dict(probe_dict:dict, column_key:str):
    '''Get values under the column_key from a probe dictionary.
    Return a list of values.
    '''
    values = []
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            values += list(probe_dict[gk][tk][column_key])
    
    return values

def plot_hist(probe_dict:dict, column_key:str, y_max=None, bins=30):
    '''Plot histogram of values under the column_key'''
    plt.hist(get_values_from_probe_dict(probe_dict, column_key), bins=bins)

    if y_max is not None:
        bottom, top = plt.ylim()
        plt.ylim(bottom, y_max)

    plt.xlabel(column_key)
    plt.ylabel('Count')
    plt.show()

def plot_sequence_coverage(df:pd.core.frame.DataFrame, seq_length:int):
    '''Plot the sequence coverage of a sequences.'''
    coverage = np.zeros(seq_length, dtype=int)

    shifts = list(df['shift'])
    target_seqs = list(df['target_sequence'])

    for i in range(len(shifts)):
        shift = shifts[i]
        t_len = len(target_seqs[i])
        coverage[shift : shift + t_len] += 1

    plt.plot(np.arange(seq_length), coverage)
    plt.xlabel('Sequence position')
    plt.ylabel('Coverage')
    plt.show()

