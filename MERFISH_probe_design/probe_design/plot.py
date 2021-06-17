#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt


def get_values_from_probe_dict(probe_dict:pd.core.frame.DataFrame, column_key:str):
    '''Get values under the column_key from a probe dictionary.
    Return a list of values.
    '''
    values = []
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            values += list(probe_dict[gk][tk][column_key])
    
    return values

def plot_gc_content(probe_dict:pd.core.frame.DataFrame):
    '''Plot the GC content distribution for all probes.'''
    plt.hist(get_values_from_probe_dict(probe_dict, 'target_GC'))
    plt.xlabel('GC content (%)')
    plt.ylabel('Count')
    plt.show()
