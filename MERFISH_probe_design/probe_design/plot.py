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

def plot_hist(probe_dict:pd.core.frame.DataFrame, column_key:str, y_max=None):
    '''Plot histogram of values under the column_key'''
    plt.hist(get_values_from_probe_dict(probe_dict, column_key))

    if y_max is not None:
        bottom, top = plt.ylim()
        plt.ylim(bottom, y_max)

    plt.xlabel(column_key)
    plt.ylabel('Count')
    plt.show()

def plot_gc_content(probe_dict:pd.core.frame.DataFrame):
    '''Plot the GC content distribution for all probes.'''
    plt.hist(get_values_from_probe_dict(probe_dict, 'target_GC'))
    plt.xlabel('GC content (%)')
    plt.ylabel('Count')
    plt.show()

def plot_tm(probe_dict:pd.core.frame.DataFrame):
    '''Plot the Tm distribution for all probes.'''
    plt.hist(get_values_from_probe_dict(probe_dict, 'target_Tm'))
    plt.xlabel('Tm (°C)')
    plt.ylabel('Count')
    plt.show()

