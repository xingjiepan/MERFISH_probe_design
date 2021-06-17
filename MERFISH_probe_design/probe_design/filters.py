#!/usr/bin/env python3

import pandas as pd
from Bio.SeqUtils import GC


def calc_gc_for_probe_dict(probe_dict:pd.core.frame.DataFrame):
    '''Calculate GC content of the target sequences of the probe dictionary.
    The GC content is reported in percentile.
    '''
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            
            gcs = []
            for seq in probe_dict[gk][tk]['target_sequence']:
                gcs.append(GC(seq))

            probe_dict[gk][tk]['target_GC'] = pd.Series(gcs, index=probe_dict[gk][tk].index)

def filter_gc_content(probe_dict:pd.core.frame.DataFrame, lower_bound:float, upper_bound:float):
    '''Filter the GC content of the target sequences.'''
    for gk in probe_dict.keys():
        print(gk)
        for tk in probe_dict[gk].keys():
            new_df = probe_dict[gk][tk][probe_dict[gk][tk]['target_GC'].between(
                lower_bound, upper_bound, inclusive=True)]
            
            print(f'{tk}: {new_df.shape[0]} / {probe_dict[gk][tk].shape[0]} probes passed GC filter.')
            probe_dict[gk][tk] = new_df
