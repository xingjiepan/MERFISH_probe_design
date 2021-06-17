#!/usr/bin/env python3

import pandas as pd
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp


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
            
            print(f'\t{tk}: {new_df.shape[0]} / {probe_dict[gk][tk].shape[0]} probes passed GC filter.')
            probe_dict[gk][tk] = new_df

def calc_tm_for_probe_dict(probe_dict:pd.core.frame.DataFrame, Na_conc:float, fmd_percentile:float):
    '''Calculate melting temperatures of the target sequences of the probe dictionary.
    Arguments:
        Na_conc: concentration of the the Na+ ion in mM.
        fmd_percentile: the percentile of formamide.
    '''
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            
            tms = []
            for seq in probe_dict[gk][tk]['target_sequence']:
                tm_raw = MeltingTemp.Tm_NN(seq, nn_table=MeltingTemp.DNA_NN4, Na=Na_conc)
                tms.append(MeltingTemp.chem_correction(tm_raw, fmd=fmd_percentile))

            probe_dict[gk][tk]['target_Tm'] = pd.Series(tms, index=probe_dict[gk][tk].index)

def filter_tm(probe_dict:pd.core.frame.DataFrame, lower_bound:float, upper_bound:float):
    '''Filter the melting temperatures of the target sequences.'''
    for gk in probe_dict.keys():
        print(gk)
        for tk in probe_dict[gk].keys():
            new_df = probe_dict[gk][tk][probe_dict[gk][tk]['target_Tm'].between(
                lower_bound, upper_bound, inclusive=True)]
            
            print(f'\t{tk}: {new_df.shape[0]} / {probe_dict[gk][tk].shape[0]} probes passed Tm filter.')
            probe_dict[gk][tk] = new_df

