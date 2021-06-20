#!/usr/bin/env python3

import numpy as np
import pandas as pd
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp


def calc_gc_for_probe_dict(probe_dict:pd.core.frame.DataFrame, 
        column_key_seq:str='target_sequence', column_key_write='target_GC'):
    '''Calculate GC content of sequences under the column_key_seq column in the probe dictionary.
    The GC content is reported in percentile.
    '''
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            
            gcs = []
            for seq in probe_dict[gk][tk][column_key_seq]:
                gcs.append(GC(seq))

            probe_dict[gk][tk][column_key_write] = pd.Series(gcs, index=probe_dict[gk][tk].index)

def calc_tm_for_probe_dict(probe_dict:pd.core.frame.DataFrame, Na_conc:float, fmd_percentile:float, probe_conc:float=1,
        column_key_seq:str='target_sequence', column_key_write='target_Tm'):
    '''Calculate melting temperatures of the target sequences of the probe dictionary.
    Arguments:
        Na_conc: concentration of the the Na+ ion in mM.
        fmd_percentile: the percentile of formamide.
        probe_conc: concentration of the individual probes in nM.
    '''
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            
            tms = []
            for seq in probe_dict[gk][tk][column_key_seq]:
                tm_raw = MeltingTemp.Tm_NN(seq, nn_table=MeltingTemp.DNA_NN4, Na=Na_conc,
                        dnac1=probe_conc, dnac2=0)
                
                tms.append(MeltingTemp.chem_correction(tm_raw, fmd=fmd_percentile))

            probe_dict[gk][tk][column_key_write] = pd.Series(tms, index=probe_dict[gk][tk].index)

def filter_probe_dict_by_metric(probe_dict:pd.core.frame.DataFrame, column_key:str, 
        lower_bound:float=-np.Inf, upper_bound:float=np.Inf):
    '''Filter the probe dictionary by a metric.'''
    for gk in probe_dict.keys():
        print(gk)
        for tk in probe_dict[gk].keys():
            new_df = probe_dict[gk][tk][probe_dict[gk][tk][column_key].between(
                lower_bound, upper_bound, inclusive=False)]
            
            print(f'\t{tk}: {new_df.shape[0]} / {probe_dict[gk][tk].shape[0]} probes passed the filter {lower_bound} < {column_key} <  {upper_bound}.')
            probe_dict[gk][tk] = new_df


