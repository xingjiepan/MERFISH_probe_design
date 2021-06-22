#!/usr/bin/env python3

import numpy as np
import pandas as pd
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp


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

def calc_tm_JM(sequence:str, monovalentSalt:float=0.3, probeConc:float=5e-9):
    '''TM calculation used in JM's original MATLAB code.
    Adapted by Rongxin Fang.
    '''
    intSeq = np.array([["A", "C", "G", "T"].index(s) for s in list(sequence)])
    nnID = intSeq[:(len(intSeq)-1)] * 4 + intSeq[1:]
    end = len(intSeq)

    isValid = np.array([
        False not in (a, b, c, d) for (a, b, c, d) in zip(intSeq[:(end-1)] <= 3,
        intSeq[1:end] <= 3,
        intSeq[:(end-1)] >= 0,
        intSeq[1:end] >=0) ])

    H = np.array([
        -7.6, -8.4, -7.8, -7.2, -8.5, -8.0, -10.6, -7.8,
        -8.2, -9.8, -8.0, -8.4, -7.2, -8.2, -8.5, -7.6])   
    S = np.array([
        -21.3, -22.4, -21.0, -20.4, -22.7, -19.9, -27.2, -21.0,
        -22.2, -24.4, -19.9, -22.4, -21.3, -22.2, -22.7, -21.3])

    dG = np.zeros([2, len(intSeq)-1])
    dG[0,:] = H[nnID[isValid]]
    dG[1,:] = S[nnID[isValid]]

    H = np.cumsum(dG[0,:])[-1]
    S = np.cumsum(dG[1,:])[-1]

    # Determine ends
    fivePrimeAT = (intSeq[0] == 0) | (intSeq[0]  == 3);
    threePrimeAT = (intSeq[-1] == 0) | (intSeq[-1] == 3);

    H = H + 0.2 + 2.2*fivePrimeAT + 2.2*threePrimeAT;
    S = S + -5.7 + 6.9*fivePrimeAT + 6.9*threePrimeAT;

    S = S + 0.368*(len(sequence)-1)*np.log(monovalentSalt);

    return H*1000 / (S + 1.9872 * np.log(probeConc)) - 273.15;

def calc_tm_JM_for_probe_dict(probe_dict:pd.core.frame.DataFrame, monovalentSalt:float, probe_conc:float=1,
        column_key_seq:str='target_sequence', column_key_write='target_Tm'):
    '''Calculate melting temperatures of the target sequences of the probe dictionary.
    Use the TM calculation method in JM's original MATLAB code.
    Arguments:
        Na_conc: concentration of the the Na+ ion in mM.
        fmd_percentile: the percentile of formamide.
        probe_conc: concentration of the individual probes in nM.
    '''
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            
            tms = []
            for seq in probe_dict[gk][tk][column_key_seq]:
                tms.append(calc_tm_JM(seq, monovalentSalt, probe_conc))

            probe_dict[gk][tk][column_key_write] = pd.Series(tms, index=probe_dict[gk][tk].index)

