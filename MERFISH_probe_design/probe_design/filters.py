#!/usr/bin/env python3

from multiprocessing import Pool

import numpy as np
import pandas as pd

try:
    from Bio.SeqUtils import gc_fraction
    def GC(sequence):
        return 100 * gc_fraction(sequence, ambiguous="ignore")
except ImportError:
    # Older versions have this:
    from Bio.SeqUtils import GC

from Bio.SeqUtils import MeltingTemp

# GLOBAL VARIABLES
from . import seq2int

def filter_probe_dict_by_metric(probe_dict:pd.core.frame.DataFrame, column_key:str, 
        lower_bound:float=-np.inf, upper_bound:float=np.inf):
    '''Filter the probe dictionary by a metric.'''
    for gk in probe_dict.keys():
        print(gk)
        for tk in probe_dict[gk].keys():
            new_df= probe_dict[gk][tk][
                probe_dict[gk][tk][column_key].gt(lower_bound) & 
                probe_dict[gk][tk][column_key].lt(upper_bound)
            ]
            
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

def calc_tm_JM_for_transcript(df, monovalentSalt, probe_conc, column_key_seq, column_key_write):
    tms = []
    for seq in df[column_key_seq]:
        tms.append(calc_tm_JM(seq, monovalentSalt, probe_conc))

    df[column_key_write] = pd.Series(tms, index=df.index)
    return df

def calc_tm_JM_for_probe_dict(probe_dict:dict, monovalentSalt:float, probe_conc:float=1,
        column_key_seq:str='target_sequence', column_key_write='target_Tm', n_threads=1):
    '''Calculate melting temperatures of the target sequences of the probe dictionary.
    Use the TM calculation method in JM's original MATLAB code.
    Arguments:
        Na_conc: concentration of the the Na+ ion in mM.
        fmd_percentile: the percentile of formamide.
        probe_conc: concentration of the individual probes in nM.
    '''
    # Iterate through all genes and get the arguments for parallel processing
    ks = []
    args = []
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            ks.append((gk, tk))
            args.append([probe_dict[gk][tk], monovalentSalt, probe_conc, column_key_seq, column_key_write])
    # Add readout probes in parallel
    with Pool(n_threads) as p:
        results = p.starmap(calc_tm_JM_for_transcript, args)
    
    # Update the probe dictionary
    for i, kk in enumerate(ks):
        gk, tk = kk
        probe_dict[gk][tk] = results[i]

## Add functions related to Flex probe design:
# calculate GC for left 25 bases:

def calc_GC_for_left_probe_dict(
    probe_dict:pd.core.frame.DataFrame, 
    left_length:int=25,
    column_key_seq:str='target_sequence',
    column_key_write='target_GC_left'):
    """Calculate GC content of the left part of the sequences under the column_key_seq column in the probe dictionary.
    The GC content is reported in percentile.
    """
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            gcs = []
            for seq in probe_dict[gk][tk][column_key_seq]:
                gcs.append(GC(seq[:left_length]))

            probe_dict[gk][tk][column_key_write+str(left_length)] = pd.Series(gcs, index=probe_dict[gk][tk].index)
    
def calc_GC_for_right_probe_dict(
    probe_dict:pd.core.frame.DataFrame, 
    right_length:int=25,
    column_key_seq:str='target_sequence',
    column_key_write='target_GC_right'):
    """Calculate GC content of the right part of the sequences under the column_key_seq column in the probe dictionary.
    The GC content is reported in percentile.
    """
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            gcs = []
            for seq in probe_dict[gk][tk][column_key_seq]:
                gcs.append(GC(seq[-right_length:]))

            probe_dict[gk][tk][column_key_write+str(right_length)] = pd.Series(gcs, index=probe_dict[gk][tk].index)  

def calc_base_at_location_probe_dict(
    probe_dict:pd.core.frame.DataFrame, 
    location:int=25,
    column_key_seq:str='target_sequence',
    column_key_write='target_base_at_location'):
    """Calculate the base at a specific location in the sequences under the column_key_seq column in the probe dictionary.
    The base is reported as a string.
    """
    
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            bases = []
            for seq in probe_dict[gk][tk][column_key_seq]:
                bases.append(seq2int[seq[location].upper()])
            probe_dict[gk][tk][column_key_write+str(location)] = pd.Series(bases, index=probe_dict[gk][tk].index)

# function to remove probes with homopolymers:
def find_homopolymer(seq, max_length=4):
    """Check if the sequence has a homopolymer longer than max_length."""
    for base in 'ATCG':
        if base * max_length in seq.upper():
            return seq2int[base]  # Return the base if a homopolymer is found
    return 0
def detect_homopolymer_probe_dict(
    probe_dict:pd.core.frame.DataFrame, 
    column_key_seq:str='target_sequence',    
    min_length:int=5,
    column_key_write='has_homopolymer',
    ):
    """Remove probes with homopolymers longer than min_length."""
    
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            seqs = probe_dict[gk][tk][column_key_seq]
            probe_dict[gk][tk][column_key_write+str(min_length)] = seqs.apply(
                lambda seq: find_homopolymer(seq, min_length))