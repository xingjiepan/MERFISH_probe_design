#!/usr/bin/env python3

import numpy as np
import pandas as pd


def randomly_select_primers_with_lowest_OT(primers:pd.core.frame.DataFrame):
    total_OTs = np.array(primers['sequence_OT']) + np.array(primers['sequence_rc_OT'])
    min_OT = min(total_OTs)

    candidate_ids = [i for i in range(len(total_OTs)) if total_OTs[i] == min_OT]
    selected_id = np.random.choice(candidate_ids)

    return primers.iloc[[selected_id]]

def add_primer_sequences(probe_dict:dict, seq_upstream:str, seq_downstream:str, 
        input_column:str='target_readout_sequence', output_column:str='target_readout_primer_sequence'):
    '''Add the primer sequences to a probe dictionary.'''
    for gk in probe_dict.keys():
        for tk in probe_dict[gk].keys():
            output_seqs = []

            for seq in probe_dict[gk][tk][input_column]:
                output_seqs.append(seq_upstream + seq + seq_downstream)

            probe_dict[gk][tk][output_column] = pd.Series(output_seqs, index=probe_dict[gk][tk].index)
