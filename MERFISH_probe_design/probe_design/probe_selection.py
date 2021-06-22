#!/usr/bin/env python3


from multiprocessing import Pool
import numpy as np
import pandas as pd


def select_probes_greedy_stochastic_one_df(df:pd.core.frame.DataFrame, N_probes_per_transcript:int, N_on_bits:int):
    '''A greedy stochastic method to select probes from one data frame.'''
    if N_probes_per_transcript >= df.shape[0]:
        print(f'There are only {df.shape[0]} probes while {N_probes_per_transcript} are required! Just return everything!')
        return df

    # Get the on-bits of the transcript
    on_bits = set()
    for bc in df['probe_barcode']:
        for i in range(len(bc)):
            if bc[i] == '1':
                on_bits.add(i)
    on_bits = list(on_bits)

    if len(on_bits) != N_on_bits:
        raise Exception(f'Probes for {df.iloc[0]["gene_id"]}:{df.iloc[0]["transcript_id"]} have {len(on_bits)} on-bits instead of {N_on_bits}!')

    # Create a array to track the coverage of the transcript
    target_length = len(df.iloc[0]['target_sequence'])
    max_targetable_length = np.max(df['shift']) + target_length
    transcript_coverage = np.zeros(max_targetable_length)

    # Create a dictionary to track the coverage of the on-bits 
    on_bit_coverage = {ob:0 for ob in on_bits}

    # Select probes
    selected_indices = []
    rest_indices = [i for i in range(df.shape[0])]
    for i in range(N_probes_per_transcript):
        
        # Calculate the scores if a probe is to be added
        #trial_scores = [score_probe_subset(df, selected_indices + [r_id], on_bits) for r_id in rest_indices]
        trial_scores = []
        for r_id in rest_indices:
            # Calculate the number of overlaps
            shift = df.iloc[r_id]['shift']
            N_new_overlaps = np.sum([transcript_coverage[pos] for pos in range(shift, shift + target_length)])

            # Calculate the number of bit coverage
            trial_on_bit_coverage = on_bit_coverage.copy()
            bc = df.iloc[r_id]['probe_barcode']
            for i in range(len(bc)):
                if bc[i] == '1':
                    trial_on_bit_coverage[i] += 1
    
            # Define a score to penalize overlaps and unevenly distributed probes
            score = N_new_overlaps + max(trial_on_bit_coverage.values()) / (1 + min(trial_on_bit_coverage.values()))
            trial_scores.append(score)

        # Get the indices with the lowest score
        score_min = min(trial_scores)
        lowest_score_ids = [rest_indices[j] for j in range(len(rest_indices)) if trial_scores[j] == score_min]

        # Randomly select an ID
        selected_id = np.random.choice(lowest_score_ids)
        selected_indices.append(selected_id)
        rest_indices.remove(selected_id)

        # Update the tracking records
        # Calculate the number of overlaps
        shift = df.iloc[selected_id]['shift']
        transcript_coverage[shift:shift + target_length] += 1
        
        # Calculate the number of bit coverage
        bc = df.iloc[selected_id]['probe_barcode']
        for i in range(len(bc)):
            if bc[i] == '1':
                on_bit_coverage[i] += 1

    print(f'{df.iloc[0]["gene_id"]}:{df.iloc[0]["transcript_id"]}: selected {N_probes_per_transcript}/{df.shape[0]} probes with N_overlapping_bases={np.sum(transcript_coverage * (transcript_coverage - 1)  / 2)} and on-bit_coverage={on_bit_coverage}.')

    # Return a data frame with the selected indices
    return df.iloc[selected_indices]


def select_probes_greedy_stochastic(probe_dict:dict, N_probes_per_transcript:int, N_on_bits:int=4, N_threads:int=1):
    '''A greedy stochastic method to select probes.
    Arguments:
        probe_dict: The dictionary of probes.
        N_on_bits: The number of on bits each probe should have.
    '''
    keys = []
    args = []
    for gk in probe_dict.keys(): 
        for tk in probe_dict[gk].keys():
            keys.append((gk, tk))
            args.append([probe_dict[gk][tk], N_probes_per_transcript, N_on_bits]) 
    
    with Pool(N_threads) as p:
        results = p.starmap(select_probes_greedy_stochastic_one_df, args)

    for i in range(len(keys)):
        gk, tk = keys[i]
        probe_dict[gk][tk] = results[i]
