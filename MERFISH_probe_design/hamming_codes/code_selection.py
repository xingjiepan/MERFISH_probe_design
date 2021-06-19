#!/usr/bin/python3


import numpy as np

def find_subset_code_even_coverage_single_round(code_list:list, code_length:int, N_codes:int):
    '''Find a subset of N_codes code that evenly cover each bits
    by a greedy method. The method is stochastic.
    Return:
        selected_codes. rest_codes, variance of the coverage
    '''
    assert(len(code_list) >= N_codes)

    bits_coverage = np.zeros(code_length)

    rest_codes = code_list
    selected_codes = [] 

    # Choose codes iteratively
    for i in range(N_codes):
        
        # Test the variances of the bit-coverage for each possible code-addition
        trial_vars = []
        for trial_code in rest_codes:
            
            bits_coverage_trial = bits_coverage.copy()
            for b in trial_code:
                bits_coverage_trial[b] += 1

            trial_vars.append(np.var(bits_coverage_trial))

        # Find the codes that gives the lowest variance
        min_var = np.min(trial_vars)
        candidate_ids = [j for j in range(len(rest_codes)) if trial_vars[j] == min_var]
        selected_id = np.random.choice(candidate_ids)

        # Update the code lists
        rest_codes_new = []
        for j in range(len(rest_codes)):
            if j == selected_id:
                selected_codes.append(rest_codes[j])

                for b in rest_codes[j]:
                    bits_coverage[b] += 1

            else:
                rest_codes_new.append(rest_codes[j])

        rest_codes = rest_codes_new

    return selected_codes, rest_codes, np.var(bits_coverage)


def find_subset_code_even_coverage(code_list:list, code_length:int, N_codes:int, N_rounds:int=10):
    '''Find a subset of N_codes code that evenly cover each bits
    by a greedy method. The method is stochastic. Run the underlying
    method N_rounds times and select the one that gives the least variance.
    '''
    best_var = np.Inf
    
    for i in range(N_rounds):
        result = find_subset_code_even_coverage_single_round(code_list, code_length, N_codes)
       
        print(f'Round = {i + 1}, variance = {result[2]}')
        if best_var > result[2]:
            selected_codes = result[0]
            rest_codes = result[1]
            best_var = result[2]

    return selected_codes, rest_codes
