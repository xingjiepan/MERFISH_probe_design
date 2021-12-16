#!/usr/bin/env python3
'''Generate modified Hamming code by randomly choosing codes.
Multiple people from the Zhuang lab contributed to this version
including Bogdan Bintu, Po Zheng, Will Allen and Xingjie Pan.
'''
import numpy as np
from itertools import combinations
from time import time
from multiprocessing import Pool


def H_dist(code1:set, code2:set):
    return len(code1.union(code2)) - len(code1.intersection(code2))  

def bit_coverage(code_list:list, code_length):
    coverages = np.zeros(code_length)

    for code in code_list:
        for bit in code:
            coverages[bit] += 1

    return coverages

def generate_one_code_set(n_bits, n_on_bits, min_hamming_distance, randomize, verbose):
    '''Randomly generate one set of codes.'''
    # Re-seed the generator for each thread
    np.random.seed()

    # Get all possible codes with n_bits and n_on_bits 
    candidate_codes = list(combinations(list(range(n_bits)), n_on_bits))
    if randomize:
        np.random.shuffle(candidate_codes)

    # Iterate through all codes
    chosen_codes = []
    while len(candidate_codes) > 0:
        c = set(candidate_codes.pop(0))
        
        # Keep the code if compatible
        keep = True
        for chosen_c in chosen_codes:
            if H_dist(c, chosen_c) < min_hamming_distance:
                keep = False
                break
        if keep:
            chosen_codes.append(c)
    
    coverages = bit_coverage(chosen_codes, n_bits)

    if verbose:
        print(f'Generated code set: len={len(chosen_codes)}, var={np.var(coverages)}')
    
    return chosen_codes


def generate_modified_hamming_codes(n_bits:int, n_on_bits:int, min_hamming_distance:int, 
        min_codebook_size:int=0, n_rand_repeats:int=500, n_threads:int=1, verbose=True):
    '''Generate code sets and pick the best one.
    The code sets include a regular one plus n_rand_repeats random sets.
    '''

    def code_set_score(size, var):
        '''Score the code set based on its size and variance.
        The higher the better.
        '''
        return size - 5 * var 

    start = time() 

    best_coding = []
    best_var = np.inf

    # Generate code sets in parallel
    n_repeats = 1 + n_rand_repeats

    args = [[n_bits, n_on_bits, min_hamming_distance, True, verbose] for i in range(n_repeats)]
    args[0][3] = False # Don't randomize the first set

    with Pool(n_threads) as p:
        all_code_sets = p.starmap(generate_one_code_set, args)

    # Find the best code set
    for i in range(n_repeats):

        chosen_codes = all_code_sets[i]

        # Count the coverage for each bit

        coverages = bit_coverage(chosen_codes, n_bits)
   
        # Update the best code set if the current set is better in both size and variance

        if len(chosen_codes) > min_codebook_size and \
                code_set_score(len(chosen_codes), np.var(coverages)) > code_set_score(len(best_coding), best_var):
            best_coding = chosen_codes
            best_var = np.var(coverages)

    print(f'The best code set has {len(best_coding)} codes. The coverage variance is {best_var}.')

    end = time()
    print("-- Duration: ", end-start)
    
    return best_coding


if __name__ == '__main__':
    code_list = generate_modified_hamming_codes(20, 4, 4, n_rand_repeats=50, n_threads=2)
    print(code_list)

