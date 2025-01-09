#!/usr/bin/env python3
'''This script contains a method to design a MERFISH
codebook to evenly distribute the on-bits.
The method was originally written by Will Alen and then
modified by Xingjie Pan.
'''

from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt


def on_bits_to_binary_code(on_bits:list, code_length:int):
    '''Convert a list of on-bits to a list of binary code.'''
    return [1 if i in on_bits else 0 for i in range(code_length)]

def calc_dot_distribution(ct_expr:np.ndarray, binary_codes:np.ndarray):
    '''Calculate the dot distribution in cell types at different bits.'''
    return np.matmul(ct_expr, binary_codes)

def calc_assignment_score(ct_expr:np.ndarray, ct_weights:np.ndarray, binary_codes:np.ndarray):
    '''Calculate the assignment score.
    The score is the weight variance. The lower the better.
    '''
    dot_distribution = calc_dot_distribution(ct_expr, binary_codes)
    return np.sum(np.var(dot_distribution, axis=1) * ct_weights)

def random_swap_two_rows(a:np.ndarray):
    '''Randomly swap two rows of a numpy array.'''
    i, j = np.random.choice(range(a.shape[0]), 2, replace=False)
    a[[i, j]] = a[[j, i]]

def optimize_bit_assignments_simulated_annealing(ct_expr:np.ndarray, 
        ct_weights:np.ndarray, binary_codes:np.ndarray, N_rounds:int=5, N_iter:int=2000):
    ''' Optimize the bit assignments by simulated annealing.
    Arguments:
        ct_expr: N_cell_types x N_genes expression matrix. Normalized such that the sum of all elements is N_cell_types.
        ct_weights: N_cell_types x 1. Weights for cell types normalized to sum one.
        binary_code: Matrix of binary codes.
        N_rounds: Number of rounds.
        N_iter: Number of iterations per round.
    NOTE: ct_expr and ct_weights must be in same order
    '''
    assert(ct_expr.shape[0] == ct_weights.shape[0])
    assert(ct_expr.shape[1] == binary_codes.shape[0])
    
    # Re-seed the generator for each thread
    np.random.seed()

    # Initialize the curent state
    current_codes = binary_codes
    current_score = calc_assignment_score(ct_expr, ct_weights, current_codes)
    best_score = current_score
    best_codes = binary_codes
    init_temperature = current_score

    # Run the simulated annealing
    for i in range(N_rounds):
        tmp = init_temperature

        for j in range(N_iter):
            tmp = max(tmp * 0.95, 1e-30) # Avoid numeric errors for too small numbers

            # Do the trial
            trial_codes = current_codes.copy()
            random_swap_two_rows(trial_codes)
            trial_score = calc_assignment_score(ct_expr, ct_weights, trial_codes)
            
            # Determine if accept the trial
            score_diff = trial_score - current_score

            if score_diff < 0 or np.random.uniform() < np.exp( - score_diff / tmp):
                current_codes = trial_codes
                current_score = trial_score

                # Update the best assignment
                if current_score < best_score:
                    best_score = current_score
                    best_codes = current_codes

    print(f'Simulated annealing finished with best score = {best_score}.')
    return best_score, best_codes

def optimize_bit_assignments_simulated_annealing_parallel(ct_expr:np.ndarray, 
        ct_weights:np.ndarray, binary_codes:np.ndarray, N_test, N_threads, N_rounds:int=5, N_iter:int=2000):
    '''Run the simulated annealing in parallel.'''
    args = [[ct_expr, ct_weights, binary_codes, N_rounds, N_iter] for i in range(N_test)]

    # Run simulated annealing for N_test times 
    if N_threads == 1:
        print("Not multi-processing if N_threads==1. ")
        results = [optimize_bit_assignments_simulated_annealing(*_arg) for _arg in args]
    with Pool(N_threads) as p:
        results = p.starmap(optimize_bit_assignments_simulated_annealing, args)

    # Return the best result
    best_id = np.argmin([r[0] for r in results]) 
    return results[best_id]

def plot_dot_distribution(ct_expr:np.ndarray, binary_codes:np.ndarray):
    '''Plot the distribution of numbers of dots in all cell types at all bits.'''
    dot_distribution = calc_dot_distribution(ct_expr, binary_codes)

    fig, ax = plt.subplots()
    im = ax.imshow(np.transpose(dot_distribution), vmin=0)
    cbar = ax.figure.colorbar(im)

    plt.xlabel('Cell type')
    plt.ylabel('Bit')
    plt.show()
