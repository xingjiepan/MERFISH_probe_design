#!/usr/bin/env python3
'''Functions for analyzing the qualities
of Hamming codes.
'''

import numpy as np
import matplotlib.pyplot as plt


def H_dist(code1:set, code2:set):
    return len(code1.union(code2)) - len(code1.intersection(code2))  


def plot_pairwise_hamming_distance_distribution(code_list:list):
    # Make sure that the on-bits are in sets for H_dist calculation
    code_list = [set(c) for c in code_list]

    pairwise_dists = []

    for i in range(len(code_list)):
        for j in range(i + 1, len(code_list)):
            pairwise_dists.append(H_dist(code_list[i], code_list[j]))

    plt.hist(pairwise_dists)
    plt.xlabel('Hamming distance')
    plt.ylabel('Count')
    plt.show()

def plot_bit_coverage(code_list:list, code_length):
    coverages = np.zeros(code_length)

    for code in code_list:
        for bit in code:
            coverages[bit] += 1

    plt.scatter(np.arange(code_length), coverages)
    plt.xlabel('Bit')
    plt.ylabel('Count')
    plt.show()
