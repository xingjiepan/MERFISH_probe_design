#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO


def load_fasta_into_df(fasta_file):
    '''Load a fasta file into a pandas data frame.'''
    d = {'id':[], 'description':[], 'sequence':[]}
    
    for record in SeqIO.parse(fasta_file, 'fasta'):
        d['id'].append(record.id)
        d['description'].append(record.description)
        d['sequence'].append(str(record.seq))

    return pd.DataFrame.from_dict(d)

