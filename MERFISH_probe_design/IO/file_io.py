#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO


def load_fasta_into_df(fasta_file:str):
    '''Load a fasta file into a pandas data frame.'''
    d = {'id':[], 'description':[], 'sequence':[]}
    
    for record in SeqIO.parse(fasta_file, 'fasta'):
        d['id'].append(record.id)
        d['description'].append(record.description)
        d['sequence'].append(str(record.seq))

    return pd.DataFrame.from_dict(d)

def load_merlin_codebook(codebook_file:str):
    '''Load the MERlin style codebook.'''
    version = ''
    codebook_name = ''
    bit_names = []
    barcode_dict = {'name':[], 'id':[], 'barcode_str':[]}
    
    with open(codebook_file, 'r') as f:
        lines = f.readlines()
        is_header = True
    
        for l in lines:
            sl = l.split(',')
            
            # Skip blank lines
            if len(sl) == 0:
                continue
            
            # Load the header
            if is_header:
                if sl[0].strip() == 'version':
                    version = sl[1].strip()
                elif sl[0].strip() == 'codebook_name':
                    codebook_name = sl[1].strip()
                elif sl[0].strip() == 'bit_names':
                    bit_names = [sl[i].strip() for i in range(1, len(sl))]
                elif sl[0].strip() == 'name':
                    is_header = False
                    continue
                
            # Load the barcode table
            else:
                barcode_dict['name'].append(sl[0].strip())
                barcode_dict['id'].append(sl[1].strip())
                barcode_dict['barcode_str'].append(sl[2].strip())
    
    return version, codebook_name, bit_names, pd.DataFrame.from_dict(barcode_dict)
