#!/usr/bin/env python3
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import reverse_complement
import re
ensembl_full_regexp = r'([a-zA-Z0-9\.]+) ([a-zA-Z]+) (chromosome|scaffold):([a-zA-Z0-9\.\:]+)\:(1|\-1) gene:([a-zA-Z0-9\.]+) gene_biotype:([a-zA-Z0-9\.\_]+) (.+)'
gencode_regexp = r'(?P<transcript_id>[A-Z0-9\.]+)\|(?P<gene_id>[A-Z0-9\.]+)\|(?P<havana_gene>[A-Z0-9\.-]+)\|(?P<havana_transcript>[A-Z0-9\.\-]+)\|(?P<transcript_name>[A-Za-z0-9\.\-\_\(\)]+)\|(?P<gene_name>[A-Za-z0-9\.\-\_\(\)]+)\|(?P<length>[0-9]+)\|(?P<transcript_type>[A-Za-z0-9\.\-\_]+)\|'

def load_fasta_into_df(fasta_file:str, load_rc:bool=False):
    '''Load a fasta file into a pandas data frame.'''
    d = {'id':[], 'description':[], 'sequence':[]}
    if load_rc:
        d['sequence_rc'] = []

    for record in SeqIO.parse(fasta_file, 'fasta'):
        d['id'].append(record.id)
        d['description'].append(record.description)
        d['sequence'].append(str(record.seq))

        if load_rc:
            d['sequence_rc'].append(reverse_complement(str(record.seq)))

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

def write_merlin_codebook(codebook_file:str, version:str, codebook_name:str, bit_names:list, 
        gene_names:list, transcript_names:str, barcode_str_list:list):
    '''Write a MERlin style codebook.'''
    with open(codebook_file, 'w') as f:
        # Write the header
        f.write(f'version, {version}\n')
        f.write(f'codebook_name, {codebook_name}\n')
        f.write(f'bit_names, {", ".join(bit_names)}\n')
        f.write('name, id, barcode\n')

        # Write the genes
        for i in range(len(gene_names)):
            f.write(f'{gene_names[i]}, {transcript_names[i]}, {barcode_str_list[i]}\n')


def load_transcriptome(transcripts_fasta_file:str, fpkm_tracking_file:str=None, style='ensembl'):
    '''Load the transcriptome into a pandas data frame.
    If the FPKM tracking file is not provided, set all FPKMs to be 1.
    '''
    
    # Load the transcripts
    transcripts = load_fasta_into_df(transcripts_fasta_file)
    print(f'Loaded {transcripts.shape[0]} transcripts.')
    transcripts.rename(columns={'id':'transcript_id'}, inplace=True)
   
    if fpkm_tracking_file is not None:
        # Load the FPKMs
        fpkms = pd.read_csv(fpkm_tracking_file, sep='\t')
        print(f'Loaded FPKMs for {fpkms.shape[0]} transcripts of {len(pd.unique(fpkms["gene_id"]))} genes.')
        fpkms.rename(columns={'tracking_id':'transcript_id'}, inplace=True)
        
        # Merge the two data frames
        transcriptome = transcripts.merge(fpkms, how='inner', on='transcript_id')
        print(f'Kept {transcriptome.shape[0]} transcripts of {len(pd.unique(transcriptome["gene_id"]))} genes after merging.')
   
    else:
        transcripts['FPKM'] = [1] * transcripts.shape[0]
        
        
        # try to parse description column, if its ensembl file:
        if style == 'ensembl':
            # Extract gene_id and gene_short_name from the description column
            # The description column is a string with the format:
            # gene_id:ENSG00000123456 gene_symbol:MYC
            # We will use regex to extract the gene_id and gene_short_name
            # from the description column.
            # The regex pattern is:
            # gene_id:([a-zA-Z0-9\.]+) gene_symbol:([a-zA-Z0-9\.]+)
            # This will match the gene_id and gene_short_name in the description column.
            gene_id_list, gene_name_list = [], []
            print(len(transcripts))
            for _str in transcripts['description']:
                _gene_id_match = re.search('gene:([a-zA-Z0-9\.]+) ', _str)
                if _gene_id_match:
                    gene_id_list.append(_gene_id_match.groups()[0])
                else:
                    gene_id_list.append(None)
                _gene_name_match = re.search('gene_symbol:([a-zA-Z0-9\.]+) ', _str)
                if _gene_name_match:
                    gene_name_list.append(_gene_name_match.groups()[0])
                else:
                    gene_name_list.append(None)
                    
            transcripts['gene_id'] = pd.Series(gene_id_list, dtype=str)
            transcripts['gene_short_name'] = pd.Series(gene_name_list, dtype=str)
            transcriptome = transcripts
            return transcriptome 
        
        elif style == 'gencode':
            parsed_dicts = []
            # use default 
            for _str in transcripts['description']:
                _matches = re.search(gencode_regexp, _str)
                if _matches:
                    parsed_dicts.append(_matches.groupdict())
                else:
                    raise ValueError(f"{_str}")
        
            # assemble df
            parsed_df = pd.DataFrame(parsed_dicts)
            parsed_df['sequence'] = transcripts['sequence']
            # rename gene_name
            parsed_df.rename(columns={'gene_name':'gene_short_name'}, inplace=True)
            parsed_df['FPKM'] = transcripts['FPKM']
            print(f"return gencode style transcriptome")
            print(parsed_dicts[0])
            return parsed_df
    
    


def load_primers(forward_primer_file:str, reverse_primer_file:str):
    '''Load the primers from fasta files'''
    forward_primers = load_fasta_into_df(forward_primer_file, load_rc=True)
    reverse_primers = load_fasta_into_df(reverse_primer_file, load_rc=True)

    return forward_primers, reverse_primers 

