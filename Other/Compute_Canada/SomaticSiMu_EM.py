#%%
"""
Imports
"""
import pandas as pd
import glob
import collections
import random
import numpy as np
from multiprocessing import Pool    
import functools
from collections import defaultdict
import os
import sys
import multiprocessing
import time
from multiprocessing import Pool
import time
from functools import partial
import argparse

def abs_path(target_name, directory_level): 
    """
Returns absolute file path of target name in working directory.

Arguments:
    target_name (str): Name of file or folder to find.
    directory_level (str): Level of os search, either File or Folder.   
    """
    #Find the relative working directory of the script
    wk_dir = os.path.dirname(os.path.realpath('__file__'))
    
    if directory_level == "File":
        #Absolute file path
        for root, dirs, files in os.walk(wk_dir):
            for name in files:
                if target_name == name:
                    target_path = (os.path.abspath(os.path.join(root, name))) 
             
    #Absolute file path
    if directory_level == "Directory":
        for root, dirs, files in os.walk(wk_dir):
            for name in dirs:
                if target_name == name:
                    target_path = (os.path.abspath(os.path.join(root, name))) 
    
    return target_path



def count_kmers(data, k):
    """
Count kmers in a sequence. Outputs dictionary.

Arguments:
    data (str): Input sequence.
    k (int): Length of k-mer.
    """
    global d
    d = collections.defaultdict(int)
    for i in range(len(data)-(k-1)):
        d[data[i:i+k]] +=1
    for key in list(d.keys()):
        if "N" in key:
            del d[key]
    return d
 
    

def read_fasta(fp):
    """
Reads fasta file. Outputs sequence.

Arguments:
    fp (str): Input fasta file. 
    """  
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))





def seq_slice(sequence_abs_path, slice_start=None, slice_end=None):
    """
Set slice of genome sequence (inclusive of start and end of slice) as a string (high memory cost...)

Arguments:
    sequence_abs_path (str): Absolute path of sequence
    slice_start (int): Start of string slice
    slice_end (int): End of string slice, inclusive
    """
    with open(sequence_abs_path) as fp:
        for name, seq in read_fasta(fp):
            sequence = str(seq)
            sample = sequence[slice_start:slice_end]
    
    return sample

"""
Codon dictionary
"""
codon_dict = { 
       "CYS": ["TGT", "TGC"], 
       "ASP": ["GAT", "GAC"], 
       "SER": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"], 
       "GLN": ["CAA", "CAG"], 
       "MET": ["ATG"], 
       "ASN": ["AAC", "AAT"], 
       "PRO": ["CCT", "CCG", "CCA", "CCC"], 
       "LYS": ["AAG", "AAA"], 
       "STOP": ["TAG", "TGA", "TAA"], 
       "THR": ["ACC", "ACA", "ACG", "ACT"], 
       "PHE": ["TTT", "TTC"], 
       "ALA": ["GCA", "GCC", "GCG", "GCT"], 
       "GLY": ["GGT", "GGG", "GGA", "GGC"], 
       "ILE": ["ATC", "ATA", "ATT"], 
       "LEU": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"], 
       "HIS": ["CAT", "CAC"], 
       "ARG": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"], 
       "TRP": ["TGG"], 
       "VAL": ["GTA", "GTC", "GTG", "GTT"], 
       "GLU": ["GAG", "GAA"], 
       "TYR": ["TAT", "TAC"]} 

"List of Exogenous Signatures in PCAWG Dataset"
signatures = {"SBS1": "Age",
              "SBS4": "Tobacco",
              "SBS7b": "UV",
              "SBS11": "Temozolomide",
              "SBS18": "Reactive Oxygen Species",
              "SBS22": "Aristolochic Acid",
              "SBS24": "Aflatoxin",
              "SBS25": "Chemotherapy",
              "SBS32": "Azathioprine",
              "SBS42": "Haloalkane",
              "SBS88": "Colibactin",
              "SBS90": "Duocarmycin"
               }




def syn_codon_dict(codon_dict):
    """
Synonymous codon dictionary

Arguments:
    codon_dict (dict)
    """
    syn_codon_dict = {}
    for codon_set in list(codon_dict.values()):
        for codon in codon_set:
            syn_codon = codon_set[:]  
            syn_codon.remove(codon)
            syn_codon_dict.setdefault(codon, syn_codon)
    return syn_codon_dict
    


def ref_data(input_name, num=False, prop=False):
    """
Reads in reference datasets from PCAWG studies.

Arguments:
    input_name (str): Name of csv file 
    num (bool): Set true to return dataframe of PCAWG SigProfiler Signatures in Samples 
    prop (bool): Set true to return dataframe of SigProfiler Signature Contributions
    
Note: Must set only one of num or prop arguments to True. 
    """
    if prop == num:
        raise ValueError("Please set only one argument between the arguments num and prop to True.")
        
    if num is True:
        
        num_file_path = abs_path(input_name, "File")
        num_data = pd.read_csv(num_file_path)
        num_data.sort_values(by='Cancer Types', axis=0, inplace=True)
        num_data.set_index(keys=['Cancer Types'], drop=False,inplace=True)
        
        return num_data

    if prop is True:
            
        prop_file_path = abs_path(input_name, "File")
        prop_data = pd.read_csv(prop_file_path)
        
        return prop_data
    

#Cell 1 Loaded
print("Cell 1 of 4 Loaded")

#%%

"""
Input sequence file path
"""
input_file_path = abs_path("Homo_sapiens.GRCh38.dna.chromosome.22.fasta", "File")

"""
Read in input reference files:
data = Pandas Dataframe of number of mutations attributed to each signature for each sample.
cancer_type = List of unique cancer types.
prop_data = Pandas Dataframe of proportion of 96 type mutation classification for each signature.
"""

"""
SBS input reference files
"""
sbs_num_file_path = abs_path("PCAWG_sigProfiler_SBS_signatures_in_samples.csv","File")
sbs_num_data = pd.read_csv(sbs_num_file_path)
sbs_num_data.sort_values(by='Cancer Types', axis=0, inplace=True)
sbs_num_data.set_index(keys=['Cancer Types'], drop=False,inplace=True)

sbs_prop_file_path = abs_path("sigProfiler_SBS_signatures.csv","File")
sbs_prop_data = pd.read_csv(sbs_prop_file_path)
  


"""
Expected Frequency Absolute Folder Paths for SBS, DBS and ID
"""
sbs_freq_folder_path = abs_path("SBS_Expected_Frequency", "Directory")


#Cell 2 Loaded
print("Cell 2 of 4 Loaded")


#%%


def sequence_index_dict(input_file_path, slice_start=None, slice_end=None, count=False):
    """
Pre-processing of the sequence into dictionary of format {kmer: [index1, index2, index3]}
High memory cost, only use SomaticSiMu up to 50 million base pairs of simulation

Arguments:
    input_file_path (str): Input file path of sequence to be preprocessed
    slice_start (int): Start of sequence slice for preprocessing
    slice_end (int): End of sequence slice for preprocessing
    count (bool): True returns counts of sequence, False returns indexes of sequence
    """
    sample = seq_slice(input_file_path, slice_start, slice_end)
    
    sample_index_dict = defaultdict(list)
    
    for seq_len in range(0, len(sample)-3):
        if 'N' in sample[seq_len:seq_len+3]:
            pass
        
        else:
            sample_index_dict[str(sample[seq_len:seq_len+3])].append(seq_len)
            

    if count is True:
        sample_count_dict = {key: len(value) for key, value in sample_index_dict.items()}
        return sample_count_dict
    
    else:
        return sample_index_dict

#Cell 3 Loaded
print("Cell 3 of 4 Loaded")

#%%



def somaticsimu_em(sample, sample_index_dict, input_file_path, signature, number_of_lineages):

    mutation_rate = np.random.normal(0.0005, 0.0001, 1)[0]
    #Create directories
    try:
        seq_directory = abs_path("Sequence", "Directory")
    except:
        try:
            os.mkdir("Sequence")
        except OSError:
            print ("Creation of the directory failed.")
            
    try:
        mut_mat_directory = abs_path("Metadata", "Directory")
    except:
        try:
            os.mkdir("Metadata")
        except OSError:
            print ("Creation of the directory failed.")
    
    
    #Number of mutations to simulate
    n_mutations = len(sample) * mutation_rate
    
    #Signature probabilities
    signature_df = sbs_prop_data.iloc[:, :2]
    signature_df[signature] = sbs_prop_data[signature]/sbs_prop_data[signature].sum()
    
    #Sample mutations
    mutation_array = np.random.multinomial(n_mutations, list(signature_df[signature]))
    
    #Simualted Mutation Counts 
    mutation_df = sbs_prop_data.iloc[:, :2]
    mutation_df['Count'] = mutation_array
    
     #Apply the mutations linearly 
    
    mutation_filtered_df = mutation_df[mutation_df.Count != 0]
    
    for i in range(len(mutation_filtered_df)):
        repeat  = mutation_filtered_df.iloc[i, 2]
        for j in range(repeat):
            position = random.choice(sample_index_dict[mutation_filtered_df.iloc[i, 1]])
            if sample[position] == mutation_filtered_df.iloc[i, 0][0]:
                sample = sample[:position] + mutation_filtered_df.iloc[i, 0][2] + sample[position+1:]
            else:
                "Wrong index!"
                
    #Make directories
    try:
        os.mkdir(seq_directory + "/" + signature)
    except OSError:
        print ("Creation of the directory failed.")
    
    try:
        os.mkdir(mut_mat_directory + "/" + signature)
    except OSError:
        print ("Creation of the directory failed.")
        

    #Write mutated sequence to fasta file
    sample_sequence_file = seq_directory + "/" + signature + "/" + signature + '_Sequence_' + str(number_of_lineages) + '.fasta'
    with open(sample_sequence_file, 'w+') as fasta_file:
        fasta_file.write(">" + str(signature) + "_Sequence_" + str(number_of_lineages) + "\n")
        fasta_file.write("")
        fasta_file.write(sample)
        
    #Write SBS mutation frequency tables to csv file
    sbs_freq_path = mut_mat_directory + "/" + signature + "/" + signature + '_Sequence_' + str(number_of_lineages)+ '_sbs_freq_table.csv'
    with open(sbs_freq_path, 'w+'):
        mutation_df.to_csv(sbs_freq_path, index=False)
            

            

#Cell 4 Loaded
print("Cell 4 of 4 Loaded")

# %%

def main(): 
    
    # Initiate the parser
    parser = argparse.ArgumentParser(description="SomaticSiMu Signature Tutorial Parameters")
    
    # Add long and short argument
    parser.add_argument("--generation", "-g", help="number of simulated sequences", default=10)
    parser.add_argument("--signature", "-s", help="cancer type")
    parser.add_argument("--reference", "-r", help="full file path of reference sequence used as input for the simulation", default=abs_path("Homo_sapiens.GRCh38.dna.chromosome.22.fasta", "File"))
    # Read arguments from the command line
    args = parser.parse_args()
  
   

   # try:
    if args.signature not in signatures.keys():
        print('Signature not found in existing types.')
    else:
        iterable = range(1, int(args.generation) + 1)
        pool = multiprocessing.Pool(4)
        starttime= time.time()
        
        signature = str(args.signature)
        input_file_path = str(args.reference)
        #Reference sequence
        sample = seq_slice(input_file_path, slice_start=None, slice_end=None)
        #Map sequence
        sample_index_dict = sequence_index_dict(input_file_path, slice_start=None, slice_end=None, count=False)
                    
        func = partial(somaticsimu_em, sample, sample_index_dict, input_file_path, signature)
    
        pool.map(func , iterable )
        pool.close()
        pool.join()
        print('That took {} seconds'.format(time.time() - starttime))
   # except:
   #     print("Faile")
        
    
        
if __name__ ==  '__main__':
    main()






