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

"List of Cancer Types in PCAWG Dataset"
cancer_type_list = ['Biliary-AdenoCA',
 'Bladder-TCC',
 'Bone-Benign',
 'Bone-Epith',
 'Bone-Osteosarc',
 'Breast-AdenoCA',
 'Breast-DCIS',
 'Breast-LobularCA',
 'CNS-GBM',
 'CNS-Medullo',
 'CNS-Oligo',
 'CNS-PiloAstro',
 'Cervix-AdenoCA',
 'Cervix-SCC',
 'ColoRect-AdenoCA',
 'Eso-AdenoCA',
 'Head-SCC',
 'Kidney-ChRCC',
 'Kidney-RCC',
 'Liver-HCC',
 'Lung-AdenoCA',
 'Lung-SCC',
 'Lymph-BNHL',
 'Lymph-CLL',
 'Myeloid-AML',
 'Myeloid-MDS',
 'Myeloid-MPN',
 'Ovary-AdenoCA',
 'Panc-AdenoCA',
 'Panc-Endocrine',
 'Prost-AdenoCA',
 'Skin-Melanoma',
 'SoftTissue-Leiomyo',
 'SoftTissue-Liposarc',
 'Stomach-AdenoCA',
 'Thy-AdenoCA',
 'Uterus-AdenoCA']


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
            
        prop_file_path = abs_path(input_name, input_type )
        prop_data = pd.read_csv(prop_file_path)
        
        return prop_data
    


def kmer_ref_count(kmer_length_start, kmer_length_end, ref_seq_name = "GRCh38.p13.genome.fasta", output_directory = "kmer_ref_count"):
    """
Count number of kmers of specified length in a fasta sequence. 
Used to count number of kmers in GChr.38

Arguments:
    kmer_length_start (int): Start counting kmers at this index
    kmer_length_end (int): Stop counting kmers at this index
    ref_seq_name (str): Name of file to count kmers from, default is the reference genome "GRCh38.p13.genome.fasta"
    output_directory (str): Name of folder to output kmer counts, default is the folder "kmer_ref_count"
    """
    input_file_path = abs_path(ref_seq_name, "File")
    
    with open(input_file_path) as fp:
        for name, seq in read_fasta(fp):
            for kmer_length in list(range(kmer_length_start, kmer_length_end+1)):
         
                count_kmers((str(seq)), kmer_length)
                count_df = pd.DataFrame.from_dict(d.keys())
                count_df['count'] = d.values()
                count_df.columns = [str(kmer_length), 'count']
        
                count_file_path = abs_path(output_directory, "Directory") + "/" + str(kmer_length)+'-mer_' + name + '.csv'
            
                with open(count_file_path, 'w'):
                    count_df.to_csv(count_file_path)
                    

#Cell 1 Loaded
print("Cell 1 of 6 Loaded")

#%%

"""
Input sequence file path
"""
input_file_path = abs_path("Homo_sapiens.GRCh38.dna.chromosome.22.fasta", "File")
# print(input_sequence)
#input_sequence = seq_slice(input_file_path, 26100000, 27100000)
sig_weights =  pd.read_csv(abs_path("timelines_sigWeights.txt", "File"), sep='\t')
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
DBS input reference files
"""
dbs_num_file_path = abs_path("PCAWG_sigProfiler_DBS_signatures_in_samples.csv","File")
dbs_num_data = pd.read_csv(dbs_num_file_path)
dbs_num_data.sort_values(by='Cancer Types', axis=0, inplace=True)
dbs_num_data.set_index(keys=['Cancer Types'], drop=False,inplace=True)

dbs_prop_file_path = abs_path("sigProfiler_DBS_signatures.csv","File")
dbs_prop_data = pd.read_csv(dbs_prop_file_path)


"""
Insertion/Deletion input reference files
"""
id_num_file_path = abs_path("PCAWG_SigProfiler_ID_signatures_in_samples.csv","File")
id_num_data = pd.read_csv(id_num_file_path)
id_num_data.sort_values(by='Cancer Types', axis=0, inplace=True)
id_num_data.set_index(keys=['Cancer Types'], drop=False,inplace=True)

id_prop_file_path = abs_path("sigProfiler_ID_signatures.csv","File")
id_prop_data = pd.read_csv(id_prop_file_path)


"""
Expected Frequency Absolute Folder Paths for SBS, DBS and ID
"""
sbs_freq_folder_path = abs_path("SBS_Expected_Frequency", "Directory")
dbs_freq_folder_path = abs_path("DBS_Expected_Frequency", "Directory")
id_freq_folder_path = abs_path("ID_Expected_Frequency", "Directory")

#Cell 2 Loaded
print("Cell 2 of 6 Loaded")


#%%


def sig_proportion(cancer_type, num_data, std_outlier):
    """
Find the proportions of different combinations of signatures among all samples
Output proportion of samples with each combination of signatures. 

Arguments:
    cancer_type (str): Cancer type 
    num_data (dataframe): Signature data
    std_outlier (float): Parameter to exclude outliers by degrees of standard deviation 
    """
    #Setup dataframe with signature data
    n_mutation_df = num_data.loc[cancer_type,"Sample Names":]
    n_mutation_df.set_index("Sample Names", inplace=True)
    del n_mutation_df['Accuracy']
    
    cols = n_mutation_df.columns
    
    #Find contributing signature combinations for each sequence
    sig_combos = list(n_mutation_df.apply(lambda x: x > 0).apply(lambda x: list(n_mutation_df.columns[x.values]), axis=1))
  
    #Find the set of unique contributing signature combinations 
    unique_sig_combos = [list(x) for x in set(tuple(x) for x in sig_combos)]
    
    
    #All outliers found for each sequence in no particular order
    net_outlier_sample = []
    
    for index in unique_sig_combos:
        combo_df = n_mutation_df[index][(n_mutation_df[index] != 0).all(1)]
        filter_combo_df = combo_df[~(np.abs(combo_df-combo_df.mean()) > (std_outlier*combo_df.std()))]

        for sig in index:
            net_outlier_sample.append(filter_combo_df[filter_combo_df[sig].isnull()].index.tolist())
    
    #Set of unique outlier sequences
    set_outlier = list(set(x for l in net_outlier_sample for x in l))
    
    #print("Removed the following outliers from the signature combination probability matrix:" + str(set_outlier))

    #Drop outlier sequences 
    n_mutation_df.drop(set_outlier,inplace=True)
    
    #Dataframe of the proportion of a set of signatures contribution out of all possible combinations
    bt=n_mutation_df.apply(lambda x: x > 0)
    
    sig_prop = ( bt.apply( lambda x: (cols[x.values]) , axis=1).value_counts() ) / len(bt)
    
    return sig_prop



def outlier_detection(cancer_type, num_data, std_outlier):
    """
Identifies outliers in the dataframe, returns list of outliers

Arguments;
    cancer_type (str): Cancer type
    num_data (dataframe): Signature data
    std_outlier: Parameter to exclude outliers by degrees of standard deviation 
    """ 
    #Setup dataframe with signature data
    n_mutation_df = num_data.loc[cancer_type,"Sample Names":]
    n_mutation_df.set_index("Sample Names", inplace=True)
    del n_mutation_df['Accuracy']
    
    #Find contributing signature combinations for each sequence
    sig_combos = list(n_mutation_df.apply(lambda x: x > 0).apply(lambda x: list(n_mutation_df.columns[x.values]), axis=1))
  
    #Find the set of unique contributing signature combinations 
    unique_sig_combos = [list(x) for x in set(tuple(x) for x in sig_combos)]
    
    #All outliers found for each sequence in no particular order
    net_outlier_sample = []
    
    for index in unique_sig_combos:
        combo_df = n_mutation_df[index][(n_mutation_df[index] != 0).all(1)]
        filter_combo_df = combo_df[~(np.abs(combo_df-combo_df.mean()) > (std_outlier*combo_df.std()))]

        for sig in index:
            net_outlier_sample.append(filter_combo_df[filter_combo_df[sig].isnull()].index.tolist())
    
    #Set of unique outlier sequences
    set_outlier = list(set(x for l in net_outlier_sample for x in l))
    
    return set_outlier



def sig_freq(num_data, mut_type):
    """
Output one file for each sample. File contains expected number of 96, 78 and 24 type mutations for 
single base substitution (SBS), double base substitution (DBS) and insertion/deletion (ID) in the sample. 

Arguments:
    num_data (dataframe): Signature data
    mut_type (str): "SBS", "DBS", "ID"
    """
    
    if mut_type == "SBS":
        
        temp_data3 = num_data.copy().loc[:, 'Accuracy':]
        del temp_data3['Accuracy']
        cols = temp_data3.columns
        for row in range(len(temp_data3)):
            new_df = sbs_prop_data.loc[:,'Type':'SubType']
            for column in range(len(cols)):
                frequency = temp_data3.iloc[row, column]
                new_df[cols[column]] = [(frequency*item) for item in sbs_prop_data[cols[column]]]
            ct = temp_data3.index[row]
            sample_id = sbs_num_data.iloc[row, 1]
            
            with open(abs_path("SBS_Expected_Frequency", "Directory") + "/SBS_" + ct + '_' + sample_id + '.csv', 'w'):
                new_df.to_csv(abs_path("SBS_Expected_Frequency", "Directory") + "/SBS_" + ct + '_' + sample_id + '.csv')

    if mut_type == "DBS":
        
        temp_data3 = num_data.loc[:, 'Accuracy':]
        del temp_data3['Accuracy']
        cols = temp_data3.columns
        for row in range(len(temp_data3)):
            new_df = dbs_prop_data.loc[:,'Mutation Type'].to_frame() 
            for column in range(len(cols)):
                frequency = temp_data3.iloc[row, column]
                new_df[cols[column]] = [(frequency*item) for item in dbs_prop_data[cols[column]]]         
            ct = temp_data3.index[row]
            sample_id = dbs_num_data.iloc[row, 1]   
            with open(abs_path("DBS_Expected_Frequency", "Directory") + "/DBS_" + ct + '_' + sample_id + '.csv', 'w'):
                new_df.to_csv(abs_path("DBS_Expected_Frequency", "Directory") + "/DBS_" + ct + '_' + sample_id + '.csv')

    if mut_type == "ID":
    
        temp_data3 = num_data.loc[:, 'Accuracy':]
        del temp_data3['Accuracy']
        cols = temp_data3.columns
        for row in range(len(temp_data3)):
            new_df = id_prop_data.loc[:,'Mutation Type':'ID_type']
            for column in range(len(cols)):
                frequency = temp_data3.iloc[row, column]
                new_df[cols[column]] = [(frequency*item) for item in id_prop_data[cols[column]]]         
            ct = temp_data3.index[row]
            sample_id = id_num_data.iloc[row, 1]   
            with open(abs_path("ID_Expected_Frequency", "Directory") + "/ID_" + ct + '_' + sample_id + '.csv', 'w'):
                new_df.to_csv(abs_path("SBS_Expected_Frequency", "Directory") + "/ID_" + ct + '_' + sample_id + '.csv')

    if mut_type not in ["SBS", "DBS", "ID"]:
        raise ValueError("It looks like you did not input an appropriate mut_type argument (SBS, DBS, ID)")

    else:
        print("Finished producing expected mutation frequences of " + str(mut_type) + "type to the " + str(mut_type) + "_Expected_Frequency directory.")



def mut_mean(cancer_type, num_data, freq_folder_path, mut_type, std_outlier):
    """
Output average number of 96 type mutations among all samples for specified cancer type.

Arguments:
    cancer_type (str): Cancer type
    num_data (dataframe): Signature data
    freq_folder_path (str): Folder with expected mutation frequencies for each sample
    mut_type (str): "SBS", "DBS", "ID"
    std_outlier (float): Parameter to exclude outliers by degrees of standard deviation 
    """
    if mut_type == "SBS":
         
        outlier_list = outlier_detection(cancer_type, num_data, std_outlier)
        outlier_paths = [abs_path("SBS_" + cancer_type + "_" + outlier + ".csv", "File") for outlier in outlier_list]
        file_paths = glob.glob(str(sbs_freq_folder_path) + "/SBS_" + cancer_type + '_*' + '.csv')
        filtered_file_paths = [x for x in file_paths if x not in outlier_paths]
        

        print("Removed " + str(len(file_paths)-len(filtered_file_paths)) + " outliers from the average mutation frequency matrix")
        print(outlier_list)
        
        og_df = pd.read_csv(filtered_file_paths[0])
        mean_df = og_df.loc[:, 'SBS1':]
    
        for file_index in range(1, len(file_paths)):
            raw_mean = pd.read_csv(file_paths[file_index])
            slice_read_df = raw_mean.loc[:, 'SBS1':]
            mean_df += slice_read_df
            
        df_averages = mean_df/len(file_paths) 
        mut_subtypes = og_df.loc[:, 'Type': 'SubType']
        average_freq_df = pd.concat([mut_subtypes, df_averages], axis=1)
        
        return average_freq_df
  
    if mut_type == "DBS":
        
        outlier_list = outlier_detection(cancer_type, num_data, std_outlier)
        outlier_paths = [abs_path("DBS_" + cancer_type + "_" + outlier + ".csv", "File") for outlier in outlier_list]
        file_paths = glob.glob(str(dbs_freq_folder_path)+ "/DBS_" + cancer_type + '_*' + '.csv')
        filtered_file_paths = [x for x in file_paths if x not in outlier_paths]
        
        print("Dropped " + str(len(file_paths)-len(filtered_file_paths)) + " outliers!")
        print("")
        print(outlier_list)
        
        og_df = pd.read_csv(filtered_file_paths[0], index_col=0)
        mean_df = og_df.loc[:, 'DBS1':]
        
        for file_index in range(1, len(file_paths)):
            raw_mean = pd.read_csv(file_paths[file_index])
            slice_read_df = raw_mean.loc[:, 'DBS1':]
            mean_df += slice_read_df
            
        dbs_df_averages = mean_df/len(file_paths) 
        mut_subtypes = og_df.loc[:,'Mutation Type'].to_frame() 
        average_freq_df = pd.concat([mut_subtypes, dbs_df_averages], axis=1)
        
        return average_freq_df

    if mut_type == "ID":
        
        outlier_list = outlier_detection(cancer_type, num_data, std_outlier)
        outlier_paths = [abs_path("ID_" + cancer_type + "_" + outlier + ".csv", "File") for outlier in outlier_list]
        file_paths = glob.glob(str(id_freq_folder_path) + "/ID_" + cancer_type + '_*' + '.csv')
        filtered_file_paths = [x for x in file_paths if x not in outlier_paths]
        
        print("Dropped " + str(len(file_paths)-len(filtered_file_paths)) + " outliers!")
        print("")
        print(outlier_list)
        
        og_df = pd.read_csv(filtered_file_paths[0], index_col=0)
        mean_df = og_df.loc[:, 'ID1':]
        
        for file_index in range(1, len(file_paths)):
            raw_mean = pd.read_csv(file_paths[file_index])
            slice_read_df = raw_mean.loc[:, 'ID1':]
            mean_df += slice_read_df
            
        id_df_averages = mean_df/len(file_paths) 
        mut_subtypes = og_df.loc[:,'Mutation Type':'ID_type']
        average_freq_df = pd.concat([mut_subtypes, id_df_averages], axis=1)
        
        return average_freq_df
  
    if mut_type not in ["SBS", "DBS", "ID"]:
        raise ValueError("It looks like you did not input an appropriate mut_type argument (SBS, DBS, ID)")

    else: 
        print("Completed calculations of average " + str(mut_type) + "mutation expected frequency for " + str(cancer_type))

 

def genome_kmer_count(kmer_length, kmer_count_folder_path = abs_path("kmer_ref_count", "Directory")):
    """
Output genome kmer counts (Chr1-22, X and Y) of reference genome GChr38

Arguments:
    kmer_length: Counts of k-mers in reference sequence at specified length
    kmer_count_folder_path: Folder path of kmer counts, default is folder of reference genome kmer counts 
    """
    chr_index = [i for i in range(kmer_length + 1 ,23)]
    chr_index.append('X')
    chr_index.append('Y')
    
    chr_1 = pd.read_csv(kmer_count_folder_path + '/' + str(kmer_length) + '-mer/' + str(kmer_length) + '-mer_chr1.csv')
    total_kmer_count = pd.read_csv(kmer_count_folder_path + '/' + str(kmer_length) + '-mer/' + str(kmer_length) + '-mer_chr1.csv').loc[:, 'count']
    
    for file_name in chr_index:
        read_in = pd.read_csv(kmer_count_folder_path + '/' + str(kmer_length) + '-mer/' + str(kmer_length) + '-mer_chr' + str(file_name) + '.csv').loc[:, 'count']
        total_kmer_count += read_in
 
    return pd.concat([chr_1.iloc[:,1], total_kmer_count], axis=1)


#Cell 3 Loaded
print("Cell 3 of 6 Loaded")

#%%

"""
Total Kmer counts for lengths 1-6 of GChr.38
"""
kmer_count_folder_path = abs_path("kmer_ref_count", "Directory")
#/Users/DavidChen/Desktop/Project/kmer_ref_count

kmer_1_count = genome_kmer_count(1, kmer_count_folder_path)
kmer_1_count.columns = ['kmer', 'count']
kmer_2_count = genome_kmer_count(2, kmer_count_folder_path)
kmer_2_count.columns = ['kmer', 'count']
kmer_3_count = genome_kmer_count(3, kmer_count_folder_path)
kmer_3_count.columns = ['kmer', 'count']
kmer_4_count = genome_kmer_count(4, kmer_count_folder_path)
kmer_4_count.columns = ['kmer', 'count']
kmer_5_count = genome_kmer_count(5, kmer_count_folder_path)
kmer_5_count.columns = ['kmer', 'count']
kmer_6_count = genome_kmer_count(6, kmer_count_folder_path)
kmer_6_count.columns = ['kmer', 'count']

#Cell 4 Loaded
print("Cell 4 of 6 Loaded")

#%%
def mut_prob(cancer_type, prop_data, num_data, freq_folder_path, mut_type, std_outlier):
    """
Output dictionary of probability of mutation for any given k-mer in whole genome sequence

Arguments:
    cancer_type (str): Cancer type
    prop_data (dataframe): Signature proportion data
    num_data (dataframe): Signature frequency data
    freq_folder_path (str): Folder with expected mutation frequencies for each sample
    mut_type (str): "SBS", "DBS", "ID"
    std_outlier (float): Parameter to exclude outliers by degrees of standard deviation 
    """
    if mut_type == "SBS":
        
        sbs_mut = mut_mean(cancer_type, num_data, freq_folder_path, mut_type, std_outlier)
        sbs_mut['total'] = sbs_mut['SubType'].map(kmer_3_count.set_index('kmer')['count'])
        sbs_mut_prob = sbs_mut.loc[:,'SBS1':'SBS60'].div(sbs_mut.total, axis=0)
        sbs_index_mut_prob = pd.concat([sbs_mut[['Type', 'SubType']], sbs_mut_prob], axis=1)
        sbs_index_mut_prob.sort_values(['SubType','Type']).set_index(['SubType','Type'],inplace=True)
    
        return sbs_index_mut_prob
    
    if mut_type == "DBS":

        dbs_mut = mut_mean(cancer_type, num_data, freq_folder_path, mut_type, std_outlier)
        dbs_mut['2mer_index'] = dbs_mut['Mutation Type'].apply(lambda x: x[:2])
        dbs_mut['total'] = dbs_mut['2mer_index'].map(kmer_2_count.set_index('kmer')['count'])
        dbs_mut_prob = dbs_mut.loc[:,'DBS1':'DBS11'].div(dbs_mut['total'], axis=0)
        dbs_index_mut_prob = pd.concat([dbs_mut[['Mutation Type', '2mer_index']], dbs_mut_prob], axis=1)
        dbs_index_mut_prob.sort_values(['Mutation Type']).set_index(['Mutation Type'],inplace=True)
   
        return dbs_index_mut_prob
    
    if mut_type == "ID":

        id_mut = mut_mean(cancer_type, num_data, freq_folder_path, mut_type, std_outlier)
        id_mut['Index'] = id_prop_data['Index']
        
        id_kmer_count_1_6 = kmer_1_count.append([kmer_2_count,
                                                 kmer_3_count,
                                                 kmer_4_count,
                                                 kmer_5_count,
                                                 kmer_6_count])
        
        id_mut['total'] = id_mut['Index'].map(id_kmer_count_1_6.set_index('kmer')['count'])
        
        id_mut_prob = id_mut.loc[:,'ID1':'ID17'].div(id_mut['total'], axis=0)
        id_index_mut_prob_1_6 = pd.concat([id_mut.loc[:, 'Mutation Type':'ID_type'], id_mut_prob], axis=1)
        id_index_mut_prob_1_6.sort_values(['Mutation Type']).set_index(['Mutation Type'],inplace=True)
      
        return id_index_mut_prob_1_6
    
    if mut_type not in ["SBS", "DBS", "ID"]:
        raise ValueError("It looks like you did not input an appropriate mut_type argument (SBS, DBS, ID)")
   
    else:
        print("Completed calculations of " + str(mut_type) + "mutation probabilities for " + str(cancer_type))
            


def sequence_index_dict(input_file_path, slice_start, slice_end, count=False):
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
    
    for seq_len in range(0, len(sample)-6):
        if 'N' in sample[seq_len:seq_len+6]:
            pass
        
        else:
            sample_index_dict[str(sample[seq_len:seq_len+6])].append(seq_len)
            

    if count is True:
        sample_count_dict = {key: len(value) for key, value in sample_index_dict.items()}
        return sample_count_dict
    
    else:
        return sample_index_dict
    
 #Cell 5 Loaded
print("Cell 5 of 6 Loaded")

#%%
def somatic_sim(cancer_type, reading_frame, std_outlier, simulation_type, sequence_abs_path, slice_start, slice_end,power, syn_rate, non_syn_rate, number_of_lineages):
 
    #INPUT CHECKS
    if reading_frame not in [1, 2, 3]:
        print("reading_frame argument must be 1, 2, or 3.")
        sys.exit()
        
    if cancer_type not in cancer_type_list:
        print("cancer_type argument must be in " + str(cancer_type_list))
        sys.exit()
        
    if isinstance(std_outlier, int) == False:
        print("std_outlier argument must be int type")
        sys.exit()
    
    if isinstance(number_of_lineages, int) == False:
        print("number_of_lineages argument must be int type")
        sys.exit()
        
    if isinstance(simulation_type, str) == False:
        print("simulation_type argument must be int type")
        sys.exit()
        
    if isinstance(slice_start, int) == False:
        print("slice_start argument must be int type")
        sys.exit()
    
    if isinstance(slice_end, int) == False:
        print("slice_end argument must be int type")
        sys.exit()
        
    if isinstance(power, int) == False:
        print("power argument must be int type")
        sys.exit()
        
    if isinstance(syn_rate, int) == False:
        print("syn_rate argument must be int type")
        sys.exit()
        
    if isinstance(non_syn_rate, int) == False:
        print("non_syn_rate argument must be int type")
        sys.exit()
        
    if sequence_abs_path.split('.')[-1].lower() != "fasta":
        print("sequence_absolute_path must be a fasta file")
        sys.exit()
        
    #Set random seed 
    #Necessary for each iteration since child worker must have different numpy seed in multiprocessing
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
    
    #Output directory paths for mutated sequence and SBS, DBS, Insertion and Deletion frequency matrices
    seq_directory = abs_path("Sample", "Directory")
    mut_mat_directory = abs_path("Frequency_Table", "Directory")
    mut_metadata_directory = abs_path("Mutation_Metadata", "Directory")
    mut_sigs_directory = abs_path("Signature_Combinations", "Directory")
    
    
    #Setup probabilities into processed dataframes
    print('Lineage ' + str(number_of_lineages))
    
    #Types of mutations for a specific lineage
    insertion_dict={}
    deletion_dict={}
    sbs_dict={}
    dbs_dict={}
    
    #Initialize one iteration of the sequence to be mutated
    sample_seq = seq_slice(sequence_abs_path, slice_start, slice_end)
    
    #Synonymous codons 
    syn_codon = syn_codon_dict(codon_dict)
    
    if simulation_type == "end":
        
        #Select combination of signatures and set up mutation probabilities for indels
        id_prop = sig_proportion(cancer_type, id_num_data, std_outlier)
        id_sig_combination = random.choices(list(id_prop.index), list(id_prop.values))
        
        #Select combination of signatures and set up mutation probabilities for single base substitution
        sbs_prop = sig_proportion(cancer_type, sbs_num_data, std_outlier)
        sbs_sig_combination = random.choices(list(sbs_prop.index), list(sbs_prop.values))  
        
        #Select combination of signatures and set up mutation probabilities for double base substitution
        dbs_prop = sig_proportion(cancer_type, dbs_num_data, std_outlier)
        dbs_sig_combination = random.choices(list(dbs_prop.index), list(dbs_prop.values))  
        
        #Index of each kmer
        sample_index_dict = defaultdict(list)
    
        for seq_len in range(0, len(sample_seq)-6):
            if 'N' in sample_seq[seq_len:seq_len+6]:
                pass
    
            else:
                sample_index_dict[str(sample_seq[seq_len:seq_len+6])].append(seq_len)
        
        #Count of each kmer
        sample_count_dict = {key: len(value) for key, value in sample_index_dict.items()}
        
        #Normalization of mutation probabilities to whole genome burden
        ref_dir = abs_path("6-mer",  "Directory")
        kmer_ref = (glob.glob(ref_dir+ "/6-mer_chr*"))
        
        kmer_count = pd.read_csv(kmer_ref[0], index_col=0)['count'].fillna(0)
        for i in kmer_ref[1:-1]:
            sample = pd.read_csv(i, index_col=0)['count'].fillna(0)
            kmer_count = kmer_count.add(sample, fill_value=0)
              
            
        kmer_reference_count_dict = dict(zip(pd.read_csv(kmer_ref[0], index_col=0)["6"], kmer_count))
        
        sample_expected_count_dict = dict(zip(pd.read_csv(kmer_ref[0], index_col=0)["6"], kmer_count))
        for key in sample_expected_count_dict:
            sample_expected_count_dict[key] = sample_expected_count_dict[key] * (len(sample_seq) / sum(list( kmer_reference_count_dict.values()))) 
        
        
        normalize_constant = {k: float(sample_count_dict[k])/sample_expected_count_dict[k] for k in sample_count_dict}
        normalized_sample_count_dict = {k: int(sample_count_dict[k]/normalize_constant[k]) for k in sample_count_dict}

        print('Simulating Tumour Stage End of Lineage ' + str(number_of_lineages))
    
        #Types of mutations for a specific generation  of a specific lineage
        gen_insertion_dict={}
        gen_deletion_dict={}
        gen_sbs_dict={}
        gen_dbs_dict={}  
        
        #Types of mutations for a specific lineage and generation to be outputted as a frequency matrix
        output_insertion_dict={}
        output_deletion_dict={}
        output_sbs_dict={}
        output_dbs_dict={}
    
        id_prob = mut_prob(cancer_type, id_prop_data, id_num_data, id_freq_folder_path, "ID", std_outlier)
        id_mut_sig_df = pd.concat([id_prob.loc[:47,'Mutation Type':'ID_type'], id_prob[id_sig_combination[0]].dropna()], axis=1)
            
        id_total_prob_df = id_mut_sig_df.iloc[:,:5]
        id_total_prob_df['total_prob']= id_mut_sig_df.loc[:, 'ID_type':].sum(axis=1)
        id_total_prob_df['total_prob'] = id_total_prob_df['total_prob'].multiply(power)
           
        id_sorted_df = id_total_prob_df.sort_values(['ID_type','Index','Size'])
        id_sorted_df.set_index(['Index'], inplace=True)
            
        deletion_df = id_sorted_df.iloc[:12, :]
        insertion_df = id_sorted_df.iloc[12:24, :]
            
        print("Insertion and Deletion mutation probability matrix set up.")
        print(id_sig_combination)
        
       
        sbs_prob = mut_prob(cancer_type, 
                           sbs_prop_data,
                           sbs_num_data, 
                           sbs_freq_folder_path,
                           "SBS",
                           std_outlier)
        
        sbs_mut_sig_df = pd.concat([sbs_prob.loc[:,'Type':'SubType'], sbs_prob[sbs_sig_combination[0]]], axis=1)
        
        sbs_total_prob_df = sbs_mut_sig_df.iloc[:,:2]
        sbs_total_prob_df['total_prob'] = sbs_mut_sig_df.iloc[:, 2:].sum(axis=1)
            
        sbs_total_prob_df['total_prob'] = sbs_total_prob_df[sbs_total_prob_df.select_dtypes(include=['number']).columns] * power
        
        sbs_sorted_df = sbs_total_prob_df.sort_values(['SubType'])
        sbs_sorted_df.set_index(['SubType'], inplace=True)
        
        print("Single Base Substitution mutation probability matrix set up.")
        print(sbs_sig_combination)
            
        dbs_prob = mut_prob(cancer_type, 
                           dbs_prop_data,
                           dbs_num_data, 
                           dbs_freq_folder_path,
                           "DBS",
                           std_outlier)
           
        dbs_mut_sig_df = pd.concat([dbs_prob.loc[:,'Mutation Type':'2mer_index'], dbs_prob[dbs_sig_combination[0]]], axis=1)
        
        dbs_total_prob_df = dbs_mut_sig_df.iloc[:,:2]   
        dbs_total_prob_df['total_prob']= dbs_mut_sig_df.iloc[:, 1:].sum(axis=1) 
        dbs_total_prob_df['total_prob'] = dbs_total_prob_df[dbs_total_prob_df.select_dtypes(include=['number']).columns] * power
        
        dbs_sorted_df =  dbs_total_prob_df.sort_values(['Mutation Type'])
        dbs_sorted_df.set_index(['2mer_index'], inplace=True)   
        
        print("Double Base Substitution mutation probability matrix set up.")
        print(dbs_sig_combination)
    
        for keys in list(sample_index_dict.keys()):
            
            count=1
            for length in list(range(1, len(keys))):
                if keys[length] == keys[0]:
                    count += 1
                else:
                    break
      
            #Insertion
            if "G" in keys[:count] or "A" in keys[:count]:
                pass
            
            else:
                temp_i_prob = [insertion_df.loc[insertion_df['Size'] == count].loc[keys[:count],'total_prob']] 
                insertion_count = np.random.binomial(normalized_sample_count_dict[keys], temp_i_prob, 1) 
            
                for m in range(sum(i_number > 0 for i_number in list(insertion_count))):
                    insertion_index = np.random.choice(sample_index_dict[keys])
                    insertion_dict[insertion_index] = keys[0]
                    gen_insertion_dict[insertion_index] = keys[0]
                    output_insertion_dict[insertion_index] = keys[:count]
                
                
            #Deletion
            if "G" in keys[:count] or "A" in keys[:count]:
                pass
            
            else:
                temp_d_prob = [deletion_df.loc[deletion_df['Size'] == count].loc[keys[:count],'total_prob']] 
                deletion_count = np.random.binomial(normalized_sample_count_dict[keys], temp_d_prob, 1)  
            
                for n in range(sum(d_number > 0 for d_number in list(deletion_count))):
                    deletion_index = np.random.choice(sample_index_dict[keys])
                    deletion_dict[deletion_index] = keys[0]
                    gen_deletion_dict[deletion_index] = keys[0]
                    output_deletion_dict[deletion_index] = keys[:count]
        
            #SBS
            if "G" in keys[1] or "A" in keys[1]:
                pass
            
            else:
                sbs_types = list(sbs_sorted_df.loc[keys[:3],'Type']) + [None]
                sbs_prob = list(sbs_sorted_df.loc[keys[:3],'total_prob'])+ [1- sum(list(sbs_sorted_df.loc[keys[:3],'total_prob']))]             
                sbs = [sbs_num for sbs_num in list(np.random.choice(sbs_types,p=sbs_prob, size = normalized_sample_count_dict[keys])) if sbs_num]
                for sbs_value in sbs:
                    single_base_index = np.random.choice(sample_index_dict[keys])
                    sbs_dict[single_base_index] = sbs_value[2]
                    gen_sbs_dict[single_base_index] = sbs_value[2]
                    output_sbs_dict[single_base_index] = [sbs_value, sample_seq[single_base_index:single_base_index+3]]
                         
            #DBS   
            if keys[:2] in dbs_sorted_df.index:
                dbs_types = list(dbs_sorted_df.loc[keys[:2],'Mutation Type']) + [None]     
                dbs_prob = list(dbs_sorted_df.loc[keys[:2],'total_prob'])+ [1- sum(list(dbs_sorted_df.loc[keys[:2],'total_prob']))]
                dbs = [dbs_num for dbs_num in list(np.random.choice(dbs_types,p=dbs_prob, size = normalized_sample_count_dict[keys])) if dbs_num]
                for dbs_value in dbs:
                    double_base_index = np.random.choice(sample_index_dict[keys])
                    dbs_dict[double_base_index] = dbs_value[3:]
                    gen_dbs_dict[double_base_index] = dbs_value[3:]
                    output_dbs_dict[double_base_index] = dbs_value
                    
            else:
                pass
            
      
        #SBS Synonymous/Non-Synonymous
        #Reading frame starts at the first base
        if reading_frame == 1:
            
            #Apply syn/non-syn mutation rate for SBS
            for mut in list(gen_sbs_dict.keys()):
            
                read_index = mut%3
            
                #First base of codon
                if read_index == 1:
                    
                    if gen_sbs_dict[mut] + sample_seq[mut+1:mut+3] == sample_seq[mut:mut+3]:
                        pass
                    
                    if gen_sbs_dict[mut] + sample_seq[mut+1:mut+3] in syn_codon[sample_seq[mut:mut+3]]:
                        syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                    
                    else:
                        syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                                                   
                #Second base of codon     
                if read_index == 2:
                    
                    if sample_seq[mut-1] + gen_sbs_dict[mut] + sample_seq[mut+1] == sample_seq[mut-1:mut+2]:
                        pass
                    
                    if sample_seq[mut-1] + gen_sbs_dict[mut] + sample_seq[mut+1] in syn_codon[sample_seq[mut-1:mut+2]]:
                        syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                                                         
                    else:
                        syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                            
                #Third base of codon     
                if read_index == 0:
                    
                    if sample_seq[mut-2:mut] + gen_sbs_dict[mut] == sample_seq[mut-2:mut+1]:
                        pass
                    
                    if sample_seq[mut-2:mut] + gen_sbs_dict[mut] in syn_codon[sample_seq[mut-2:mut+1]]:
                        syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                                           
                    else:
                        syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                            
             
        #Reading frame starts at the second base
        if reading_frame == 2:
            
            #Apply syn/non-syn mutation rate for SBS
            for mut in list(gen_sbs_dict.keys()):
            
                read_index = mut%3
            
                #First base of codon
                if read_index == 2:
                    
                    if gen_sbs_dict[mut] + sample_seq[mut+1:mut+3] == sample_seq[mut:mut+3]:
                        pass
                    
                    if gen_sbs_dict[mut] + sample_seq[mut+1:mut+3] in syn_codon[sample_seq[mut:mut+3]]:
                        syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                                           
                    else:
                        syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                                                   
                    #Second base of codon     
                if read_index == 0:
                    
                    if sample_seq[mut-1] + gen_sbs_dict[mut] + sample_seq[mut+1] == sample_seq[mut-1:mut+2]:
                        pass
                    
                    if sample_seq[mut-1] + gen_sbs_dict[mut] + sample_seq[mut+1] in syn_codon[sample_seq[mut-1:mut+2]]:
                        syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                                                         
                    else:
                        syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                            
                #Third base of codon     
                if read_index == 1:
                    
                    if sample_seq[mut-2:mut] + gen_sbs_dict[mut] == sample_seq[mut-2:mut+1]:
                        pass
                    
                    if sample_seq[mut-2:mut] + gen_sbs_dict[mut] in syn_codon[sample_seq[mut-2:mut+1]]:
                        syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                                           
                    else:
                        syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                            
                                              
        #Reading frame starts at the third base
        if reading_frame == 3:
            
            #Apply syn/non-syn mutation rate for SBS
            for mut in list(gen_sbs_dict.keys()):
            
                read_index = mut%3
            
                #First base of codon
                if read_index == 0:
                    
                    if gen_sbs_dict[mut] + sample_seq[mut+1:mut+3] == sample_seq[mut:mut+3]:
                        pass
                    
                    if gen_sbs_dict[mut] + sample_seq[mut+1:mut+3] in syn_codon[sample_seq[mut:mut+3]]:
                        syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                                           
                    else:
                        syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                                                   
                #Second base of codon     
                if read_index == 1:
                    
                    if sample_seq[mut-1] + gen_sbs_dict[mut] + sample_seq[mut+1] == sample_seq[mut-1:mut+2]:
                        pass
                    
                    if sample_seq[mut-1] + gen_sbs_dict[mut] + sample_seq[mut+1] in syn_codon[sample_seq[mut-1:mut+2]]:
                        syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                                                         
                    else:
                        syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                        
                        if syn_mutation[0] == 2:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                            
                #Third base of codon     
                if read_index == 0:
                    
                    if sample_seq[mut-2:mut] + gen_sbs_dict[mut] == sample_seq[mut-2:mut+1]:
                        pass
                    
                    if sample_seq[mut-2:mut] + gen_sbs_dict[mut] in syn_codon[sample_seq[mut-2:mut+1]]:
                        syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                                           
                    else:
                        syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del sbs_dict[mut]
                            del output_sbs_dict[mut]
                            
        
        
        
        #DBS Synonymous/Non-Synonymous
    
        #Reading frame starts at the first base
        if reading_frame == 1:
            
            #Apply syn/non-syn mutation rate for DBS
            for mut in list(gen_dbs_dict.keys()):
            
                read_index = mut%3
            
                #First + second base of codon
                if read_index == 1:
                    
                    if gen_dbs_dict[mut] + sample_seq[mut+2] == sample_seq[mut:mut+3]:
                        pass
                    
                    if gen_dbs_dict[mut] + sample_seq[mut+2] in syn_codon[sample_seq[mut:mut+3]]:
                        syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del dbs_dict[mut]
                            del output_dbs_dict[mut]
                                           
                    else:
                        syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del dbs_dict[mut]
                            del output_dbs_dict[mut]
                                                   
                #Second + third base of codon     
                if read_index == 2:
                    
                    if sample_seq[mut-1] + gen_dbs_dict[mut]  == sample_seq[mut-1:mut+2]:
                        pass
                    
                    if sample_seq[mut-1] + gen_dbs_dict[mut] in syn_codon[sample_seq[mut-1:mut+2]]:
                        syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del dbs_dict[mut]
                            del output_dbs_dict[mut]
                                                         
                    else:
                        syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del dbs_dict[mut]
                            del output_dbs_dict[mut]
                                   
                            
                #Third base of codon and first base of next adjacent codon    
                if read_index == 0:
                    
                    try:
                    
                        if sample_seq[mut-2:mut] + gen_dbs_dict[mut][0] == sample_seq[mut-2:mut+1]:
                            pass
                            
                        if gen_dbs_dict[mut][1] + sample_seq[mut+2:mut+4] == sample_seq[mut+1:mut+4]:
                            pass
                        
                        #Check if first codon of original sequence is synonymous to mutated codon
                        if sample_seq[mut-2:mut] + gen_dbs_dict[mut][0] in syn_codon[sample_seq[mut-2:mut+1]]:
                            syn_mutation_codon_1 = random.choices([1,2],[syn_rate**0.5, 1-syn_rate**0.5])
                           
                        #Check if second codon of original sequence is synonymous to mutated codon
                        if gen_dbs_dict[mut][1] + sample_seq[mut+2:mut+4] in syn_codon[sample_seq[mut+1:mut+4]]:
                            syn_mutation_codon_2 = random.choices([1,2],[syn_rate**0.5, 1-syn_rate**0.5])
                            
                            if syn_mutation_codon_1[0] and syn_mutation_codon_2[0] == 1:
                                pass
                            
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
    
                           
                        #Both mutated codons are non-synonymous to their original sequence codon                           
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                    except:
                        pass
        
        #Reading frame starts at the second base
        if reading_frame == 2:
            
            #Apply syn/non-syn mutation rate for DBS
            for mut in list(gen_dbs_dict.keys()):
            
                read_index = mut%3
            
                #First + second base of codon
                if read_index == 2:
                    
                    if gen_dbs_dict[mut] + sample_seq[mut+2] == sample_seq[mut:mut+3]:
                        pass
                    
                    if gen_dbs_dict[mut] + sample_seq[mut+2] in syn_codon[sample_seq[mut:mut+3]]:
                        syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del dbs_dict[mut]
                            del output_dbs_dict[mut]
                                           
                    else:
                        syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del dbs_dict[mut]
                            del output_dbs_dict[mut]
                                                   
                #Second + third base of codon     
                if read_index == 0:
                    
                    if sample_seq[mut-1] + gen_dbs_dict[mut]  == sample_seq[mut-1:mut+2]:
                        pass
                    
                    if sample_seq[mut-1] + gen_dbs_dict[mut] in syn_codon[sample_seq[mut-1:mut+2]]:
                        syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del dbs_dict[mut]
                            del output_dbs_dict[mut]
                                                         
                    else:
                        syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del dbs_dict[mut]
                            del output_dbs_dict[mut]
                            
                #Third base of codon and first base of next adjacent codon    
                if read_index == 1:
                    try:
                    
                        if sample_seq[mut-2:mut] + gen_dbs_dict[mut][0] == sample_seq[mut-2:mut+1]:
                            pass
                            
                        if gen_dbs_dict[mut][1] + sample_seq[mut+2:mut+4] == sample_seq[mut+1:mut+4]:
                            pass
                        
                        #Check if first codon of original sequence is synonymous to mutated codon
                        if sample_seq[mut-2:mut] + gen_dbs_dict[mut][0] in syn_codon[sample_seq[mut-2:mut+1]]:
                            syn_mutation_codon_1 = random.choices([1,2],[syn_rate**0.5, 1-syn_rate**0.5])
                           
                        #Check if second codon of original sequence is synonymous to mutated codon
                        if gen_dbs_dict[mut][1] + sample_seq[mut+2:mut+4] in syn_codon[sample_seq[mut+1:mut+4]]:
                            syn_mutation_codon_2 = random.choices([1,2],[syn_rate**0.5, 1-syn_rate**0.5])
                            
                            if syn_mutation_codon_1[0] and syn_mutation_codon_2[0] == 1:
                                pass
                                
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                           
                        #Both mutated codons are non-synonymous to their original sequence codon                           
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                    except:
                        pass
                            
        #Reading frame starts at the third base
        if reading_frame == 3:
            
            #Apply syn/non-syn mutation rate for DBS
            for mut in list(gen_dbs_dict.keys()):
            
                read_index = mut%3
            
                #First + second base of codon
                if read_index == 0:
                    
                    if gen_dbs_dict[mut] + sample_seq[mut+2] == sample_seq[mut:mut+3]:
                        pass
                    
                    if gen_dbs_dict[mut] + sample_seq[mut+2] in syn_codon[sample_seq[mut:mut+3]]:
                        syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del dbs_dict[mut]
                            del output_dbs_dict[mut]
                                           
                    else:
                        syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del dbs_dict[mut]
                            del output_dbs_dict[mut]
                                                   
                #Second + third base of codon     
                if read_index == 1:
                    
                    if sample_seq[mut-1] + gen_dbs_dict[mut]  == sample_seq[mut-1:mut+2]:
                        pass
                    
                    if sample_seq[mut-1] + gen_dbs_dict[mut] in syn_codon[sample_seq[mut-1:mut+2]]:
                        syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del dbs_dict[mut]
                            del output_dbs_dict[mut]
                                                         
                    else:
                        syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                        
                        if syn_mutation[0] == 1:
                            pass
                        else:
                            del dbs_dict[mut]
                            del output_dbs_dict[mut]
                            
                #Third base of codon and first base of next adjacent codon    
                if read_index == 2:
                    try:
                        
                        if sample_seq[mut-2:mut] + gen_dbs_dict[mut][0] == sample_seq[mut-2:mut+1]:
                            pass
                            
                        if gen_dbs_dict[mut][1] + sample_seq[mut+2:mut+4] == sample_seq[mut+1:mut+4]:
                            pass
                        
                        #Check if first codon of original sequence is synonymous to mutated codon
                        if sample_seq[mut-2:mut] + gen_dbs_dict[mut][0] in syn_codon[sample_seq[mut-2:mut+1]]:
                            syn_mutation_codon_1 = random.choices([1,2],[syn_rate**0.5, 1-syn_rate**0.5])
                           
                        #Check if second codon of original sequence is synonymous to mutated codon
                        if gen_dbs_dict[mut][1] + sample_seq[mut+2:mut+4] in syn_codon[sample_seq[mut+1:mut+4]]:
                            syn_mutation_codon_2 = random.choices([1,2],[syn_rate**0.5, 1-syn_rate**0.5])
                            
                            if syn_mutation_codon_1[0] and syn_mutation_codon_2[0] == 1:
                                pass
                                
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                           
                        #Both mutated codons are non-synonymous to their original sequence codon                           
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                   
                    except:
                        pass
        
        #96 type SBS mutation frequency matrix
        sbs_mut_freq = sbs_sorted_df.copy().iloc[:,:1].reset_index()
        sbs_mut_freq['Frequency'] = 0
        for sbs_count in list(output_sbs_dict.values()):
            row = sbs_mut_freq[(sbs_mut_freq['SubType'] == sbs_count[1]) & (sbs_mut_freq['Type'] == sbs_count[0])].index
            sbs_mut_freq.iloc[row[0], 2] += 1
    
        #78 type DBS mutation frequency matrix
        dbs_mut_freq = dbs_sorted_df.copy().iloc[:,:1].reset_index()
        dbs_mut_freq['Frequency'] = 0
        for dbs_count in list(output_dbs_dict.values()):
            row = dbs_mut_freq[(dbs_mut_freq['Mutation Type'] == dbs_count)].index
            dbs_mut_freq.iloc[row[0], 2] += 1
            
        #12 type Single Base Insertion mutation frequency matrix
        ins_mut_freq = insertion_df.copy().iloc[:,:3].reset_index()
        ins_mut_freq['Frequency'] = 0
        for insertion_count in list(output_insertion_dict.values()):
            row = ins_mut_freq[(ins_mut_freq['Index'] == insertion_count)].index
            ins_mut_freq.iloc[row[0], 4] += 1
    
        #12 type Single Base Deletion mutation frequency matrix
        del_mut_freq = deletion_df.copy().iloc[:,:3].reset_index()
        del_mut_freq['Frequency'] = 0
        for deletion_count in list(output_deletion_dict.values()):
            row = del_mut_freq[(del_mut_freq['Index'] == deletion_count)].index
            del_mut_freq.iloc[row[0], 4] += 1
            
            
        #SBS Metadata
        sbs_mut_metadata = pd.DataFrame(index = range(len(output_sbs_dict)), columns = ["SubType", "Type", "Index"])
        for mutation in range(len(output_sbs_dict)):
            sbs_mut_metadata.loc[mutation, "SubType"] = list(output_sbs_dict.values())[mutation][1]
            sbs_mut_metadata.loc[mutation, "Type"] = list(output_sbs_dict.values())[mutation][0]
            sbs_mut_metadata.loc[mutation, "Index"] = list(output_sbs_dict.keys())[mutation]
            
        #DBS Metadata
        dbs_mut_metadata = pd.DataFrame(index = range(len(output_dbs_dict)), columns = ["2-mer Index", "Type", "Index"])
        for mutation in range(len(output_dbs_dict)):
            dbs_mut_metadata.loc[mutation, "2-mer Index"] = list(output_dbs_dict.values())[mutation][:2]
            dbs_mut_metadata.loc[mutation, "Type"] = list(output_dbs_dict.values())[mutation]
            dbs_mut_metadata.loc[mutation, "Index"] = list(output_dbs_dict.keys())[mutation]
            
        #Insertion Metadata
        ins_mut_metadata = pd.DataFrame(index = range(len(insertion_dict)), columns = ["Type", "Index"])
        for mutation in range(len(insertion_dict)):
            ins_mut_metadata.loc[mutation, "Type"] = list(insertion_dict.values())[mutation]
            ins_mut_metadata.loc[mutation, "Index"] = list(insertion_dict.keys())[mutation]
                    
        #Deletion Metadata
        del_mut_metadata = pd.DataFrame(index = range(len(deletion_dict)), columns = ["Type", "Index"])
        for mutation in range(len(deletion_dict)):
            del_mut_metadata.loc[mutation, "Type"] = list(deletion_dict.values())[mutation]
            del_mut_metadata.loc[mutation, "Index"] = list(deletion_dict.keys())[mutation]
            
            
        #Apply the mutations linearly 
        for index_sbs in list(sbs_dict.keys()):
            sample_seq = sample_seq[:index_sbs] + sbs_dict[index_sbs] + sample_seq[index_sbs+1:]
    
        for index_dbs in list(dbs_dict.keys()):
            sample_seq = sample_seq[:index_dbs] + dbs_dict[index_dbs] + sample_seq[index_dbs+2:]
    
        #Index of mutation adjusts each time an insertion or deletion is added linearly
        offset_index = 0
        
        for index_insertion in sorted(list(insertion_dict.keys())):
            sample_seq = sample_seq[:int(index_insertion) + offset_index] + insertion_dict[index_insertion] + sample_seq[int(index_insertion) + offset_index:]
            offset_index += 1
            
        for index_deletion in sorted(list(deletion_dict.keys())):
            sample_seq = sample_seq[:int(index_deletion) + offset_index] + sample_seq[int(index_deletion) + 1 + offset_index:]
            offset_index -= 1
      
        
      
        #Write mutated sequence to fasta file
        try:
            os.mkdir(seq_directory + "/" + cancer_type)
        except OSError:
            print ("Creation of the directory %s failed" % (seq_directory + "/" + cancer_type))
        else:
            print ("Successfully created the directory %s " % (seq_directory + "/" + cancer_type))
                
        #Write mutated sequence to fasta file
        sample_sequence_file = seq_directory + "/" + cancer_type + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages) + '.fasta'
        with open(sample_sequence_file, 'w+') as fasta_file:
            fasta_file.write(">" + str(cancer_type) + "_Lineage_" + str(number_of_lineages) + "_end_stage \n")
            fasta_file.write("")
            fasta_file.write(sample_seq)
            
        #Write SBS mutation frequency tables to csv file
        sbs_freq_path = mut_mat_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_sbs_freq_table.csv'
        with open(sbs_freq_path, 'w+'):
            sbs_mut_freq.to_csv(sbs_freq_path, index=False)
            
        #Write DBS mutation frequency tables to csv file
        dbs_freq_path = mut_mat_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_dbs_freq_table.csv'
        with open(dbs_freq_path, 'w+'):
            dbs_mut_freq.to_csv(dbs_freq_path, index=False)
          
        #Write Insertion mutation frequency tables to csv file
        insertion_freq_path = mut_mat_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_ins_freq_table.csv'
        with open(insertion_freq_path, 'w+'):
            ins_mut_freq.to_csv(insertion_freq_path, index=False)
            
        #Write Deletion mutation frequency tables to csv file
        deletion_freq_path = mut_mat_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_del_freq_table.csv'
        with open(deletion_freq_path, 'w+'):
            del_mut_freq.to_csv(deletion_freq_path, index=False)
            
        
        
        #Write SBS mutation index tables to csv file
        sbs_metadata_path = mut_metadata_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_sbs_index_table.csv'
        with open(sbs_metadata_path, 'w+'):
            sbs_mut_metadata.to_csv(sbs_metadata_path, index=False)
            
        #Write DBS mutation index tables to csv file
        dbs_metadata_path = mut_metadata_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_dbs_index_table.csv'
        with open(dbs_metadata_path, 'w+'):
            dbs_mut_metadata.to_csv(dbs_metadata_path, index=False)
          
        #Write Insertion mutation index tables to csv file
        insertion_metadata_path = mut_metadata_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_ins_index_table.csv'
        with open(insertion_metadata_path, 'w+'):
            ins_mut_metadata.to_csv(insertion_metadata_path, index=False)
            
        #Write Deletion mutation index tables to csv file
        deletion_metadata_path = mut_metadata_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_del_index_table.csv'
        with open(deletion_metadata_path, 'w+'):
            del_mut_metadata.to_csv(deletion_metadata_path, index=False)
 
        #SBS signatures simulated
        sbs_sigs_path = mut_sigs_directory + "/" + cancer_type + '_End_Stage_Lineage_sbs_sigs.txt'
        with open(sbs_sigs_path , 'a+') as outfile:
            outfile.write(str(sbs_sig_combination[0].tolist()))
            outfile.write("\n")
        #DBS signatures simulated
        dbs_sigs_path = mut_sigs_directory + "/" + cancer_type + '_End_Stage_Lineage_dbs_sigs.txt'
        with open(dbs_sigs_path , 'a+') as outfile:
            outfile.write(str(dbs_sig_combination[0].tolist()))
            outfile.write("\n")
        #ID signatures simulated
        id_sigs_path = mut_sigs_directory + "/" + cancer_type + '_End_Stage_Lineage_id_sigs.txt'
        with open(id_sigs_path , 'a+') as outfile:
            outfile.write(str(id_sig_combination[0].tolist()))
            outfile.write("\n")
            
    if simulation_type == "temporal":
        
        #Select combination of signatures and set up mutation probabilities for indels
        id_prop = sig_proportion(cancer_type, id_num_data, std_outlier)
        id_sig_combination = random.choices(list(id_prop.index), list(id_prop.values))
        #id_sig_combination = [x for x in id_sig_combination[0] if x in sig_weights[sig_weights['cancer_type'] == cancer_type]['signatures'].tolist()]
        
        #Select combination of signatures and set up mutation probabilities for single base substitution
        sbs_prop = sig_proportion(cancer_type, sbs_num_data, std_outlier)
        sbs_sig_combination = random.choices(list(sbs_prop.index), list(sbs_prop.values))  
        sbs_sig_combination = [x for x in sbs_sig_combination[0] if x in sig_weights[sig_weights['cancer_type'] == cancer_type]['signatures'].tolist()]
        
        #Select combination of signatures and set up mutation probabilities for double base substitution
        dbs_prop = sig_proportion(cancer_type, dbs_num_data, std_outlier)
        dbs_sig_combination = random.choices(list(dbs_prop.index), list(dbs_prop.values))  
        #dbs_sig_combination = [x for x in dbs_sig_combination[0] if x in sig_weights[sig_weights['cancer_type'] == cancer_type]['signatures'].tolist()]
         
        for i in ["Early", "Late"]:
            
            #Index of each kmer
            sample_index_dict = defaultdict(list)
        
            for seq_len in range(0, len(sample_seq)-6):
                if 'N' in sample_seq[seq_len:seq_len+6]:
                    pass
        
                else:
                    sample_index_dict[str(sample_seq[seq_len:seq_len+6])].append(seq_len)
            
            #Count of each kmer
            sample_count_dict = {key: len(value) for key, value in sample_index_dict.items()}
           
            #Normalization of mutation probabilities to whole genome burden
            
            ref_dir = abs_path("6-mer",  "Directory")
            kmer_ref = (glob.glob(ref_dir+ "/6-mer_chr*"))
            kmer_count = pd.read_csv(kmer_ref[0], index_col=0)['count'].fillna(0)
            for i in kmer_ref[1:-1]:
                sample = pd.read_csv(i, index_col=0)['count'].fillna(0)
                kmer_count = kmer_count.add(sample, fill_value=0)
              
                
            kmer_reference_count_dict = dict(zip(pd.read_csv(kmer_ref[0], index_col=0)["6"], kmer_count))
            
            sample_expected_count_dict = dict(zip(pd.read_csv(kmer_ref[0], index_col=0)["6"], kmer_count))
            for key in sample_expected_count_dict:
                sample_expected_count_dict[key] = sample_expected_count_dict[key] * (len(sample_seq) / sum(list(kmer_reference_count_dict.values()))) 
            
            normalize_constant = {k: float(sample_count_dict[k])/sample_expected_count_dict[k] for k in sample_count_dict}
            normalized_sample_count_dict = {k: int(sample_count_dict[k]/normalize_constant[k]) for k in sample_count_dict}

       
            print('Simulating Tumour Stage ' + str(i)+ ' of Lineage ' + str(number_of_lineages))
        
            #Types of mutations for a specific generation  of a specific lineage
            gen_insertion_dict={}
            gen_deletion_dict={}
            gen_sbs_dict={}
            gen_dbs_dict={}  
            
            #Types of mutations for a specific lineage and generation to be outputted as a frequency matrix
            output_insertion_dict={}
            output_deletion_dict={}
            output_sbs_dict={}
            output_dbs_dict={}
        
            id_prob = mut_prob(cancer_type, id_prop_data, id_num_data, id_freq_folder_path, "ID", std_outlier)
            id_mut_sig_df = pd.concat([id_prob.loc[:47,'Mutation Type':'ID_type'], id_prob[id_sig_combination[0]].dropna()], axis=1)
                
            id_total_prob_df = id_mut_sig_df.iloc[:,:5]
            id_total_prob_df['total_prob']= id_mut_sig_df.loc[:, 'ID_type':].sum(axis=1)
            id_total_prob_df['total_prob'] = id_total_prob_df['total_prob'].multiply(power)
               
            id_sorted_df = id_total_prob_df.sort_values(['ID_type','Index','Size'])
            id_sorted_df.set_index(['Index'], inplace=True)
                
            deletion_df = id_sorted_df.iloc[:12, :]
            insertion_df = id_sorted_df.iloc[12:24, :]
                
            print("Insertion and Deletion mutation probability matrix set up.")
            print(id_sig_combination)
            
           
            sbs_prob = mut_prob(cancer_type, 
                               sbs_prop_data,
                               sbs_num_data, 
                               sbs_freq_folder_path,
                               "SBS",
                               std_outlier)
            
            #Constant multiplier for SBS signature activity that differs between early and late stages
            weights = sig_weights[sig_weights['cancer_type'] == cancer_type].set_index("signatures")
            for j in sbs_sig_combination:
                sbs_prob[j] = sbs_prob[j]*weights.loc[j,i.lower()]
            
            sbs_mut_sig_df = pd.concat([sbs_prob.loc[:,'Type':'SubType'], sbs_prob[sbs_sig_combination[0]]], axis=1)
            
            sbs_total_prob_df = sbs_mut_sig_df.iloc[:,:2]
            sbs_total_prob_df['total_prob'] = sbs_mut_sig_df.iloc[:, 2:].sum(axis=1)
                
            sbs_total_prob_df['total_prob'] = sbs_total_prob_df[sbs_total_prob_df.select_dtypes(include=['number']).columns] * power
            
            sbs_sorted_df = sbs_total_prob_df.sort_values(['SubType'])
            sbs_sorted_df.set_index(['SubType'], inplace=True)
            
            print("Single Base Substitution mutation probability matrix set up.")
            print(sbs_sig_combination)
                
            dbs_prob = mut_prob(cancer_type, 
                               dbs_prop_data,
                               dbs_num_data, 
                               dbs_freq_folder_path,
                               "DBS",
                               std_outlier)
               
            dbs_mut_sig_df = pd.concat([dbs_prob.loc[:,'Mutation Type':'2mer_index'], dbs_prob[dbs_sig_combination[0]]], axis=1)
            
            dbs_total_prob_df = dbs_mut_sig_df.copy()    
            dbs_total_prob_df['total_prob']= dbs_mut_sig_df.iloc[:, 1:].sum(axis=1) 
            dbs_total_prob_df['total_prob'] = dbs_total_prob_df[dbs_total_prob_df.select_dtypes(include=['number']).columns] * power
            
            dbs_sorted_df =  dbs_total_prob_df.sort_values(['Mutation Type'])
            dbs_sorted_df.set_index(['2mer_index'], inplace=True)   
            
            print("Double Base Substitution mutation probability matrix set up.")
            print(dbs_sig_combination)
        
            for keys in list(sample_index_dict.keys()):
                
                count=1
                for length in list(range(1, len(keys))):
                    if keys[length] == keys[0]:
                        count += 1
                    else:
                        break
          
                #Insertion
                if "G" in keys[:count] or "A" in keys[:count]:
                    pass
                
                else:
                    temp_i_prob = [insertion_df.loc[insertion_df['Size'] == count].loc[keys[:count],'total_prob']] 
                    insertion_count = np.random.binomial(normalized_sample_count_dict[keys], temp_i_prob, 1) 
                
                    for m in range(sum(i_number > 0 for i_number in list(insertion_count))):
                        insertion_index = np.random.choice(sample_index_dict[keys])
                        insertion_dict[insertion_index] = keys[0]
                        gen_insertion_dict[insertion_index] = keys[0]
                        output_insertion_dict[insertion_index] = keys[:count]
                    
                    
                #Deletion
                if "G" in keys[:count] or "A" in keys[:count]:
                    pass
                
                else:
                    temp_d_prob = [deletion_df.loc[deletion_df['Size'] == count].loc[keys[:count],'total_prob']] 
                    deletion_count = np.random.binomial(normalized_sample_count_dict[keys], temp_d_prob, 1)  
                
                    for n in range(sum(d_number > 0 for d_number in list(deletion_count))):
                        deletion_index = np.random.choice(sample_index_dict[keys])
                        deletion_dict[deletion_index] = keys[0]
                        gen_deletion_dict[deletion_index] = keys[0]
                        output_deletion_dict[deletion_index] = keys[:count]
            
                #SBS
                if "G" in keys[1] or "A" in keys[1]:
                    pass
                
                else:
                    sbs_types = list(sbs_sorted_df.loc[keys[:3],'Type']) + [None]
                    sbs_prob = list(sbs_sorted_df.loc[keys[:3],'total_prob'])+ [1- sum(list(sbs_sorted_df.loc[keys[:3],'total_prob']))]             
                    sbs = [sbs_num for sbs_num in list(np.random.choice(sbs_types,p=sbs_prob, size = normalized_sample_count_dict[keys])) if sbs_num]
                    for sbs_value in sbs:
                        single_base_index = np.random.choice(sample_index_dict[keys])
                        sbs_dict[single_base_index] = sbs_value[2]
                        gen_sbs_dict[single_base_index] = sbs_value[2]
                        output_sbs_dict[single_base_index] = [sbs_value, sample_seq[single_base_index:single_base_index+3]]
                             
                #DBS   
                if keys[:2] in dbs_sorted_df.index:
                    dbs_types = list(dbs_sorted_df.loc[keys[:2],'Mutation Type']) + [None]     
                    dbs_prob = list(dbs_sorted_df.loc[keys[:2],'total_prob'])+ [1- sum(list(dbs_sorted_df.loc[keys[:2],'total_prob']))]
                    dbs = [dbs_num for dbs_num in list(np.random.choice(dbs_types,p=dbs_prob, size = normalized_sample_count_dict[keys])) if dbs_num]
                    for dbs_value in dbs:
                        double_base_index = np.random.choice(sample_index_dict[keys])
                        dbs_dict[double_base_index] = dbs_value[3:]
                        gen_dbs_dict[double_base_index] = dbs_value[3:]
                        output_dbs_dict[double_base_index] = dbs_value
                        
                else:
                    pass
                
            #SBS Synonymous/Non-Synonymous
            #Reading frame starts at the first base
            if reading_frame == 1:
                
                #Apply syn/non-syn mutation rate for SBS
                for mut in list(gen_sbs_dict.keys()):
                
                    read_index = mut%3
                
                    #First base of codon
                    if read_index == 1:
                        
                        if gen_sbs_dict[mut] + sample_seq[mut+1:mut+3] == sample_seq[mut:mut+3]:
                            pass
                        
                        if gen_sbs_dict[mut] + sample_seq[mut+1:mut+3] in syn_codon[sample_seq[mut:mut+3]]:
                            syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                        
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                                       
                    #Second base of codon     
                    if read_index == 2:
                        
                        if sample_seq[mut-1] + gen_sbs_dict[mut] + sample_seq[mut+1] == sample_seq[mut-1:mut+2]:
                            pass
                        
                        if sample_seq[mut-1] + gen_sbs_dict[mut] + sample_seq[mut+1] in syn_codon[sample_seq[mut-1:mut+2]]:
                            syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                                             
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                
                    #Third base of codon     
                    if read_index == 0:
                        
                        if sample_seq[mut-2:mut] + gen_sbs_dict[mut] == sample_seq[mut-2:mut+1]:
                            pass
                        
                        if sample_seq[mut-2:mut] + gen_sbs_dict[mut] in syn_codon[sample_seq[mut-2:mut+1]]:
                            syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                               
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                
            #Reading frame starts at the second base
            if reading_frame == 2:
                
                #Apply syn/non-syn mutation rate for SBS
                for mut in list(gen_sbs_dict.keys()):
                
                    read_index = mut%3
                
                    #First base of codon
                    if read_index == 2:
                        
                        if gen_sbs_dict[mut] + sample_seq[mut+1:mut+3] == sample_seq[mut:mut+3]:
                            pass
                        
                        if gen_sbs_dict[mut] + sample_seq[mut+1:mut+3] in syn_codon[sample_seq[mut:mut+3]]:
                            syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                               
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                                       
                        #Second base of codon     
                    if read_index == 0:
                        
                        if sample_seq[mut-1] + gen_sbs_dict[mut] + sample_seq[mut+1] == sample_seq[mut-1:mut+2]:
                            pass
                        
                        if sample_seq[mut-1] + gen_sbs_dict[mut] + sample_seq[mut+1] in syn_codon[sample_seq[mut-1:mut+2]]:
                            syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                                             
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                
                    #Third base of codon     
                    if read_index == 1:
                        
                        if sample_seq[mut-2:mut] + gen_sbs_dict[mut] == sample_seq[mut-2:mut+1]:
                            pass
                        
                        if sample_seq[mut-2:mut] + gen_sbs_dict[mut] in syn_codon[sample_seq[mut-2:mut+1]]:
                            syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                               
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                
            #Reading frame starts at the third base
            if reading_frame == 3:
                
                #Apply syn/non-syn mutation rate for SBS
                for mut in list(gen_sbs_dict.keys()):
                
                    read_index = mut%3
                
                    #First base of codon
                    if read_index == 0:
                        
                        if gen_sbs_dict[mut] + sample_seq[mut+1:mut+3] == sample_seq[mut:mut+3]:
                            pass
                        
                        if gen_sbs_dict[mut] + sample_seq[mut+1:mut+3] in syn_codon[sample_seq[mut:mut+3]]:
                            syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                               
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                                       
                    #Second base of codon     
                    if read_index == 1:
                        
                        if sample_seq[mut-1] + gen_sbs_dict[mut] + sample_seq[mut+1] == sample_seq[mut-1:mut+2]:
                            pass
                        
                        if sample_seq[mut-1] + gen_sbs_dict[mut] + sample_seq[mut+1] in syn_codon[sample_seq[mut-1:mut+2]]:
                            syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                                             
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 2:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                
                    #Third base of codon     
                    if read_index == 0:
                        
                        if sample_seq[mut-2:mut] + gen_sbs_dict[mut] == sample_seq[mut-2:mut+1]:
                            pass
                        
                        if sample_seq[mut-2:mut] + gen_sbs_dict[mut] in syn_codon[sample_seq[mut-2:mut+1]]:
                            syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                               
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del sbs_dict[mut]
                                del output_sbs_dict[mut]
                                
            #DBS Synonymous/Non-Synonymous
            
            #Reading frame starts at the first base
            if reading_frame == 1:
                
                #Apply syn/non-syn mutation rate for DBS
                for mut in list(gen_dbs_dict.keys()):
                
                    read_index = mut%3
                
                    #First + second base of codon
                    if read_index == 1:
                        
                        if gen_dbs_dict[mut] + sample_seq[mut+2] == sample_seq[mut:mut+3]:
                            pass
                        
                        if gen_dbs_dict[mut] + sample_seq[mut+2] in syn_codon[sample_seq[mut:mut+3]]:
                            syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                                               
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                                                       
                    #Second + third base of codon     
                    if read_index == 2:
                        
                        if sample_seq[mut-1] + gen_dbs_dict[mut]  == sample_seq[mut-1:mut+2]:
                            pass
                        
                        if sample_seq[mut-1] + gen_dbs_dict[mut] in syn_codon[sample_seq[mut-1:mut+2]]:
                            syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                                                             
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                                
                    #Third base of codon and first base of next adjacent codon    
                    if read_index == 0:
                        try:
                        
                            if sample_seq[mut-2:mut] + gen_dbs_dict[mut][0] == sample_seq[mut-2:mut+1]:
                                pass
                                
                            if gen_dbs_dict[mut][1] + sample_seq[mut+2:mut+4] == sample_seq[mut+1:mut+4]:
                                pass
                            
                            #Check if first codon of original sequence is synonymous to mutated codon
                            if sample_seq[mut-2:mut] + gen_dbs_dict[mut][0] in syn_codon[sample_seq[mut-2:mut+1]]:
                                syn_mutation_codon_1 = random.choices([1,2],[syn_rate**0.5, 1-syn_rate**0.5])
                               
                            #Check if second codon of original sequence is synonymous to mutated codon
                            if gen_dbs_dict[mut][1] + sample_seq[mut+2:mut+4] in syn_codon[sample_seq[mut+1:mut+4]]:
                                syn_mutation_codon_2 = random.choices([1,2],[syn_rate**0.5, 1-syn_rate**0.5])
                                
                                if syn_mutation_codon_1[0] and syn_mutation_codon_2[0] == 1:
                                    pass
                                    
                                else:
                                    del dbs_dict[mut]
                                    del output_dbs_dict[mut]
                               
                            #Both mutated codons are non-synonymous to their original sequence codon                           
                            else:
                                syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                                
                                if syn_mutation[0] == 1:
                                    pass
                                else:
                                    del dbs_dict[mut]
                                    del output_dbs_dict[mut]
                        
                        except:
                            pass
            
            #Reading frame starts at the second base
            if reading_frame == 2:
                
                #Apply syn/non-syn mutation rate for DBS
                for mut in list(gen_dbs_dict.keys()):
                
                    read_index = mut%3
                
                    #First + second base of codon
                    if read_index == 2:
                        
                        if gen_dbs_dict[mut] + sample_seq[mut+2] == sample_seq[mut:mut+3]:
                            pass
                        
                        if gen_dbs_dict[mut] + sample_seq[mut+2] in syn_codon[sample_seq[mut:mut+3]]:
                            syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                                               
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                                                       
                    #Second + third base of codon     
                    if read_index == 0:
                        
                        if sample_seq[mut-1] + gen_dbs_dict[mut]  == sample_seq[mut-1:mut+2]:
                            pass
                        
                        if sample_seq[mut-1] + gen_dbs_dict[mut] in syn_codon[sample_seq[mut-1:mut+2]]:
                            syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                                                             
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                                
                    #Third base of codon and first base of next adjacent codon    
                    if read_index == 1:
                        
                        try:
                        
                            if sample_seq[mut-2:mut] + gen_dbs_dict[mut][0] == sample_seq[mut-2:mut+1]:
                                pass
                                
                            if gen_dbs_dict[mut][1] + sample_seq[mut+2:mut+4] == sample_seq[mut+1:mut+4]:
                                pass
                            
                            #Check if first codon of original sequence is synonymous to mutated codon
                            if sample_seq[mut-2:mut] + gen_dbs_dict[mut][0] in syn_codon[sample_seq[mut-2:mut+1]]:
                                syn_mutation_codon_1 = random.choices([1,2],[syn_rate**0.5, 1-syn_rate**0.5])
                               
                            #Check if second codon of original sequence is synonymous to mutated codon
                            if gen_dbs_dict[mut][1] + sample_seq[mut+2:mut+4] in syn_codon[sample_seq[mut+1:mut+4]]:
                                syn_mutation_codon_2 = random.choices([1,2],[syn_rate**0.5, 1-syn_rate**0.5])
                                
                                if syn_mutation_codon_1[0] and syn_mutation_codon_2[0] == 1:
                                    pass
                                    
                                else:
                                    del dbs_dict[mut]
                                    del output_dbs_dict[mut]
                               
                            #Both mutated codons are non-synonymous to their original sequence codon                           
                            else:
                                syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                                
                                if syn_mutation[0] == 1:
                                    pass
                                else:
                                    del dbs_dict[mut]
                                    del output_dbs_dict[mut]
                        except:
                            pass
                                
            #Reading frame starts at the third base
            if reading_frame == 3:
                
                #Apply syn/non-syn mutation rate for DBS
                for mut in list(gen_dbs_dict.keys()):
                
                    read_index = mut%3
                
                    #First + second base of codon
                    if read_index == 0:
                        
                        if gen_dbs_dict[mut] + sample_seq[mut+2] == sample_seq[mut:mut+3]:
                            pass
                        
                        if gen_dbs_dict[mut] + sample_seq[mut+2] in syn_codon[sample_seq[mut:mut+3]]:
                            syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                                               
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                                                       
                    #Second + third base of codon     
                    if read_index == 1:
                        
                        if sample_seq[mut-1] + gen_dbs_dict[mut]  == sample_seq[mut-1:mut+2]:
                            pass
                        
                        if sample_seq[mut-1] + gen_dbs_dict[mut] in syn_codon[sample_seq[mut-1:mut+2]]:
                            syn_mutation = random.choices([1,2],[syn_rate, 1-syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                                                             
                        else:
                            syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                            
                            if syn_mutation[0] == 1:
                                pass
                            else:
                                del dbs_dict[mut]
                                del output_dbs_dict[mut]
                                
                    #Third base of codon and first base of next adjacent codon    
                    if read_index == 2:
                        try:
                            if sample_seq[mut-2:mut] + gen_dbs_dict[mut][0] == sample_seq[mut-2:mut+1]:
                                pass
                                
                            if gen_dbs_dict[mut][1] + sample_seq[mut+2:mut+4] == sample_seq[mut+1:mut+4]:
                                pass
                            
                            #Check if first codon of original sequence is synonymous to mutated codon
                            if sample_seq[mut-2:mut] + gen_dbs_dict[mut][0] in syn_codon[sample_seq[mut-2:mut+1]]:
                                syn_mutation_codon_1 = random.choices([1,2],[syn_rate**0.5, 1-syn_rate**0.5])
                               
                            #Check if second codon of original sequence is synonymous to mutated codon
                            if gen_dbs_dict[mut][1] + sample_seq[mut+2:mut+4] in syn_codon[sample_seq[mut+1:mut+4]]:
                                syn_mutation_codon_2 = random.choices([1,2],[syn_rate**0.5, 1-syn_rate**0.5])
                                
                                if syn_mutation_codon_1[0] and syn_mutation_codon_2[0] == 1:
                                    pass
                                    
                                else:
                                    del dbs_dict[mut]
                                    del output_dbs_dict[mut]
                               
                            #Both mutated codons are non-synonymous to their original sequence codon                           
                            else:
                                syn_mutation = random.choices([1,2],[non_syn_rate, 1-non_syn_rate])
                                
                                if syn_mutation[0] == 1:
                                    pass
                                else:
                                    del dbs_dict[mut]
                                    del output_dbs_dict[mut]
                        except:
                            pass
                   
            
            #96 type SBS mutation frequency matrix
            sbs_mut_freq = sbs_sorted_df.copy().iloc[:,:1].reset_index()
            sbs_mut_freq['Frequency'] = 0
            for sbs_count in list(output_sbs_dict.values()):
                row = sbs_mut_freq[(sbs_mut_freq['SubType'] == sbs_count[1]) & (sbs_mut_freq['Type'] == sbs_count[0])].index
                sbs_mut_freq.iloc[row[0], 2] += 1
        
            #78 type DBS mutation frequency matrix
            dbs_mut_freq = dbs_sorted_df.copy().iloc[:,:1].reset_index()
            dbs_mut_freq['Frequency'] = 0
            for dbs_count in list(output_dbs_dict.values()):
                row = dbs_mut_freq[(dbs_mut_freq['Mutation Type'] == dbs_count)].index
                dbs_mut_freq.iloc[row[0], 2] += 1
                
            #12 type Single Base Insertion mutation frequency matrix
            ins_mut_freq = insertion_df.copy().iloc[:,:3].reset_index()
            ins_mut_freq['Frequency'] = 0
            for insertion_count in list(output_insertion_dict.values()):
                row = ins_mut_freq[(ins_mut_freq['Index'] == insertion_count)].index
                ins_mut_freq.iloc[row[0], 4] += 1
          
            #12 type Single Base Deletion mutation frequency matrix
            del_mut_freq = deletion_df.copy().iloc[:,:3].reset_index()
            del_mut_freq['Frequency'] = 0
            for deletion_count in list(output_deletion_dict.values()):
                row = del_mut_freq[(del_mut_freq['Index'] == deletion_count)].index
                del_mut_freq.iloc[row[0], 4] += 1
                
            #SBS Metadata
            sbs_mut_metadata = pd.DataFrame(index = range(len(output_sbs_dict)), columns = ["SubType", "Type", "Index"])
            for mutation in range(len(output_sbs_dict)):
                sbs_mut_metadata.loc[mutation, "SubType"] = list(output_sbs_dict.values())[mutation][1]
                sbs_mut_metadata.loc[mutation, "Type"] = list(output_sbs_dict.values())[mutation][0]
                sbs_mut_metadata.loc[mutation, "Index"] = list(output_sbs_dict.keys())[mutation]
            
            #DBS Metadata
            dbs_mut_metadata = pd.DataFrame(index = range(len(output_dbs_dict)), columns = ["2-mer Index", "Type", "Index"])
            for mutation in range(len(output_dbs_dict)):
                dbs_mut_metadata.loc[mutation, "2-mer Index"] = list(output_dbs_dict.values())[mutation][:2]
                dbs_mut_metadata.loc[mutation, "Type"] = list(output_dbs_dict.values())[mutation]
                dbs_mut_metadata.loc[mutation, "Index"] = list(output_dbs_dict.keys())[mutation]
            
            #Insertion Metadata
            ins_mut_metadata = pd.DataFrame(index = range(len(insertion_dict)), columns = ["Type", "Index"])
            for mutation in range(len(insertion_dict)):
                ins_mut_metadata.loc[mutation, "Type"] = list(insertion_dict.values())[mutation]
                ins_mut_metadata.loc[mutation, "Index"] = list(insertion_dict.keys())[mutation]
                    
            #Deletion Metadata
            del_mut_metadata = pd.DataFrame(index = range(len(deletion_dict)), columns = ["Type", "Index"])
            for mutation in range(len(deletion_dict)):
                del_mut_metadata.loc[mutation, "Type"] = list(deletion_dict.values())[mutation]
                del_mut_metadata.loc[mutation, "Index"] = list(deletion_dict.keys())[mutation]
                
                
            #Apply the mutations linearly 
            for index_sbs in list(sbs_dict.keys()):
                sample_seq = sample_seq[:index_sbs] + sbs_dict[index_sbs] + sample_seq[index_sbs+1:]
        
            for index_dbs in list(dbs_dict.keys()):
                sample_seq = sample_seq[:index_dbs] + dbs_dict[index_dbs] + sample_seq[index_dbs+2:]
        
            #Index of mutation adjusts each time an insertion or deletion is added linearly
            offset_index = 0
            
            for index_insertion in sorted(list(insertion_dict.keys())):
                sample_seq = sample_seq[:int(index_insertion) + offset_index] + insertion_dict[index_insertion] + sample_seq[int(index_insertion) + offset_index:]
                offset_index += 1
                
            for index_deletion in sorted(list(deletion_dict.keys())):
                sample_seq = sample_seq[:int(index_deletion) + offset_index] + sample_seq[int(index_deletion) + 1 + offset_index:]
                offset_index -= 1
          
            #Write mutated sequence to fasta file
            try:
                os.mkdir(seq_directory + "/" + cancer_type)
            except OSError:
                print ("Creation of the directory %s failed" % (seq_directory + "/" + cancer_type))
            else:
                print ("Successfully created the directory %s " % (seq_directory + "/" + cancer_type))
                    
            #Write mutated sequence to fasta file
            sample_sequence_file = seq_directory + "/" + cancer_type + "/" + cancer_type + '_' + str(i) + '_Stage_Lineage_' + str(number_of_lineages) + '.fasta'
            with open(sample_sequence_file, 'w+') as fasta_file:
                fasta_file.write(">" + str(cancer_type) + "_Lineage_" + str(number_of_lineages) + str(i) + "_stage \n")
                fasta_file.write("")
                fasta_file.write(sample_seq)
                
                
                
            #Write SBS mutation frequency tables to csv file
            sbs_freq_path = mut_mat_directory + "/" + cancer_type + '_'+ str(i) + '_Stage_Lineage_' + str(number_of_lineages)+ '_sbs_freq_table.csv'
            with open(sbs_freq_path, 'w+'):
                sbs_mut_freq.to_csv(sbs_freq_path, index=False)
                
            #Write DBS mutation frequency tables to csv file
            dbs_freq_path = mut_mat_directory + "/" + cancer_type + '_'+ str(i) + '_Stage_Lineage_' + str(number_of_lineages)+ '_dbs_freq_table.csv'
            with open(dbs_freq_path, 'w+'):
                dbs_mut_freq.to_csv(dbs_freq_path, index=False)
              
            #Write Insertion mutation frequency tables to csv file
            insertion_freq_path = mut_mat_directory + "/" + cancer_type + '_'+ str(i) + '_Stage_Lineage_' + str(i) + '_Lineage_' + str(number_of_lineages)+ '_ins_freq_table.csv'
            with open(insertion_freq_path, 'w+'):
                ins_mut_freq.to_csv(insertion_freq_path, index=False)
                
            #Write Deletion mutation frequency tables to csv file
            deletion_freq_path = mut_mat_directory + "/" + cancer_type + '_'+ str(i) + '_Stage_Lineage_' + str(number_of_lineages)+ '_del_freq_table.csv'
            with open(deletion_freq_path, 'w+'):
                del_mut_freq.to_csv(deletion_freq_path, index=False)
                
                
                
            #Write SBS mutation index tables to csv file
            sbs_metadata_path = mut_metadata_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_sbs_index_table.csv'
            with open(sbs_metadata_path, 'w+'):
                sbs_mut_metadata.to_csv(sbs_metadata_path, index=False)
                
            #Write DBS mutation index tables to csv file
            dbs_metadata_path = mut_metadata_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_dbs_index_table.csv'
            with open(dbs_metadata_path, 'w+'):
                dbs_mut_metadata.to_csv(dbs_metadata_path, index=False)
              
            #Write Insertion mutation index tables to csv file
            insertion_metadata_path = mut_metadata_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_ins_index_table.csv'
            with open(insertion_metadata_path, 'w+'):
                ins_mut_metadata.to_csv(insertion_metadata_path, index=False)
                
            #Write Deletion mutation index tables to csv file
            deletion_metadata_path = mut_metadata_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_del_index_table.csv'
            with open(deletion_metadata_path, 'w+'):
                del_mut_metadata.to_csv(deletion_metadata_path, index=False)
                
   
            #SBS signatures simulated
            sbs_sigs_path = mut_sigs_directory + "/" + cancer_type + "_" + str(i) + '_Stage_Lineage_sbs_sigs.txt'
            with open(sbs_sigs_path , 'a+') as outfile:
                outfile.write(str(sbs_sig_combination[0].tolist()))
                outfile.write("\n")
            #DBS signatures simulated
            dbs_sigs_path = mut_sigs_directory + "/" + cancer_type + "_" + str(i) + '_Stage_Lineage_dbs_sigs.txt'
            with open(dbs_sigs_path , 'a+') as outfile:
                outfile.write(str(dbs_sig_combination[0].tolist()))
                outfile.write("\n")
            #ID signatures simulated
            id_sigs_path = mut_sigs_directory + "/" + cancer_type + "_" + str(i) + '_Stage_Lineage_id_sigs.txt'
            with open(id_sigs_path , 'a+') as outfile:
                outfile.write(str(id_sig_combination[0].tolist()))
                outfile.write("\n")
                
        
#Cell 6 Loaded
print("Cell 6 (SomaticSiMu) of 6 Loaded")

# %%

def main(): 
    
    # Initiate the parser
    parser = argparse.ArgumentParser(description="SomaticSiMu Parameters")
    
    # Add long and short argument
    parser.add_argument("--generation", "-g", help="number of simulated sequences", default=10)
    parser.add_argument("--cancer", "-c", help="cancer type")
    parser.add_argument("--reading_frame", "-f", help="index start of reading frame", default=1)
    parser.add_argument("--std", "-s", help="exclude signature data outside of n std from the mean", default=3)
    parser.add_argument("--simulation_type", "-v", help="simulation type", default="end")
    parser.add_argument("--slice_start", "-a", help="start of the slice of the input sequence")
    parser.add_argument("--slice_end", "-b", help="end of the slice of the input sequence")
    parser.add_argument("--power", "-p", help="multiplier of mutation burden from burden observed in in vivo samples", default=1)
    parser.add_argument("--syn_rate", "-x", help="proportion of synonymous mutations out of all simulated mutations kept in the output simulated sequence", default=1)
    parser.add_argument("--non_syn_rate", "-y", help="proportion of non-synonymous mutations out of all simulated mutations kept in the output simulated sequence", default=1)
    parser.add_argument("--reference", "-r", help="full file path of reference sequence used as input for the simulation")
    # Read arguments from the command line
    args = parser.parse_args()
  
   
        
    if str(args.cancer) == "all":
        
        for item in cancer_type_list:
            try:
            
                iterable = range(1, int(args.generation) + 1)
                pool = multiprocessing.Pool(multiprocessing.cpu_count())
                starttime= time.time()
                
                cancer_type = item
                reading_frame = int(args.reading_frame)
                std_outlier = args.std
                simulation_type = str(args.simulation_type)
                sequence_abs_path = str(args.reference)
                slice_start = int(args.slice_start)
                slice_end = int(args.slice_end)
                power = args.power
                syn_rate = args.syn_rate
                non_syn_rate = args.non_syn_rate
                            
                func = partial(somatic_sim, cancer_type, reading_frame, std_outlier, simulation_type, sequence_abs_path, slice_start, slice_end, power, syn_rate, non_syn_rate)
            
                pool.map(func , iterable )
                pool.close()
                pool.join()
                print('That took {} seconds'.format(time.time() - starttime))
            
            except:
                pass
    else:
        try:
            if args.cancer not in cancer_type_list:
                print('Cancer type not found in existing types.')
            else:
                iterable = range(1, int(args.generation) + 1)
                pool = multiprocessing.Pool(multiprocessing.cpu_count())
                starttime= time.time()
                
                cancer_type = str(args.cancer)
                reading_frame = int(args.reading_frame)
                std_outlier = args.std
                simulation_type = str(args.simulation_type)
                sequence_abs_path = str(args.reference)
                slice_start = int(args.slice_start)
                slice_end = int(args.slice_end)
                power = args.power
                syn_rate = args.syn_rate
                non_syn_rate = args.non_syn_rate
                            
                func = partial(somatic_sim, cancer_type, reading_frame, std_outlier, simulation_type, sequence_abs_path, slice_start, slice_end, power, syn_rate, non_syn_rate)
            
                pool.map(func , iterable )
                pool.close()
                pool.join()
                print('That took {} seconds'.format(time.time() - starttime))
        except:
            pass
        
    
        
if __name__ ==  '__main__':
    main()


#%%


#for i in [80, 82, 83, 84, 94, 95]:
#    try:
#        somatic_sim(cancer_type="Bone-Benign", reading_frame=1, std_outlier=3, number_of_lineages=i, simulation_type="end", sequence_abs_path=input_file_path, slice_start=1, slice_end=50818467,power=1, syn_rate=1, non_syn_rate=1)
#    except:
#        print(i)
        
#for i in [20, 26, 27, 28, 29, 30]:
#    try:
#        somatic_sim(cancer_type="Myeloid-MPN", reading_frame=1, std_outlier=3, number_of_lineages=i, simulation_type="end", sequence_abs_path=input_file_path, slice_start=1, slice_end=50818467,power=1, syn_rate=1, non_syn_rate=1)
#    except:
#        print(i)
        
#for i in [74, 93, 95,  97]:
#    try:
#        somatic_sim(cancer_type="CNS-PiloAstro", reading_frame=1, std_outlier=3, number_of_lineages=i, simulation_type="end", sequence_abs_path=input_file_path, slice_start=1, slice_end=50818467,power=1, syn_rate=1, non_syn_rate=1)
#    except:
#        print("ERROR" + i)







