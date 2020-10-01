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
print("Cell 1 of 8 Loaded")

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
print("Cell 2 of 8 Loaded")


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
print("Cell 3 of 8 Loaded")

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
print("Cell 4 of 8 Loaded")

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
print("Cell 5 of 8 Loaded")

#%%
def somatic_sim(cancer_type, reading_frame, std_outlier, number_of_lineages, simulation_type, sequence_abs_path, slice_start, slice_end,power=1, syn_rate=1, non_syn_rate=1):
 
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
        kmer_ref = (glob.glob("/Users/davidchen/Documents/GitHub/SomaticSiMu/kmer_ref_count/6-mer/6-mer_chr*"))
        kmer_count = pd.read_csv(kmer_ref[0], index_col=0)['count'].fillna(0)
        for i in kmer_ref[1:-1]:
            sample = pd.read_csv(i, index_col=0)['count'].fillna(0)
            kmer_count = kmer_count + sample
            
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
        for mutation in range(len(output_sbs_dict)):
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
        for mutation in range(len(deletions_dict)):
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
        sample_sequence_file = seq_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages) + '.fasta'
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
            
        
        
        #Write SBS mutation frequency tables to csv file
        sbs_metadata_path = mut_metadata_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_sbs_index_table.csv'
        with open(sbs_freq_path, 'w+'):
            sbs_mut_metadata.to_csv(sbs_metadata_path, index=False)
            
        #Write DBS mutation frequency tables to csv file
        dbs_metadata_path = mut_metadata_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_dbs_index_table.csv'
        with open(dbs_freq_path, 'w+'):
            dbs_mut_metadata.to_csv(dbs_metadata_path, index=False)
          
        #Write Insertion mutation frequency tables to csv file
        insertion_metadata_path = mut_metadata_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_ins_index_table.csv'
        with open(insertion_freq_path, 'w+'):
            ins_mut_metadata.to_csv(insertion_metadata_path, index=False)
            
        #Write Deletion mutation frequency tables to csv file
        deletion_metadata_path = mut_metadata_directory + "/" + cancer_type + '_End_Stage_Lineage_' + str(number_of_lineages)+ '_del_index_table.csv'
        with open(deletion_freq_path, 'w+'):
            del_mut_metadata.to_csv(deletion_metadata_path, index=False)
    
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
            kmer_ref = (glob.glob("/Users/davidchen/Documents/GitHub/SomaticSiMu/kmer_ref_count/6-mer/6-mer_chr*"))
            kmer_count = pd.read_csv(kmer_ref[0], index_col=0)['count'].fillna(0)
            for i in kmer_ref[1:-1]:
                sample = pd.read_csv(i, index_col=0)['count'].fillna(0)
                kmer_count = kmer_count + sample
                
            kmer_reference_count_dict = dict(zip(pd.read_csv(kmer_ref[0], index_col=0)["6"], kmer_count))
            
            sample_expected_count_dict = dict(zip(pd.read_csv(kmer_ref[0], index_col=0)["6"], kmer_count))
            for key in sample_expected_count_dict:
                sample_expected_count_dict[key] = sample_expected_count_dict[key] * (len(sample_seq) / sum(list(kmer_count_dict.values()))) 
            
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
            sample_sequence_file = seq_directory + "/" + cancer_type + '_' + str(i) + '_Stage_Lineage_' + str(number_of_lineages) + '.fasta'
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
            
        
#Cell 6 Loaded
print("Cell 6 (SomaticSiMu) of 8 Loaded")

# %%

from tkinter import *
from tkinter import messagebox
from tkinter import ttk 
from os.path import isfile, join
from PIL import Image , ImageTk
import threading

count_progress = 0

def multiprocessing_func( p ):
    
    arg = p[0]
    number_of_lineages = p[-1]

    try:
        
        somatic_sim(cancer_type = arg[1] , 
                reading_frame = int(arg[2]) , 
                std_outlier = int(arg[3]) , 
                number_of_lineages = number_of_lineages , 
                simulation_type = arg[4] , 
                sequence_abs_path = abs_path( arg[5].strip() , "File"), 
                slice_start = int(arg[6]) , 
                slice_end = int(arg[7]) , 
                power=int(arg[8]) , 
                syn_rate= int(arg[9]) , 
                non_syn_rate=int(arg[-1]) ,  
                ) 

    except Exception as e:
        messagebox.showinfo( 'Window title'  , 'Error :- ' + '\n ' + str(e) ,  parent= screen1 )



from itertools import product
from functools import partial

def run( *lis): 
    global arg
    arg = [i.get( ) for i in lis]
    print(arg)
    if arg.count('') > 0:
        messagebox.showinfo('Window' , 'Please enter all arguments for simulation.' , parent= screen1)
    else:          
        
        starttime= time.time()
        pool = multiprocessing.Pool( 8)
        
        my_arg  = list ( product( [arg] ,  list( range( int(arg[ 0 ]) ) )  ) )

        pool.map(  multiprocessing_func , my_arg  )

        pool.close()
        pool.join()
        print('That took {} seconds'.format(time.time() - starttime))


class New_wind_open:
    def __init__(self , root , txt ):
        self.root  , self.txt = root , txt
        self.new_wind  =Toplevel( self.root )
        self.new_wind.title('SomaticSiMu')
        self.new_wind.geometry('500x500')
        self.new_wind.resizable(0,0)

        self.frame4 = Frame( self.new_wind , height=200, width=200 )
        self.frame4.config(background="DarkBlue")
        self.frame4.pack(fill='both')

        Label( self.frame4 , text = 'SomaticSiMu' , bg = 'DarkBlue',  fg = 'white' , font = ( 'verdana' , 20 ) ).pack()
        Label(self.new_wind ).pack()
        Label( self.new_wind , text =  self.txt  , fg = 'Black' , font =('Helvitka' , 12  , 'bold') ).pack()


def make(root , txt):
    my_obj = New_wind_open( root , txt )


def open_img(root , name):
    new_img_wind = Toplevel(root)
    new_img_wind.geometry('700x550')
    new_img_wind.resizable(0,0)

    frame6 = Frame( new_img_wind , height=150, width=150 )
    frame6.config(background="DarkBlue")
    frame6.pack(fill='both')
    Label(frame6 , text = 'SomaticSiMu' ,bg = 'Darkblue' , fg = 'White' , font = ('verdana' , 14) ).pack()
    
    img = Image.open( name)
    img = img.resize((600,450), Image.ANTIALIAS) 
    photo = ImageTk.PhotoImage(img)
    Label( new_img_wind  ,  image = photo).pack()
    new_img_wind.mainloop()

def gui():
    
    global screen1 

    global Number_of_lineages , Cancer_type , Reading_frame ,  Std_outlier , Simulation_type , Sequence_abs_path , Slice_start , Slice_end , Power , Syn_rate , Non_syn_rate 

    screen1 = Toplevel() 
    screen1.title('SomaticSiMu')
    screen1.geometry('600x510')

    screen1.resizable(0,0)

    frame3 = Frame( screen1 , height=200, width=200 )
    frame3.config(background="DarkBlue")
    frame3.pack(fill='both')

    Label( frame3 , text = 'SomaticSiMu' , bg = 'DarkBlue',  fg = 'white' , font = ( 'verdana' , 30 ) ).pack()

    Button( screen1 , text = 'Sequence abs path:',bg= 'White' , fg = 'Black', font = ('Helvitka' , 10 ,'bold' ) , command = lambda : open_img( screen1 ,  abs_path('sequence_abs_path.png', "File") )  ).place(x = 80 , y = 80)
    Sequence_abs_path = Entry( screen1 , font = ( 'verdana' , 10 ) ) #Box to type in string
    Sequence_abs_path.place(x = 80 , y = 100)

    n = StringVar() 
    Button(screen1 , text = 'Cancer_type'  , bg= 'White' , fg = 'Black' , font = ('Helvitka' , 10 , 'bold' ) , command = lambda : open_img( screen1 , abs_path('cancer_type.png', "File" ) )).place(x = 80 , y = 140)  
    Cancer_type = ttk.Combobox( screen1 , height = 12 , width = 14 ,textvariable = n) 
    Cancer_type['values'] = ( cancer_type_list ) 
    Cancer_type.place(x = 80 , y= 160) 
    Cancer_type.current() 
    
    Button( screen1 , text = 'Std Outlier'  ,bg= 'White' , fg = 'Black' , font = ('Helvitka' , 10 , 'bold' ) ,  command = lambda : open_img( screen1 , abs_path('std_outlier.png', "File" ) )).place(x = 80 , y = 200)
    Std_outlier = Entry( screen1 , font = ( 'verdana' , 10 ) )  #Box to type in float number
    Std_outlier.place( x = 80 , y= 220 )

    Button( screen1 , text = 'Reading frame',bg= 'White' , fg = 'Black' , font = ('Helvitka' , 10  , 'bold' ) , command = lambda : open_img(screen1 , abs_path('reading_frame.png', "File")  )).place(x = 80 , y = 260)
    n2=StringVar()
    Reading_frame = ttk.Combobox( screen1 , height = 12 , width = 14 ,textvariable = n2)   #Drop down menu for the integers 1, 2 and 3
    Reading_frame['values'] = ( [1,2,3] ) 
    Reading_frame.place(x = 80 , y = 280)
    Reading_frame.current() 

    Button( screen1 , text = 'Number_of_lineages',bg= 'White' , fg = 'Black' , font = ('Helvitka' , 10 , 'bold' ) , command = lambda: open_img( screen1 , abs_path('number_of_lineage.png', "File" )  )).place(x = 80 , y = 320)
    Number_of_lineages =  Entry( screen1 , font = ( 'verdana' , 10 ) ) #Box to type in a whole number integer
    Number_of_lineages.place( x =  80 , y = 340 )

    Button( screen1 , text = 'Simulation_type' , bg= 'White' , fg = 'Black' , font = ('Helvitka' , 10 , 'bold' ) ,  command = lambda : open_img( screen1 , abs_path('simulation_type.png', "File" ) ) ).place(x = 80 , y = 380)
    n3=StringVar()
    Simulation_type = ttk.Combobox( screen1 , height = 12 , width = 14 ,textvariable = n3) #Drop down menu for the simulation of temporal or end stage mutation signature activities
    Simulation_type['values'] = (["temporal", "end"])
    Simulation_type.place( x=80, y=400)
    Simulation_type.current()
  
    Button( screen1 , text = 'Slice_start', bg= 'White' , fg = 'Black' , font = ('Helvitka' , 10 , 'bold' ) ,  command = lambda : open_img( screen1 , abs_path('slice_start.png', "File" ) )  ).place(x = 360  , y = 80)
    Slice_start = Entry( screen1 , font = ( 'verdana' , 10 ) )  #Box to type in integer number
    Slice_start.place(x = 360  , y = 100)

    Button( screen1 , text = 'Slice_end' ,bg= 'White' , fg = 'Black' , font = ('Helvitka' , 10 , 'bold' ),   command = lambda : open_img( screen1 , abs_path('slice_end.png', "File" ) ) ).place(x = 360  , y = 140)
    Slice_end = Entry( screen1 , font = ( 'verdana' , 10 ) )  #Box to type in integer number
    Slice_end.place(x = 360  , y = 160)

    Button( screen1 , text = 'Power' , bg= 'White' , fg = 'Black', font = ('Helvitka' , 10 , 'bold' ) ,  command = lambda : open_img( screen1 , abs_path('power.png', "File" ) )  ).place(x = 360  , y = 200)
    Power =  Entry( screen1 , font = ( 'verdana' , 10 ) )  #Box to type in float number
    Power.place(x = 360 , y = 220)

    Button( screen1 , text = 'Syn_rate',bg= 'White' , fg = 'Black'  , font = ('Helvitka' , 10 , 'bold' ) ,  command = lambda : open_img( screen1 , abs_path('syn_rate.png',"File" )  ) ).place(x = 360  , y = 260)
    Syn_rate =  Entry( screen1 , font = ( 'verdana' , 10 ) ) #Box to type in float number
    Syn_rate.place(x = 360 , y  = 280)

    Button( screen1 , text = 'Non_syn_rate:' , bg= 'White' , fg = 'Black' , font = ('Helvitka' , 10 , 'bold' ) ,  command = lambda : open_img( screen1 , abs_path('non_syn_rate.png', "File")  ) ).place(x = 360 , y = 320)
    Non_syn_rate = Entry( screen1 , font = ( 'verdana' , 10 ) )  #Box to type in float number
    Non_syn_rate.place(x = 360 , y = 340)

    Button( screen1 , text = 'Simulate' , bg='White' , fg=  'Black'  , font = ('Helvitka' , 18 , 'bold') , command = lambda : run( Number_of_lineages , Cancer_type , Reading_frame ,  Std_outlier , Simulation_type , Sequence_abs_path , Slice_start , Slice_end , Power , Syn_rate , Non_syn_rate  ) ).place(x = 360 , y = 400 )
    
#Cell 7 Loaded
print("Cell 7 (Simulation) of 8 Loaded")
#%%


import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.patches as patches
from matplotlib.figure import Figure
matplotlib.use("TkAgg")

def gui2():
    global screen2
    screen2 = Toplevel() 
    screen2.title('SomaticSiMu')
    screen2.resizable(0,0)
    screen2.geometry('500x400')
    frame2 = Frame( screen2 , height=200, width=200 )
    frame2.config(background="DarkBlue")
    frame2.pack(fill='both')
    
    Label( frame2 , text = 'SomaticSiMu' , bg = 'DarkBlue' , fg = 'white' , font = ( 'verdana' , 30 ) ).pack()

    n3 = StringVar() 
    Button(screen2 , text = 'Cancer_type' , bg= 'White' , fg = 'Black'  , font = ('Helvitka' , 12 , 'bold' ) ,  command = lambda : open_img( screen2 , abs_path('viz_cancer_type.png', "File" ) )).place(x = 60 , y = 80)  
    Cancer_type2 = ttk.Combobox( screen2 , height = 12 , width = 14 ,textvariable = n3) 
    Cancer_type2['values'] = ( cancer_type_list ) 
    Cancer_type2.place(x =60 , y= 100) 
    Cancer_type2.current() 

    n5 = StringVar() 
    Button(screen2 , text = 'Visualization_type',bg= 'White' , fg = 'Black'  , font = ('Helvitka' , 12 , 'bold' ) , command = lambda : open_img( screen2 , abs_path('visualization_type.png', "File" ))).place(x = 290 , y = 80)  
    Sim_type2 = ttk.Combobox( screen2 , height = 12 , width = 14 ,textvariable = n5)
    Sim_type2['values'] = ( ['Early', 'Late' , 'End'] )  
    Sim_type2.place(x = 290 , y= 100) 
    Sim_type2.current() 

    Button( screen2 , text = 'Gen_start', bg= 'White' , fg = 'Black' , font = ('Helvitka' , 12 , 'bold' ),command = lambda : open_img( screen2 , abs_path('viz_gen_start.png', "File" ) ) ).place(x = 60 , y = 150)
    Gen_start = Entry( screen2 , font = ( 'verdana' , 10 ) )  #Box to type in integer number
    Gen_start.place(x = 60 , y = 170)


    Button( screen2 , text = 'Gen_end', bg= 'White' , fg = 'Black'  , font = ('Helvitka' , 12 , 'bold' ), command = lambda : open_img( screen2 , abs_path('viz_gen_end.png', "File" ) ) ).place(x = 290 , y = 150)
    Gen_end = Entry( screen2 , font = ( 'verdana' , 10 ) )  #Box to type in integer number
    Gen_end.place(x = 290 , y = 170)

    n4 = StringVar() 
    Button(screen2 , text = 'Mut_type',bg= 'White' , fg = 'Black'  , font = ('Helvitka' , 12 , 'bold' ) , command = lambda : open_img( screen2 , abs_path('viz_mut_type.png', "File" ))).place(x = 60 , y = 220)  
    Mut_type = ttk.Combobox( screen2 , height = 12 , width = 14 ,textvariable = n4)
    Mut_type['values'] = ( ['SBS', 'DBS' , 'Insertion' , 'Deletion', 'Mutation Burden'] )  
    Mut_type.place(x = 60 , y= 240) 
    Mut_type.current() 

    Button( screen2 , text = 'Visualize Graph' , bg='White' , fg=  'Black'  , font = ('Helvitka' , 18 ,'bold') , command = lambda : show( screen2 ,  Cancer_type2 , Sim_type2 , Gen_start , Gen_end , Mut_type )  ).place(x = 290 , y = 240 )

    screen2.mainloop()



def mut_catalog(cancer_type, simulation_type, gen_start, gen_end, mut_type):
    try:
        myapp = Tk()
        myapp.geometry('1600x900')
        
        myapp.title('SomaticSiMu')

        frame1 = Frame( myapp  )
      
        frame1.pack( fill='both')

        if mut_type.lower() == "sbs":
            
            graph_data = pd.DataFrame(index=(range(96)), columns=range(gen_start, gen_end+1))
            
            for gen in range(gen_start, gen_end+1):
                try:
                    df_read_in = pd.read_csv(abs_path(cancer_type + '_' + str(simulation_type) + '_Stage_Lineage_' + str(gen) + '_' + mut_type.lower() + '_freq_table.csv', 'File'))
                    graph_data[gen] = df_read_in['Frequency']
                except:
                    print(cancer_type + '_' + str(simulation_type) + '_Stage_Lineage_' + str(gen) + '_' + mut_type.lower() + '_freq_table.csv does not exist!')
                    
            graph_data['SubType'] = df_read_in['SubType']
            graph_data['Type'] = df_read_in['Type']
            
            for gen in range(gen_start, gen_end+1):
                 graph_data[gen] = graph_data[gen].div(graph_data[gen].sum())
            
            graph_data.sort_values(by=["Type","SubType"],inplace=True )
            
            if len(range(gen_start, gen_end+1)) > 1:
                
                upper_error = [graph_data.iloc[mut,:-2].quantile(0.95) for mut in range(96)]
                per_95 = []
                for i in upper_error:
                    per_95.append(i - graph_data.iloc[upper_error.index(i), :-2].mean(axis=0))
                
                
                lower_error = [graph_data.iloc[mut,:-2].quantile(0.05) for mut in range(96)]
                per_5 = []
                counter = 0
                for i in lower_error:
                    if i == 0:
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0))
                    else:  
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0) - i)
                    counter += 1
                
            
                fig = plt.figure(figsize=(50, 10))
                ax = fig.add_subplot(111)
                
                N = 96
                ind = np.arange(96)   
                
                color_list = [(0.416, 0.733, 0.918), (0,0,0), (0.765, 0.172, 0.157), (0.785, 0.785, 0.785), (0.678, 0.808, 0.412), (0.878, 0.773, 0.769)]
                
                p1 = ax.bar(range(0,16), graph_data.iloc[0:16, :-2].mean(axis=1), bottom=0, color = color_list[0], yerr = (per_5[0:16], per_95[0:16]))
                p2 = ax.bar(range(16,32), graph_data.iloc[16:32, :-2].mean(axis=1), bottom=0, color = color_list[1], yerr = (per_5[16:32], per_95[16:32]))
                p3 = ax.bar(range(32,48), graph_data.iloc[32:48, :-2].mean(axis=1), bottom=0, color = color_list[2], yerr = (per_5[32:48], per_95[32:48]))
                p4 = ax.bar(range(48,64), graph_data.iloc[48:64, :-2].mean(axis=1), bottom=0, color = color_list[3], yerr = (per_5[48:64], per_95[48:64]))
                p5 = ax.bar(range(64,80), graph_data.iloc[64:80, :-2].mean(axis=1), bottom=0, color = color_list[4], yerr = (per_5[64:80], per_95[64:80]))
                p6 = ax.bar(range(80,96), graph_data.iloc[80:96, :-2].mean(axis=1), bottom=0, color = color_list[5], yerr = (per_5[80:96], per_95[80:96]))
                  
            else:
                upper_error = [graph_data.iloc[mut,:-2].quantile(0.95) for mut in range(96)]
                lower_error = [graph_data.iloc[mut,:-2].quantile(0.05) for mut in range(96)]
            
                fig = plt.figure(figsize=(50, 10))
                ax = fig.add_subplot(111)
                
                N = 96
                ind = np.arange(96)   
                
                color_list = [(0.416, 0.733, 0.918), (0,0,0), (0.765, 0.172, 0.157), (0.785, 0.785, 0.785), (0.678, 0.808, 0.412), (0.878, 0.773, 0.769)]
                
                p1 = ax.bar(range(0,16), graph_data.iloc[0:16, :-2].mean(axis=1), bottom=0, color = color_list[0])
                p2 = ax.bar(range(16,32), graph_data.iloc[16:32, :-2].mean(axis=1), bottom=0, color = color_list[1])
                p3 = ax.bar(range(32,48), graph_data.iloc[32:48, :-2].mean(axis=1), bottom=0, color = color_list[2])
                p4 = ax.bar(range(48,64), graph_data.iloc[48:64, :-2].mean(axis=1), bottom=0, color = color_list[3])
                p5 = ax.bar(range(64,80), graph_data.iloc[64:80, :-2].mean(axis=1), bottom=0, color = color_list[4])
                p6 = ax.bar(range(80,96), graph_data.iloc[80:96, :-2].mean(axis=1), bottom=0, color = color_list[5])
                
                
                
            y_limit = ax.get_ylim()[1]
            
            rect1 = patches.Rectangle((1, y_limit + y_limit/50), 15, y_limit/20, color = color_list[0], clip_on=False) 
            rect2 = patches.Rectangle((17, y_limit + y_limit/50), 15, y_limit/20, color = color_list[1], clip_on=False) 
            rect3 = patches.Rectangle((33, y_limit + y_limit/50), 15, y_limit/20, color = color_list[2], clip_on=False) 
            rect4 = patches.Rectangle((49, y_limit + y_limit/50), 15, y_limit/20, color = color_list[3], clip_on=False) 
            rect5 = patches.Rectangle((65, y_limit + y_limit/50), 15, y_limit/20, color = color_list[4], clip_on=False) 
            rect6 = patches.Rectangle((81, y_limit + y_limit/50), 15, y_limit/20, color = color_list[5], clip_on=False) 
            
            plt.text(7, y_limit + y_limit/9, "C>A", fontsize=12, weight="bold")
            plt.text(23, y_limit + y_limit/9, "C>G", fontsize=12, weight="bold")
            plt.text(39, y_limit + y_limit/9, "C>T", fontsize=12, weight="bold")
            plt.text(55, y_limit + y_limit/9, "T>A", fontsize=12, weight="bold")
            plt.text(71, y_limit + y_limit/9, "T>C", fontsize=12, weight="bold")
            plt.text(87, y_limit + y_limit/9, "T>G", fontsize=12, weight="bold")
            
            ax.add_patch(rect1)
            ax.add_patch(rect2)
            ax.add_patch(rect3)
            ax.add_patch(rect4)
            ax.add_patch(rect5)
            ax.add_patch(rect6)
            
            ax.set_xlabel('SubType', fontsize=15, weight="bold")
            ax.yaxis.set_label_text('Percentage of Single Base Substitutions', fontsize=15, weight="bold")
            
            ax.set_xticks(ind)
            ax.set_xticklabels(graph_data['SubType'], fontsize=6)
            
            plt.xticks(rotation=90)
            ax.margins(x=0)
            ax.set_title(" ", pad=140)
            
            ax.grid(axis = 'y', color=color_list[3], linestyle='-', linewidth = 1)
            ax.set_axisbelow(True)
            ax.xaxis.labelpad= 10
            ax.yaxis.labelpad= 10
            
            canvas = FigureCanvasTkAgg(fig , frame1)
            
            canvas.get_tk_widget().pack( fill=BOTH, expand=True)
  
        if mut_type.lower() == "dbs":
            
            graph_data = pd.DataFrame(index=(range(78)), columns=range(gen_start, gen_end+1))
            
            for gen in range(gen_start, gen_end+1):
                try:
                    df_read_in = pd.read_csv(abs_path(cancer_type + '_' + str(simulation_type) + '_Stage_Lineage_' + str(gen) + '_' + mut_type.lower() + '_freq_table.csv', 'File'))
                    graph_data[gen] = df_read_in['Frequency']
                except:
                    print(cancer_type + '_' + str(simulation_type) + '_Stage_Lineage_' + str(gen) + '_' + mut_type.lower() + '_freq_table.csv does not exist!')
                    
            graph_data['2mer_index'] = df_read_in['2mer_index']
            graph_data['Mutation Type'] = df_read_in['Mutation Type']
            
            for gen in range(gen_start, gen_end+1):
                 graph_data[gen] = graph_data[gen].div(graph_data[gen].sum())
            
            graph_data.sort_values(by=["Mutation Type", "2mer_index"],inplace=True)
            

            if len(range(gen_start, gen_end+1)) > 1:
                upper_error = [graph_data.iloc[mut,:-2].quantile(0.95) for mut in range(78)]
                per_95 = []
                counter = 0
                for i in upper_error:
                    if i == 1:
                        per_95.append(0)
                    else:
                        per_95.append(i - graph_data.iloc[counter, :-2].mean(axis=0))
                        counter += 1
                per_95 =  [0 if i < 0 else i for i in per_95]
                       
                lower_error = [graph_data.iloc[mut,:-2].quantile(0.05) for mut in range(78)]
                per_5 = []
                counter = 0
                for i in lower_error:
                    if i == 0:
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0))
                    else:  
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0) - i)
                    counter += 1

                fig = plt.figure(figsize=(50, 10))
                ax = fig.add_subplot(111)
                
                N = 78
                ind = np.arange(78)   
                
                color_list = [(0.416, 0.733, 0.918), 
                              (0.227, 0.404, 0.773), 
                              (0.682, 0.796, 0.420), 
                              (0.224, 0.396, 0.0941), 
                              (0.906, 0.592, 0.592), 
                              (0.761, 0.173, 0.170),
                              (0.918, 0.686, 0.424),
                              (0.89, 0.502, 0.118),
                              (0.737, 0.588, 0.984),
                              (0.263, 0.0118, 0.584)]
                
                p1 = ax.bar(range(0,9), graph_data.iloc[0:9, :-2].mean(axis=1), bottom=0, color = color_list[0], yerr = (per_5[0:9], per_95[0:9]))
                p2 = ax.bar(range(9,15), graph_data.iloc[9:15, :-2].mean(axis=1), bottom=0, color = color_list[1], yerr = (per_5[9:15], per_95[9:15]))
                p3 = ax.bar(range(15,24), graph_data.iloc[15:24, :-2].mean(axis=1), bottom=0, color = color_list[2], yerr = (per_5[15:24], per_95[15:24]))
                p4 = ax.bar(range(24,30), graph_data.iloc[24:30, :-2].mean(axis=1), bottom=0, color = color_list[3], yerr = (per_5[24:30], per_95[24:30]))
                p5 = ax.bar(range(30,39), graph_data.iloc[30:39, :-2].mean(axis=1), bottom=0, color = color_list[4], yerr = (per_5[30:39], per_95[30:39]))
                p6 = ax.bar(range(39,45), graph_data.iloc[39:45, :-2].mean(axis=1), bottom=0, color = color_list[5], yerr = (per_5[39:45], per_95[39:45]))
                p7 = ax.bar(range(45,51), graph_data.iloc[45:51, :-2].mean(axis=1), bottom=0, color = color_list[6], yerr = (per_5[45:51], per_95[45:51]))
                p8 = ax.bar(range(51,60), graph_data.iloc[51:60, :-2].mean(axis=1), bottom=0, color = color_list[7], yerr = (per_5[51:60], per_95[51:60]))
                p9 = ax.bar(range(60,69), graph_data.iloc[60:69, :-2].mean(axis=1), bottom=0, color = color_list[8], yerr = (per_5[60:69], per_95[60:69]))
                p10 = ax.bar(range(69,78), graph_data.iloc[69:78, :-2].mean(axis=1), bottom=0, color = color_list[9], yerr = (per_5[69:78], per_95[69:78]))
                
            else:
                fig = plt.figure(figsize=(50, 10))
                ax = fig.add_subplot(111)
                
                N = 78
                ind = np.arange(78)   
                
                color_list = [(0.416, 0.733, 0.918), 
                              (0.227, 0.404, 0.773), 
                              (0.682, 0.796, 0.420), 
                              (0.224, 0.396, 0.0941), 
                              (0.906, 0.592, 0.592), 
                              (0.761, 0.173, 0.170),
                              (0.918, 0.686, 0.424),
                              (0.89, 0.502, 0.118),
                              (0.737, 0.588, 0.984),
                              (0.263, 0.0118, 0.584)]
                
                p1 = ax.bar(range(0,9), graph_data.iloc[0:9, 2], bottom=0, color = color_list[0])
                p2 = ax.bar(range(9,15), graph_data.iloc[9:15, 2], bottom=0, color = color_list[1])
                p3 = ax.bar(range(15,24), graph_data.iloc[15:24, 2], bottom=0, color = color_list[2])
                p4 = ax.bar(range(24,30), graph_data.iloc[24:30, 2], bottom=0, color = color_list[3])
                p5 = ax.bar(range(30,39), graph_data.iloc[30:39, 2], bottom=0, color = color_list[4])
                p6 = ax.bar(range(39,45), graph_data.iloc[39:45, 2], bottom=0, color = color_list[5])
                p7 = ax.bar(range(45,51), graph_data.iloc[45:51, 2], bottom=0, color = color_list[6])
                p8 = ax.bar(range(51,60), graph_data.iloc[51:60, 2], bottom=0, color = color_list[7])
                p9 = ax.bar(range(60,69), graph_data.iloc[60:69, 2], bottom=0, color = color_list[8])
                p10 = ax.bar(range(69,78), graph_data.iloc[69:78, 2], bottom=0, color = color_list[9])

            y_limit = ax.get_ylim()[1]
            
            rect1 = patches.Rectangle((0, y_limit + y_limit/50), 8.5, y_limit/12, color = color_list[0], clip_on=False) 
            rect2 = patches.Rectangle((9, y_limit + y_limit/50), 5.5, y_limit/12, color = color_list[1], clip_on=False) 
            rect3 = patches.Rectangle((15, y_limit + y_limit/50), 8.5, y_limit/12, color = color_list[2], clip_on=False) 
            rect4 = patches.Rectangle((24, y_limit + y_limit/50), 5.5, y_limit/12, color = color_list[3], clip_on=False) 
            rect5 = patches.Rectangle((30, y_limit + y_limit/50), 8.5, y_limit/12, color = color_list[4], clip_on=False) 
            rect6 = patches.Rectangle((39, y_limit + y_limit/50), 5.5, y_limit/12, color = color_list[5], clip_on=False) 
            rect7 = patches.Rectangle((45, y_limit + y_limit/50), 5.5, y_limit/12, color = color_list[6], clip_on=False) 
            rect8 = patches.Rectangle((51, y_limit + y_limit/50), 8.5, y_limit/12, color = color_list[7], clip_on=False) 
            rect9 = patches.Rectangle((60, y_limit + y_limit/50), 8.5, y_limit/12, color = color_list[8], clip_on=False) 
            rect10 = patches.Rectangle((69, y_limit + y_limit/50), 8.5, y_limit/12, color = color_list[9], clip_on=False) 


            plt.text(2, y_limit + y_limit/9, "AC>NN", fontsize=10, weight="bold")
            plt.text(9.75, y_limit + y_limit/9, "AT>NN", fontsize=10, weight="bold")
            plt.text(17, y_limit + y_limit/9, "CC>NN", fontsize=10, weight="bold")
            plt.text(25, y_limit + y_limit/9, "CG>NN", fontsize=10, weight="bold")
            plt.text(32, y_limit + y_limit/9, "CT>NN", fontsize=10, weight="bold")
            plt.text(39.625, y_limit + y_limit/9, "GC>NN", fontsize=10, weight="bold")
            plt.text(46, y_limit + y_limit/9, "TA>NN", fontsize=10, weight="bold")
            plt.text(53, y_limit + y_limit/9, "TC>NN", fontsize=10, weight="bold")
            plt.text(61.75, y_limit + y_limit/9, "TG>NN", fontsize=10, weight="bold")
            plt.text(71.25, y_limit + y_limit/9, "TT>NN", fontsize=10, weight="bold")        
        

            ax.add_patch(rect1)
            ax.add_patch(rect2)
            ax.add_patch(rect3)
            ax.add_patch(rect4)
            ax.add_patch(rect5)
            ax.add_patch(rect6)
            ax.add_patch(rect7)
            ax.add_patch(rect8)
            ax.add_patch(rect9)
            ax.add_patch(rect10)
            
            ax.set_xlabel('SubType', fontsize=15, weight="bold")
            ax.yaxis.set_label_text('Percentage of Double Base Substitutions', fontsize=15, weight="bold")
            
            ax.set_xticks(ind)
            ax.set_xticklabels(graph_data['Mutation Type'].str[3:], fontsize=8)
            
            plt.xticks(rotation=90)
            ax.margins(x=0)
            ax.set_title(" ", pad=140)
            
            ax.grid(axis = 'y', color=(0.785, 0.785, 0.785), linestyle='-', linewidth = 1)
            ax.set_axisbelow(True)
            ax.xaxis.labelpad= 10
            ax.yaxis.labelpad= 10
            
            canvas = FigureCanvasTkAgg(fig , master=frame1)
            
            canvas.get_tk_widget().pack( fill=BOTH, expand=True)
        
        if mut_type.lower() == "insertion":
            
            graph_data = pd.DataFrame(index=(range(12)), columns=range(gen_start, gen_end+1))
            
            for gen in range(gen_start, gen_end+1):
                try:
                    df_read_in = pd.read_csv(abs_path(cancer_type + '_' + str(simulation_type) + '_Stage_Lineage_' + str(gen) + '_ins_freq_table.csv', 'File'))
                    graph_data[gen] = df_read_in['Frequency']
                except:
                    print(cancer_type + '_' + str(simulation_type) + '_Stage_Lineage_' + str(gen) + '_ins_freq_table.csv does not exist!')
                    
            graph_data['Index'] = df_read_in['Index']
            graph_data['Mutation Type'] = df_read_in['Mutation Type']
            
            for gen in range(gen_start, gen_end+1):
                 graph_data[gen] = graph_data[gen].div(graph_data[gen].sum())
            
            graph_data.sort_values(by=["Index","Mutation Type"],inplace=True )
            

            if len(range(gen_start, gen_end+1)) > 1:

                upper_error = [graph_data.iloc[mut,:-2].quantile(0.95) for mut in range(12)]
                per_95 = []
                counter = 0
                for i in upper_error:
                    if i == 1:
                        per_95.append(0)
                    else:
                        per_95.append(i - graph_data.iloc[counter, :-2].mean(axis=0))
                        counter += 1
                per_95 =  [0 if i < 0 else i for i in per_95]
                       
                lower_error = [graph_data.iloc[mut,:-2].quantile(0.05) for mut in range(12)]
                per_5 = []
                counter = 0
                for i in lower_error:
                    if i == 0:
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0))
                    else:  
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0) - i)
                    counter += 1

                fig = plt.figure(figsize=(50, 10))
                ax = fig.add_subplot(111)
                
                ind = np.arange(12)   
                
                color_list = [(0.949, 0.753, 0.478), 
                              (0.937, 0.512, 0.2)]
                
                p1 = ax.bar(range(0,6), graph_data.iloc[0:6, :-2].mean(axis=1), bottom=0, color = color_list[0],  yerr = (per_5[0:6], per_95[0:6]))
                p2 = ax.bar(range(6,12), graph_data.iloc[6:12, :-2].mean(axis=1), bottom=0, color = color_list[1],  yerr = (per_5[6:12], per_95[6:12]))
                
            else:
                fig = plt.figure(figsize=(50, 10))
                ax = fig.add_subplot(111)
                
                ind = np.arange(12)   
                
                color_list = [(0.949, 0.753, 0.478), 
                              (0.937, 0.512, 0.2)]
                
                p1 = ax.bar(range(0,6), graph_data.iloc[0:6, :-2].mean(axis=1), bottom=0, color = color_list[0])
                p2 = ax.bar(range(6,12), graph_data.iloc[6:12, :-2].mean(axis=1), bottom=0, color = color_list[1])
                
               
            y_limit = ax.get_ylim()[1]
            
            rect1 = patches.Rectangle((-0.25, y_limit + y_limit/50), 5.5, y_limit/10, color = color_list[0], clip_on=False) 
            rect2 = patches.Rectangle((5.75, y_limit + y_limit/50), 5.5, y_limit/10, color = color_list[1], clip_on=False) 
            
            
            plt.text(2.5, y_limit + y_limit/20, "C", fontsize=30, weight="bold", color="black")
            plt.text(8.5, y_limit + y_limit/20, "T", fontsize=30, weight="bold", color="black")
            
            ax.add_patch(rect1)
            ax.add_patch(rect2)
            
            ax.set_xlabel('SubType', fontsize=15, weight="bold")
            ax.yaxis.set_label_text('Percentage of Insertions', fontsize=15, weight="bold")
            
            ax.set_xticks(ind)
            ax.set_xticklabels(graph_data['Index'], fontsize=8)
            
            plt.xticks(rotation=90)
            ax.margins(x=0)
            ax.set_title("1 Base Pair Deletion at Repeats", fontsize=15, weight="bold",  pad=100)
            
            ax.grid(axis = 'y', color=(0.785, 0.785, 0.785), linestyle='-', linewidth = 1)
            ax.set_axisbelow(True)
            ax.xaxis.labelpad= 10
            ax.yaxis.labelpad= 10

            canvas = FigureCanvasTkAgg(fig , master=frame1)
            canvas.get_tk_widget().pack( fill=BOTH, expand=True)
            
        if mut_type.lower() == "deletion":
            
            graph_data = pd.DataFrame(index=(range(12)), columns=range(gen_start, gen_end+1))
            
            for gen in range(gen_start, gen_end+1):
                try:
                    df_read_in = pd.read_csv(abs_path(cancer_type + '_' + str(simulation_type) + '_Stage_Lineage_' + str(gen) + '_del_freq_table.csv', 'File'))
                    graph_data[gen] = df_read_in['Frequency']
                except:
                    print(cancer_type + '_' + str(simulation_type) + '_Stage_Lineage_' + str(gen) + '_del_freq_table.csv does not exist!')
                    
            graph_data['Index'] = df_read_in['Index']
            graph_data['Mutation Type'] = df_read_in['Mutation Type']
            
            for gen in range(gen_start, gen_end+1):
                 graph_data[gen] = graph_data[gen].div(graph_data[gen].sum())
            
            graph_data.sort_values(by=["Index","Mutation Type"],inplace=True )
            

            if len(range(gen_start, gen_end+1)) > 1:

                upper_error = [graph_data.iloc[mut,:-2].quantile(0.95) for mut in range(12)]
                per_95 = []
                counter = 0
                for i in upper_error:
                    if i == 1:
                        per_95.append(0)
                    else:
                        per_95.append(i - graph_data.iloc[counter, :-2].mean(axis=0))
                        counter += 1
                per_95 =  [0 if i < 0 else i for i in per_95]
                       
                lower_error = [graph_data.iloc[mut,:-2].quantile(0.05) for mut in range(12)]
                per_5 = []
                counter = 0
                for i in lower_error:
                    if i == 0:
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0))
                    else:  
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0) - i)
                    counter += 1
            
                fig = plt.figure(figsize=(50, 10))
                ax = fig.add_subplot(111)
            
                ind = np.arange(12)   
            
                color_list = [(0.727, 0.855, 0.569), 
                          (0.345, 0.612, 0.247)]
            
                p1 = ax.bar(range(0,6), graph_data.iloc[0:6, :-2].mean(axis=1), bottom=0, color = color_list[0],  yerr = (per_5[0:6], per_95[0:6]))
                p2 = ax.bar(range(6,12), graph_data.iloc[6:12, :-2].mean(axis=1), bottom=0, color = color_list[1],  yerr = (per_5[6:12], per_95[6:12]))
                
            else:
                fig = plt.figure(figsize=(50, 10))
                ax = fig.add_subplot(111)
            
                ind = np.arange(12)   
            
                color_list = [(0.727, 0.855, 0.569), 
                          (0.345, 0.612, 0.247)]
            
                p1 = ax.bar(range(0,6), graph_data.iloc[0:6, :-2].mean(axis=1), bottom=0, color = color_list[0])
                p2 = ax.bar(range(6,12), graph_data.iloc[6:12, :-2].mean(axis=1), bottom=0, color = color_list[1])
                
               
            y_limit = ax.get_ylim()[1]
            
            rect1 = patches.Rectangle((-0.25, y_limit + y_limit/50), 5.5, y_limit/10, color = color_list[0], clip_on=False) 
            rect2 = patches.Rectangle((5.75, y_limit + y_limit/50), 5.5, y_limit/10, color = color_list[1], clip_on=False) 
            
            
            plt.text(2.5, y_limit + y_limit/20, "C", fontsize=30, weight="bold", color="black")
            plt.text(8.5, y_limit + y_limit/20, "T", fontsize=30, weight="bold", color="black")
            
            ax.add_patch(rect1)
            ax.add_patch(rect2)
            
            ax.set_xlabel('SubType', fontsize=15, weight="bold")
            ax.yaxis.set_label_text('Percentage of Deletions', fontsize=15, weight="bold")
            
            ax.set_xticks(ind)
            ax.set_xticklabels(graph_data['Index'], fontsize=8)
            
            plt.xticks(rotation=90)
            ax.margins(x=0)
            ax.set_title("1 Base Pair Deletion at Repeats", fontsize=15, weight="bold",  pad=100)
            
            ax.grid(axis = 'y', color=(0.785, 0.785, 0.785), linestyle='-', linewidth = 1)
            ax.set_axisbelow(True)
            ax.xaxis.labelpad= 10
            ax.yaxis.labelpad= 10

            canvas = FigureCanvasTkAgg(fig , master=frame1)
            canvas.get_tk_widget().pack(fill=BOTH, expand=True)


        if mut_type.lower() == "mutation burden":
           
            mutation_burden = pd.DataFrame(index = range(gen_start, gen_end+1), columns = ['sbs', 'dbs', 'ins', 'del'])
            
            for gen in range(gen_start, gen_end+1):
                try:
                    sbs_burden = pd.read_csv(abs_path(cancer_type + '_' + str(simulation_type) + '_Stage_Lineage_' + str(gen) + '_' + "sbs" + '_freq_table.csv', 'File'))['Frequency'].sum()
                    dbs_burden = pd.read_csv(abs_path(cancer_type + '_' + str(simulation_type) + '_Stage_Lineage_' + str(gen) + '_' + "dbs" + '_freq_table.csv', 'File'))['Frequency'].sum()
                    ins_burden = pd.read_csv(abs_path(cancer_type + '_' + str(simulation_type) + '_Stage_Lineage_' + str(gen) + '_' + "ins" + '_freq_table.csv', 'File'))['Frequency'].sum()
                    del_burden = pd.read_csv(abs_path(cancer_type + '_' + str(simulation_type) + '_Stage_Lineage_' + str(gen) + '_' + "del" + '_freq_table.csv', 'File'))['Frequency'].sum()
                
                    mutation_burden.iloc[list(range(gen_start, gen_end+1)).index(gen), 0] = sbs_burden
                    mutation_burden.iloc[list(range(gen_start, gen_end+1)).index(gen), 1] = dbs_burden
                    mutation_burden.iloc[list(range(gen_start, gen_end+1)).index(gen), 2] = ins_burden
                    mutation_burden.iloc[list(range(gen_start, gen_end+1)).index(gen), 3] = del_burden
                    
                except:
                    pass
                
            plt.style.use('seaborn-whitegrid')
            #fig, ax = plt.subplots(figsize=(10, 10))
            ax = mutation_burden.plot.bar(stacked=True, color= ["#E74C3C", "#3498DB","#27AE60", "#9B59B6"], figsize=(10,7))
            ax.set_ylim(mutation_burden.min().sum() - 0.2*mutation_burden.max().sum() - (mutation_burden.min().sum() - 0.1*mutation_burden.max().sum())%10, mutation_burden.max().sum() + 0.15*mutation_burden.max().sum()+10 - (mutation_burden.max().sum() + 0.2*mutation_burden.max().sum())%10)
            ax.tick_params(axis='x', rotation=0)
            ax.set_xlabel("Sequence ID")
            ax.set_ylabel("Number of simulated mutations")
            ax.xaxis.labelpad= 10
            ax.yaxis.labelpad= 10

            #canvas = FigureCanvasTkAgg(fig , master=frame1)
            #canvas.get_tk_widget().pack(fill=BOTH, expand=True)

        if mut_type.lower() not in ["sbs", "dbs", "insertion", "deletion"]:
            print("Choose a mutation type from sbs, dbs, insertion or deletion.")
        print(  'Visualization is loading...' )    
        
        myapp.mainloop()
    except Exception as e:
        myapp.destroy()
        messagebox.showinfo('Window' , str(e) , parent= screen2 )

def show(root , *lis):
    ag = [i.get()  for i in lis]
    if ag.count( '' ) > 0:
        messagebox.showinfo('Error' , 'Please enter all arguments for visualization. ' , parent= screen2)
    else:    
        mut_catalog( ag[0] , ag[1] , int(ag[2]) , int(ag[3]) , ag[-1] )



class Main_window:
    def __init__(self):
        self.main = Tk()
        self.main.geometry('500x400')
        self.main.resizable(0,0)
        self.main.title( 'SomaticSiMu' )
        self.frame = Frame( self.main , height=200, width=200 )
        self.frame.config(background="Dark Blue")
        self.frame.pack(fill='both')
        Label( self.frame , text = 'SomaticSiMu' , bg = 'Dark Blue' , fg = 'White' , font = ( 'verdana' , 30 ) ).pack()
       
        Button(self.main , text = '  Simulation ' , bg= 'White' , fg = 'Black' , font = ('Helvitka' , 30) , command = self.func_1).place(x = 160 , y = 120 )
        Button(self.main , text = '  Visualization  ', bg = 'White' , fg = 'Black' , font = ('Helvitka' , 30) , command = self.func_2).place( x= 145, y = 220 )

        Label(self.main , text  = 'Version 2.0 by David Chen' , fg = 'Black' , font =('verdana' , 10 )).place( x = 320 , y = 370 )
                
        self.main.mainloop()

    def func_1(self):
        gui()
    def func_2(self):
        gui2()        

if __name__ ==  '__main__':
    #Cell 8 Loaded
    print("Cell 8 (Visualization) of 8 Loaded")
    print("SomaticSiMu is set up.")
    multiprocessing.set_start_method('spawn')
    obj = Main_window()
    

for i in [300, 345]:
    somatic_sim(cancer_type="Biliary-AdenoCA", reading_frame=1, std_outlier=3, number_of_lineages=i, simulation_type="end", sequence_abs_path=input_file_path, slice_start=1, slice_end=50818467,power=1, syn_rate=1, non_syn_rate=1)
    
