#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 23:26:12 2020

@author: DavidChen
"""
#%%
"""
Imports
"""
import pandas as pd
from Bio import SeqIO
import glob
import collections
import random
import numpy as np
from multiprocessing import Pool    
import functools
from collections import defaultdict
import os

#%%
"""
Helper Functions
"""



"""
Finds absolute file path from file name in working directory.
"""
def abs_path(metrics_file_name):
    #Find the relative working directory of the script
    wk_dir = os.path.dirname(os.path.realpath('__file__'))
    
    #Absolute input file path
    for root, dirs, files in os.walk(wk_dir):
        for name in files:
            if metrics_file_name == name:
                metrics_file_path = (os.path.abspath(os.path.join(root, name))) 
    
    return metrics_file_path

"""
Count kmers. Outputs dictionary.
"""
def count_kmers(data, k):
    global d
    d = collections.defaultdict(int)
    for i in range(len(data)-(k-1)):
        d[data[i:i+k]] +=1
    for key in list(d.keys()):
        if "N" in key:
            del d[key]
    return d

"""
Read fasta. Outputs sequence.
"""   
def read_fasta(fp):
        name, seq = None, []
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if name: yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name: yield (name, ''.join(seq))

"""
Set slice of genome sequence (inclusive of start and end of slice) as a string (high memory cost...)
"""
def seq_slice(sequence_abs_path, slice_start, slice_end):
    with open(input_file_path) as fp:
        for name, seq in read_fasta(fp):
            sequence = str(seq)
            sample = sequence[slice_start:slice_end+1]
    
    

#%%
"""
Input sequence file path
"""
input_file_path = abs_path("Homo_sapiens.GRCh38.dna.chromosome.22.fasta")


seq_slice(input_file_path, 26100000, 27100000)

[26100000:27100000]


#%%





#%%
"""
Read in input reference files:
data = Pandas Dataframe of number of mutations attributed to each signature for each sample.
cancer_type = List of unique cancer types.
prop_data = Pandas Dataframe of proportion of 96 type mutation classification for each signature.
"""



"""
SBS input reference files
"""
sbs_num_file_path = input('Input file path of PCAWG_sigProfiler_SBS_signatures_in_samples.csv: ')
#/Users/DavidChen/Desktop/Project/Reference/PCAWG_sigProfiler_SBS_signatures_in_samples.csv
sbs_prop_file_path = input('Input file path of sigProfiler_SBS_signatures.csv: ')
#/Users/DavidChen/Desktop/Project/Reference/sigProfiler_SBS_signatures.csv

sbs_num_data = pd.read_csv(sbs_num_file_path)
sbs_num_data.sort_values(by='Cancer Types', axis=0, inplace=True)
sbs_num_data.set_index(keys=['Cancer Types'], drop=False,inplace=True)

sbs_cancer_type = sbs_num_data['Cancer Types'].unique().tolist()
sbs_prop_data = pd.read_csv(sbs_prop_file_path)
  


"""
DBS input reference files
"""
dbs_num_file_path = input('Input file path of PCAWG_sigProfiler_DBS_signatures_in_samples.csv: ')
#/Users/DavidChen/Desktop/Project/Reference/PCAWG_sigProfiler_DBS_signatures_in_samples.csv
dbs_prop_file_path = input('Input file path of sigProfiler_DBS_signatures.csv: ')
#/Users/DavidChen/Desktop/Project/Reference/sigProfiler_DBS_signatures.csv

dbs_num_data = pd.read_csv(dbs_num_file_path)
dbs_num_data.sort_values(by='Cancer Types', axis=0, inplace=True)
dbs_num_data.set_index(keys=['Cancer Types'], drop=False,inplace=True)

dbs_cancer_type = dbs_num_data['Cancer Types'].unique().tolist()
dbs_prop_data = pd.read_csv(dbs_prop_file_path)



"""
Insertion/Deletion input reference files
"""
id_num_file_path = input('Input file path of PCAWG_sigProfiler_ID_signatures_in_samples.csv: ')
#/Users/DavidChen/Desktop/Project/Reference/PCAWG_SigProfiler_ID_signatures_in_samples.csv
id_prop_file_path = input('Input file path of sigProfiler_ID_signatures.csv: ')
#/Users/DavidChen/Desktop/Project/Reference/sigProfiler_ID_signatures.csv

id_num_data = pd.read_csv(id_num_file_path)
id_num_data.sort_values(by='Cancer Types', axis=0, inplace=True)
id_num_data.set_index(keys=['Cancer Types'], drop=False,inplace=True)

id_cancer_type = id_num_data['Cancer Types'].unique().tolist()
id_prop_data = pd.read_csv(id_prop_file_path)

#%%

"""
Find the proportions of different combinations of signatures among all samples
Input cancer_type_name as the name of cancer type to be simulated
Input input_data as *_num_data 
Output proportion of samples with each combination of signatures. 
"""

def sbs_sig_proportion(cancer_type, num_data):
    global sbs_prop
    temp_data = num_data.loc[num_data.index==cancer_type]
    temp_data2 = temp_data.copy().loc[:, 'Accuracy':]
    
    del temp_data2['Accuracy']
    cols = temp_data2.columns
    bt = temp_data2.apply(lambda x: x > 0)
    sbs_prop = (bt.apply(lambda x: list(cols[x.values]), axis=1).value_counts())/len(temp_data2)
    return sbs_prop
 
    

"""
DBS sig proportion
"""
def dbs_sig_proportion(cancer_type, num_data):
    global dbs_prop
    temp_data = num_data.loc[num_data.index==cancer_type]
    temp_data2 = temp_data.copy().loc[:, 'Accuracy':]
    
    del temp_data2['Accuracy']
    cols = temp_data2.columns
    bt = temp_data2.apply(lambda x: x > 0)
    dbs_prop = (bt.apply(lambda x: list(cols[x.values]), axis=1).value_counts())/len(temp_data2)
    return dbs_prop



"""
ID sig proportion
"""
def id_sig_proportion(cancer_type, num_data):
    global id_prop
    temp_data = num_data.loc[num_data.index==cancer_type]
    temp_data2 = temp_data.copy().loc[:, 'Accuracy':]
    
    del temp_data2['Accuracy']
    cols = temp_data2.columns
    bt = temp_data2.apply(lambda x: x > 0)
    id_prop = (bt.apply(lambda x: list(cols[x.values]), axis=1).value_counts())/len(temp_data2)
    return id_prop



"""
Function that runs sbs_sig_proportion, dbs_sig_proportion, id_sig_proportion
"""
cancer_type = input('Cancer type: ')
sbs_sig_proportion(cancer_type, sbs_num_data)
dbs_sig_proportion(cancer_type, dbs_num_data)
id_sig_proportion(cancer_type, id_num_data)

#%%

"""
Run once (already run).
Input *_num_data.
Output one file for each sample. File contains expected number of 96 type mutations for 
each type of signature in the sample. 
"""
def sbs_sig_freq(input_data):
    temp_data3 = input_data.copy().loc[:, 'Accuracy':]
    del temp_data3['Accuracy']
    cols = temp_data3.columns
    for row in range(len(temp_data3)):
        new_df = sbs_prop_data.loc[:,'Type':'SubType']
        for column in range(len(cols)):
            frequency = temp_data3.iloc[row, column]
            new_df[cols[column]] = [(frequency*item) for item in sbs_prop_data[cols[column]]]
        ct = temp_data3.index[row]
        sample_id = sbs_num_data.iloc[row, 1]
        with open('/Users/DavidChen/Desktop/Project/' + ct + '_' + sample_id + '.csv', 'w'):
            new_df.to_csv('/Users/DavidChen/Desktop/Project/' + ct + '_' + sample_id + '.csv')



"""
DBS version of sig_freq
"""
def dbs_sig_freq(input_data):
    temp_data3 = input_data.loc[:, 'Accuracy':]
    del temp_data3['Accuracy']
    cols = temp_data3.columns
    for row in range(len(temp_data3)):
        new_df = dbs_prop_data.loc[:,'Mutation Type'].to_frame() 
        for column in range(len(cols)):
            frequency = temp_data3.iloc[row, column]
            new_df[cols[column]] = [(frequency*item) for item in dbs_prop_data[cols[column]]]         
        ct = temp_data3.index[row]
        sample_id = dbs_num_data.iloc[row, 1]   
        with open('/Users/DavidChen/Desktop/Project/' + ct + '_' + sample_id + '.csv', 'w'):
            new_df.to_csv('/Users/DavidChen/Desktop/Project/' + ct + '_' + sample_id + '.csv')

"""
ID version of sig_freq
"""
def id_sig_freq(input_data):
    temp_data3 = input_data.loc[:, 'Accuracy':]
    del temp_data3['Accuracy']
    cols = temp_data3.columns
    for row in range(len(temp_data3)):
        new_df = id_prop_data.loc[:,'Mutation Type':'ID_type']
        for column in range(len(cols)):
            frequency = temp_data3.iloc[row, column]
            new_df[cols[column]] = [(frequency*item) for item in id_prop_data[cols[column]]]         
        ct = temp_data3.index[row]
        sample_id = id_num_data.iloc[row, 1]   
        with open('/Users/DavidChen/Desktop/Project/ID_' + ct + '_' + sample_id + '.csv', 'w'):
            new_df.to_csv('/Users/DavidChen/Desktop/Project/ID_' + ct + '_' + sample_id + '.csv')



#%%



"""
Calculate average number of 96 type mutations among all samples for specified cancer type.

Input the specified cancer type. 
Output the average number of 96 type mutations among all samples for specified cancer type.
"""

"""
SBS
"""
def sbs_mean(cancer_type, sbs_freq_folder_path):
    global sbs_average_freq_df
    global df_averages
    file_paths = glob.glob(str(sbs_freq_folder_path) + cancer_type + '_*' + '.csv')
    og_df = pd.read_csv(file_paths[0])
    mean_df = og_df.loc[:, 'SBS1':]
    for file_index in range(1, len(file_paths)):
        raw_mean = pd.read_csv(file_paths[file_index])
        slice_read_df = raw_mean.loc[:, 'SBS1':]
        mean_df += slice_read_df
    df_averages = mean_df/len(file_paths) 
    mut_subtypes = og_df.loc[:, 'Type': 'SubType']
    sbs_average_freq_df = pd.concat([mut_subtypes, df_averages], axis=1)
    return sbs_average_freq_df

"""
DBS
"""
def dbs_mean(cancer_type, dbs_freq_folder_path):
    global dbs_average_freq_df
    global dbs_df_averages
    file_paths = glob.glob(str(dbs_freq_folder_path) + cancer_type + '_*' + '.csv')
    
    og_df = pd.read_csv(file_paths[0], index_col=0)
    mean_df = og_df.loc[:, 'DBS1':]
    for file_index in range(1, len(file_paths)):
        raw_mean = pd.read_csv(file_paths[file_index])
        slice_read_df = raw_mean.loc[:, 'DBS1':]
        mean_df += slice_read_df
    dbs_df_averages = mean_df/len(file_paths) 
    mut_subtypes = og_df.loc[:,'Mutation Type'].to_frame() 
    dbs_average_freq_df = pd.concat([mut_subtypes, dbs_df_averages], axis=1)
    return dbs_average_freq_df


"""
ID
Input cancer type.
Output the average number of 72 type mutations among all samples within one cancer type.
"""
def id_mean(cancer_type, id_freq_folder_path):
    global id_average_freq_df
    global id_df_averages
    file_paths = glob.glob(str(id_freq_folder_path) + cancer_type + '_*' + '.csv')
    og_df = pd.read_csv(file_paths[0], index_col=0)
    mean_df = og_df.loc[:, 'ID1':]
    for file_index in range(1, len(file_paths)):
        raw_mean = pd.read_csv(file_paths[file_index])
        slice_read_df = raw_mean.loc[:, 'ID1':]
        mean_df += slice_read_df
    id_df_averages = mean_df/len(file_paths) 
    mut_subtypes = og_df.loc[:,'Mutation Type':'ID_type']
    id_average_freq_df = pd.concat([mut_subtypes, id_df_averages], axis=1)
    return id_average_freq_df



"""
Function that runs sbs_mean, dbs_mean, id_mean.
"""
sbs_freq_folder_path = input('Input absolute folder path for SBS_Expected_Frequency: ')
#/Users/DavidChen/Desktop/Project/SBS_Expected_Frequency/
dbs_freq_folder_path = input('Input absolute folder path for DBS_Expected_Frequency: ')
#/Users/DavidChen/Desktop/Project/DBS_Expected_Frequency/
id_freq_folder_path = input('Input absolute folder path for ID_Expected_Frequency: ')
#/Users/DavidChen/Desktop/Project/ID_Expected_Frequency/ID_

sbs_mean(cancer_type, sbs_freq_folder_path)
dbs_mean(cancer_type, dbs_freq_folder_path)
id_mean(cancer_type, id_freq_folder_path)


#%%



"""
Input 96 type average_freq_df for one cancer type.
Output 192 mutation type classification table averaged for one cancer type. 
Frequency of pyrmidine single base substitution duplicated for their Watson-Crick complementary base pair
"""
def sbs_pur_pyr_table(sbs_average_freq_df):
    global sbs_df
    raw_df = sbs_average_freq_df
    purines = ('G','A')
    pyrimidines = ('C','T')
        
    wc = {'G': 'C', 'C':'G', 'A':'T', 'T':'A'}
    non_wc = {'G': 'T', 'C':'A', 'A':'C', 'T':'G'}
    purines_bp = {'G': 'A', 'A': 'G'}
    pyrimidines_bp = {'C': 'T', 'T': 'C'}

    df_copy = raw_df.copy()
    sbs_df = raw_df.append(df_copy)
    sbs_df = sbs_df.reset_index(drop=True)
    
    for num in list(range(len(raw_df))):
        snp = sbs_df.loc[num]['Type']
        
        if snp[0] in purines and snp[2] in purines:
            x = wc.get(snp[0])
            sub = x + '>' + str(pyrimidines_bp.get(x))
            sbs_df.iat[(num + len(raw_df)), 0] = sub
           
            string = sbs_df.loc[num + len(raw_df)]['SubType']
            split_string = [character for character in string]
            split_string[1] = wc.get(snp[0])
            subtype = ''.join(split_string)
            sbs_df.iat[(num + len(raw_df)), 1] = subtype
             
        if snp[0] in pyrimidines and snp[2] in pyrimidines:
            x = wc.get(snp[0])
            sub = x + '>' + str(purines_bp.get(x))
            sbs_df.iat[(num + len(raw_df)), 0] = sub
            
            string = sbs_df.loc[num + len(raw_df)]['SubType']
            split_string = [character for character in string]
            split_string[1] = wc.get(snp[0])
            subtype = ''.join(split_string)
            sbs_df.iat[(num + len(raw_df)), 1] = subtype
            
        if wc.get(snp[0]) == snp[2]:
            sub = snp[2] + '>' + str(snp[0])
            sbs_df.iat[(num + len(raw_df)), 0] = sub
            
            string = sbs_df.loc[num + len(raw_df)]['SubType']
            split_string = [character for character in string]
            split_string[1] = wc.get(snp[0])
            subtype = ''.join(split_string)
            sbs_df.iat[(num + len(raw_df)), 1] = subtype
            
        if non_wc.get(snp[0]) == snp[2]:
            x = wc.get(snp[0]) 
            sub = x + '>' + str(non_wc.get(x))
            sbs_df.iat[(num + len(raw_df)), 0] = sub
            
            string = sbs_df.loc[num + len(raw_df)]['SubType']
            split_string = [character for character in string]
            split_string[1] = wc.get(snp[0])
            subtype = ''.join(split_string)
            sbs_df.iat[(num + len(raw_df)), 1] = subtype
            
        else:
            pass  
        
    return sbs_df  



"""
DBS pur_pyr_table
"""
def dbs_pur_pyr_table(dbs_average_freq_df):
    global dbs_df
    
    df_copy = dbs_average_freq_df.copy()
    raw_df = df_copy.append(df_copy)
    dbs_df = raw_df.reset_index(drop=True)
        
    opp_dict = {'G':'C', 'A':'T', 'C':'G','T':'A'}
    
    for row in list(range(len(dbs_df.iloc[78:,0]), 156)):
        mut_type = dbs_df.iloc[row, 0]
        
        for index in list(range(len(mut_type))):
            if mut_type[index] in list(opp_dict.keys()):
                mut_type = mut_type[:index] + mut_type[index].replace(mut_type[index], opp_dict[mut_type[index]]) + mut_type[index+1:]
            else:
                pass
        dbs_df.at[row, 'Mutation Type'] = mut_type
     


"""
ID pur_pyr_table
"""
def id_pur_pyr_table(id_average_freq_df):
    global id_df
    
    df_copy = id_average_freq_df.copy()
    raw_df = df_copy.loc[:23,:].append(df_copy)
    id_df = raw_df.reset_index(drop=True)
        
    opp_dict = {'G':'C', 'A':'T', 'C':'G','T':'A'}
    
    for row in list(range(24, 24+ len(id_df.iloc[:24,0]))):
        mut_type = id_df.iloc[row, 0]
        opp_mut_type = mut_type[:4] + mut_type[4].replace(mut_type[4], opp_dict[mut_type[4]]) + mut_type[5:]
        id_df.at[row, 'Mutation Type'] = opp_mut_type
        id_df.at[row,'Index'] = opp_dict[mut_type[4]]
     
        
sbs_pur_pyr_table(sbs_average_freq_df)        
dbs_pur_pyr_table(dbs_average_freq_df)    
id_pur_pyr_table(id_average_freq_df) 
       
#%% 
"""
Uses count_kmers(data,k) as a helper function.
No input.
Output the count of 3-mers as a dictionary format. 
Files already output in Folder named Ref_count.
"""

"""
SBS ref_count for 3_mer contexts.
"""
def sbs_ref_count():
    for seq_record in SeqIO.parse(input_file_path, "fasta"):
        name, seq = seq_record.id, seq_record.seq
        count_kmers(str(seq),3)
        count_df = pd.DataFrame.from_dict(d.keys())
        count_df['count'] = d.values()
        count_df.columns = ['trinucleotide', 'count']
        count_file_path= '/Users/DavidChen/Desktop/Ref_count/' + name + '.csv'
        with open(count_file_path, 'w'):
            count_df.to_csv(count_file_path)
      
"""
DBS ref_count for 2-mers.
"""
def dbs_ref_count(kmer):
    for seq_record in SeqIO.parse(input_file_path, "fasta"):
        name, seq = seq_record.id, seq_record.seq
        count_kmers(str(seq),kmer)
        count_df = pd.DataFrame.from_dict(d.keys())
        count_df['count'] = d.values()
        count_df.columns = ['dinucleotide', 'count']
        count_file_path= '/Users/DavidChen/Desktop/DBS_Ref_count/' + name + '.csv'
        with open(count_file_path, 'w'):
            count_df.to_csv(count_file_path)

"""
Ref_count for any input kmer length
"""
def kmer_ref_count(kmer_length):
    for seq_record in SeqIO.parse(input_file_path, "fasta"):
        name, seq = seq_record.id, seq_record.seq
        count_kmers(str(seq),kmer_length)
        count_df = pd.DataFrame.from_dict(d.keys())
        count_df['count'] = d.values()
        count_df.columns = [str(kmer_length), 'count']
        count_file_path= '/Users/DavidChen/Desktop/kmer_ref_count/'+ str(kmer_length)+'-mer_' + name + '.csv'
        with open(count_file_path, 'w'):
            count_df.to_csv(count_file_path)
     
"""
Repeater function for kmer counts for range of kmer_length = start:end(inclusive)
"""
def repeat_kmer_generator(kmer_length_start, kmer_length_end):
    for i in list(range(kmer_length_start, kmer_length_end+1)):
        kmer_ref_count(i)
        
#%%

"""
Add total counts for every kmer in Chromosome 1 to Chromosome 22, X and Y.
"""
def sbs_kmer_count():
    global sbs_genome_kmer_count
    chr_index = [i for i in range(2,23)]
    chr_index.append('X')
    chr_index.append('Y')
    
    chr_1 = pd.read_csv('/Users/DavidChen/Desktop/Ref_count/chr1.csv')
    total_kmer_count = pd.read_csv('/Users/DavidChen/Desktop/Ref_count/chr1.csv').loc[:, 'count']
    
    for file_name in chr_index:
        read_in = pd.read_csv('/Users/DavidChen/Desktop/Ref_count/chr' + str(file_name) + '.csv').loc[:, 'count']
        total_kmer_count += read_in
        
    sbs_genome_kmer_count = pd.concat([chr_1.loc[:,'trinucleotide'], total_kmer_count], axis=1)

"""
DBS gen_kmer_count of 2mers.
"""
def dbs_kmer_count():
    global dbs_genome_kmer_count
    chr_index = [i for i in range(2,23)]
    chr_index.append('X')
    chr_index.append('Y')
    
    chr_1 = pd.read_csv('/Users/DavidChen/Desktop/DBS_Ref_count/chr1.csv')
    total_kmer_count = pd.read_csv('/Users/DavidChen/Desktop/DBS_Ref_count/chr1.csv').loc[:, 'count']
    
    for file_name in chr_index:
        read_in = pd.read_csv('/Users/DavidChen/Desktop/DBS_Ref_count/chr' + str(file_name) + '.csv').loc[:, 'count']
        total_kmer_count += read_in
        
    dbs_genome_kmer_count = pd.concat([chr_1.loc[:,'Mutation Type'], total_kmer_count], axis=1)
    
"""
Count genome kmer counts (Chr1-22, X and Y)
"""
def genome_kmer_count(kmer_length, kmer_count_folder_path):
    chr_index = [i for i in range(kmer_length + 1 ,23)]
    chr_index.append('X')
    chr_index.append('Y')
    
    chr_1 = pd.read_csv(kmer_count_folder_path + '/' + str(kmer_length) + '-mer/' + str(kmer_length) + '-mer_chr1.csv')
    total_kmer_count = pd.read_csv(kmer_count_folder_path + '/' + str(kmer_length) + '-mer/' + str(kmer_length) + '-mer_chr1.csv').loc[:, 'count']
    
    for file_name in chr_index:
        read_in = pd.read_csv(kmer_count_folder_path + '/' + str(kmer_length) + '-mer/' + str(kmer_length) + '-mer_chr' + str(file_name) + '.csv').loc[:, 'count']
        total_kmer_count += read_in
        
    return pd.concat([chr_1.iloc[:,1], total_kmer_count], axis=1)
    

"""
Total Kmer counts for lengths 1-6 of GChr.38
"""
kmer_count_folder_path = input('Input absolute folder path for kmer_ref_count: ')
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

#%%


"""
No input.
Output dictionary of average probabilities for 192 mutation type.
"""
def sbs_mut_prob():
    global sbs_index_mut_prob
    sbs_mut = sbs_df.copy()
    sbs_mut['total'] = sbs_mut['SubType'].map(kmer_3_count.set_index('kmer')['count'])
    mut_prob = sbs_df.loc[:,'SBS1':'SBS60'].div(sbs_mut.total, axis=0)
    sbs_index_mut_prob = pd.concat([sbs_df[['Type', 'SubType']], mut_prob], axis=1)
    sbs_index_mut_prob.sort_values(['SubType','Type']).set_index(['SubType','Type'],inplace=True)
    
    
"""
DBS mut_prob
Input kmer_2_count
"""   
def dbs_mut_prob():
    global dbs_index_mut_prob
    dbs_mut = dbs_df.copy()
    dbs_mut['2mer_index'] = dbs_mut['Mutation Type'].apply(lambda x: x[:2])
    dbs_mut['total'] = dbs_mut['2mer_index'].map(kmer_2_count.set_index('kmer')['count'])
    dbs_mut_prob = dbs_mut.loc[:,'DBS1':'DBS11'].div(dbs_mut['total'], axis=0)
    dbs_index_mut_prob = pd.concat([dbs_mut[['Mutation Type', '2mer_index']], dbs_mut_prob], axis=1)
    dbs_index_mut_prob.sort_values(['Mutation Type']).set_index(['Mutation Type'],inplace=True)
   
    
"""
Indel 1 base, 1-6 repeats mut_prob 
"""
    
def id_mut_prob_1_6():
    global id_index_mut_prob_1_6
    id_mut = id_df.copy()
    id_mut['total'] = id_mut['Index'].map(kmer_1_count.set_index('kmer')['count'])
    
    id_mut_prob = id_mut.loc[:,'ID1':'ID17'].div(id_mut['total'], axis=0)
    id_index_mut_prob_1_6 = pd.concat([id_mut.loc[:, 'Mutation Type':'ID_type'], id_mut_prob], axis=1)
    id_index_mut_prob_1_6.sort_values(['Mutation Type']).set_index(['Mutation Type'],inplace=True)
      
    
   
sbs_mut_prob()
dbs_mut_prob()
id_mut_prob_1_6()  
    
#%%    
    
    

"""
ID mutation for indel 1 base, 1-6 repeats mut_prob
"""
def mutation(number_of_lineages, number_of_generations):
    
    
    for i in range(1, number_of_lineages+1):
        
        print('Lineage' + str(i))
    
        insertion_dict={}
        deletion_dict={}
        sbs_dict={}
        dbs_dict={}
        
        for j in range(1, number_of_generations+1):
            
            print ('Generation' + str(j))
            
            id_sig_combination = random.choices(list(id_prop.index), list(id_prop.values))
            id_prob = id_index_mut_prob_1_6.copy()
            id_mut_sig_df = pd.concat([id_prob.loc[:47,'Mutation Type':'ID_type'], id_prob[id_sig_combination[0]].dropna()], axis=1)
            
            id_total_prob_df = id_mut_sig_df.loc[:,:'ID1']
            del id_total_prob_df['ID1']
            id_total_prob_df['total_prob']= id_mut_sig_df.loc[:, 'ID_type':].sum(axis=1)
            id_sorted_df = id_total_prob_df.sort_values(['ID_type','Index','Size'])
            id_sorted_df.set_index(['Index'], inplace=True)
            
            deletion_df = id_sorted_df.iloc[:24, :]
            insertion_df = id_sorted_df.iloc[24:, :]
            
            sbs_sig_combination = random.choices(list(sbs_prop.index), list(sbs_prop.values))
            sbs_mut_sig_df = pd.concat([sbs_index_mut_prob.loc[:,'Type':'SubType'], sbs_index_mut_prob[sbs_sig_combination[0]]], axis=1)
            sbs_total_prob_df = sbs_mut_sig_df.iloc[:,:2]
            sbs_total_prob_df['total_prob']= sbs_mut_sig_df.iloc[:, 2:].sum(axis=1)
            sbs_sorted_df = sbs_total_prob_df.sort_values(['SubType'])
            sbs_sorted_df.set_index(['SubType'], inplace=True)
            
            dbs_sig_combination = random.choices(list(dbs_prop.index), list(dbs_prop.values))    
            dbs_mut_sig_df = pd.concat([dbs_index_mut_prob.loc[:,'Mutation Type':'2mer_index'], dbs_index_mut_prob[dbs_sig_combination[0]]], axis=1)
            dbs_total_prob_df = dbs_mut_sig_df.copy()    
            dbs_total_prob_df['total_prob']= dbs_mut_sig_df.iloc[:, 1:].sum(axis=1)  
            dbs_sorted_df = dbs_total_prob_df.sort_values(['Mutation Type'])
            dbs_sorted_df.set_index(['2mer_index'], inplace=True)
               
            index = list(range(len(sample)-6))
            random.shuffle(index)
            
            for k in index:
         
                window = str(sample[k:k+6])
            
                count=1
                
                for length in list(range(1, len(window))):
                    if window[length] == window[0]:
                        count += 1
                    else:
                        break
                    
                    
                    
                    
                i_type = [count] + [None]
                temp_i_prob = [insertion_df.loc[insertion_df['Size'] == count].loc[window[0],'total_prob']] 
                i_prob = temp_i_prob + [1-temp_i_prob[0]]
                
                i_mut = np.random.choice(i_type,p=i_prob) 
                if i_mut == None:
                    pass
                else:
                    insertion_dict[str(k)] = window[0]*i_mut
                        
     
                d_type = [count] + [None]
                temp_d_prob = [deletion_df.loc[deletion_df['Size'] == count].loc[window[0],'total_prob']] 
                d_prob = temp_d_prob + [1-temp_d_prob[0]]
                
                d_mut = np.random.choice(d_type,p=d_prob) 
                if d_mut == None:
                    pass
                else:
                    deletion_dict[str(k)] = window[0]*d_mut
                            
                sbs_types = list(sbs_sorted_df.loc[window[:3],'Type']) + [None]
                sbs_prob = list(sbs_sorted_df.loc[window[:3],'total_prob'])+ [1- sum(list(sbs_sorted_df.loc[window[:3],'total_prob']))]
                sbs = np.random.choice(sbs_types,p=sbs_prob) 
                if sbs == None:
                    pass
                else:
                    sbs_dict[str(k)] = sbs[2]
                    
                    
                if window[:2] in dbs_sorted_df.index:
                    dbs_types = list(dbs_sorted_df.loc[window[:2],'Mutation Type']) + [None]     
                    dbs_prob = list(dbs_sorted_df.loc[window[:2],'total_prob'])+ [1- sum(list(dbs_sorted_df.loc[window[:2],'total_prob']))]
                    dbs = np.random.choice(dbs_types,p=dbs_prob)
                    if dbs == None:
                        pass
                    else:
                        dbs_dict[str(k)] = dbs[3:]
                else:
                    pass
             
            print(insertion_dict)
            print(deletion_dict)
            print(sbs_dict)
            print(dbs_dict)
                

#%%            
from collections import defaultdict           

  
def sample_indexer(sample):
    global sample_index_dict
    global sample_count_dict
    
    sample_index_dict = defaultdict(list)

    for seq_len in range(1, len(sample)-5):
        if 'N' in sample[seq_len:seq_len+6]:
            pass
        
        else:
            sample_index_dict[str(sample[seq_len:seq_len+6])].append(seq_len)
        
    sample_count_dict = {key: len(value) for key, value in sample_index_dict.items()}

sample_indexer(sample)


#%%
    
def somatic_sim(number_of_lineages, number_of_generations, power=1):
     
    
    #Setup probabilities into processed dataframes
    for i in range(1, number_of_lineages+1):
        
        print('Lineage' + str(i))
    
        insertion_dict={}
        deletion_dict={}
        sbs_dict={}
        dbs_dict={}
        
        id_sig_combination = random.choices(list(id_prop.index), list(id_prop.values))
        id_prob = id_index_mut_prob_1_6.copy()
        id_mut_sig_df = pd.concat([id_prob.loc[:47,'Mutation Type':'ID_type'], id_prob[id_sig_combination[0]].dropna()], axis=1)
            
        id_total_prob_df = id_mut_sig_df.iloc[:,:5]
        id_total_prob_df['total_prob']= id_mut_sig_df.loc[:, 'ID_type':].sum(axis=1)
            
        id_total_prob_df['total_prob'] = id_total_prob_df['total_prob'].multiply(power)
           
        id_sorted_df = id_total_prob_df.sort_values(['ID_type','Index','Size'])
        id_sorted_df.set_index(['Index'], inplace=True)
            
        deletion_df = id_sorted_df.iloc[:24, :]
        insertion_df = id_sorted_df.iloc[24:, :]
            
            
            
        sbs_sig_combination = random.choices(list(sbs_prop.index), list(sbs_prop.values))
        sbs_mut_sig_df = pd.concat([sbs_index_mut_prob.loc[:,'Type':'SubType'], sbs_index_mut_prob[sbs_sig_combination[0]]], axis=1)
        sbs_total_prob_df = sbs_mut_sig_df.iloc[:,:2]
        sbs_total_prob_df['total_prob'] = sbs_mut_sig_df.iloc[:, 2:].sum(axis=1)
            
        sbs_total_prob_df['total_prob'] = sbs_total_prob_df[sbs_total_prob_df.select_dtypes(include=['number']).columns] * int(power)
        
        sbs_sorted_df = sbs_total_prob_df.sort_values(['SubType'])
        sbs_sorted_df.set_index(['SubType'], inplace=True)
            
            
            
        dbs_sig_combination = random.choices(list(dbs_prop.index), list(dbs_prop.values))    
        dbs_mut_sig_df = pd.concat([dbs_index_mut_prob.loc[:,'Mutation Type':'2mer_index'], dbs_index_mut_prob[dbs_sig_combination[0]]], axis=1)
        dbs_total_prob_df = dbs_mut_sig_df.copy()    
        dbs_total_prob_df['total_prob']= dbs_mut_sig_df.iloc[:, 1:].sum(axis=1) 
            
        dbs_total_prob_df['total_prob'] = dbs_total_prob_df[dbs_total_prob_df.select_dtypes(include=['number']).columns] * power

        dbs_sorted_df =  dbs_total_prob_df.sort_values(['Mutation Type'])
        dbs_sorted_df.set_index(['2mer_index'], inplace=True)   
        
        sample_seq = sample
        
        for j in range(1, number_of_generations+1):
            
            print('Generation' + str(j))
            
            #Index of each kmer
            
            sample_index_dict = defaultdict(list)

            for seq_len in range(1, len(sample_seq)-5):
                if 'N' in sample_seq[seq_len:seq_len+6]:
                    pass
        
                else:
                    sample_index_dict[str(sample_seq[seq_len:seq_len+6])].append(seq_len)
        
            #Count of each kmer
            sample_count_dict = {key: len(value) for key, value in sample_index_dict.items()}
        
           
            for keys in list(sample_index_dict.keys()):
                
                count=1
                for length in list(range(1, len(keys))):
                    if keys[length] == keys[0]:
                        count += 1
                    else:
                        break
                
                #Insertions
                temp_i_prob = [insertion_df.loc[insertion_df['Size'] == count].loc[keys[0],'total_prob']] 
                insertion_count = np.random.binomial(sample_count_dict[keys], temp_i_prob, 1) 
                
                for m in range(sum(i_number > 0 for i_number in list(insertion_count))):
                    insertion_index = np.random.choice(sample_index_dict[keys])
                    insertion_dict[str(insertion_index)] = keys[0]
                    
                    
                #Deletion
                temp_d_prob = [deletion_df.loc[deletion_df['Size'] == count].loc[keys[0],'total_prob']] 
                deletion_count = np.random.binomial(sample_count_dict[keys], temp_d_prob, 1)  
                
                for n in range(sum(d_number > 0 for d_number in list(deletion_count))):
                    deletion_index = np.random.choice(sample_index_dict[keys])
                    deletion_dict[str(deletion_index)] = keys[0]
            
                #SBS
                sbs_types = list(sbs_sorted_df.loc[keys[:3],'Type']) + [None]
                sbs_prob = list(sbs_sorted_df.loc[keys[:3],'total_prob'])+ [1- sum(list(sbs_sorted_df.loc[keys[:3],'total_prob']))]             
                sbs = [sbs_num for sbs_num in list(np.random.choice(sbs_types,p=sbs_prob, size = sample_count_dict[keys])) if sbs_num]
                for sbs_value in sbs:
                    single_base_index = np.random.choice(sample_index_dict[keys])
                    sbs_dict[single_base_index+1] = sbs_value[2]
    
                      
                #DBS
                if keys[:2] in dbs_sorted_df.index:
                    dbs_types = list(dbs_sorted_df.loc[keys[:2],'Mutation Type']) + [None]     
                    dbs_prob = list(dbs_sorted_df.loc[keys[:2],'total_prob'])+ [1- sum(list(dbs_sorted_df.loc[keys[:2],'total_prob']))]
                    dbs = [dbs_num for dbs_num in list(np.random.choice(dbs_types,p=dbs_prob, size = sample_count_dict[keys])) if dbs_num]
                    for dbs_value in dbs:
                        double_base_index = np.random.choice(sample_index_dict[keys])
                        dbs_dict[double_base_index] = dbs_value[3:]
                else:
                    pass
                
            
            for index_insertion in list(insertion_dict.keys()):
                sample_seq = sample_seq[:int(index_insertion)] + insertion_dict[index_insertion] + sample_seq[int(index_insertion):]
                
            for index_deletion in list(deletion_dict.keys()):
                sample_seq = sample_seq[:int(index_deletion)] + sample_seq[int(index_deletion)+1:]
                
            for index_sbs in list(sbs_dict.keys()):
                sample_seq = sample_seq[:index_sbs] +  sbs_dict[index_sbs] + sample_seq[index_sbs+1:]
                
            for index_dbs in list(dbs_dict.keys()):
                sample_seq = sample_seq[:index_dbs] +  dbs_dict[index_dbs] + sample_seq[index_dbs+2:]
                
            print(insertion_dict)
            print(deletion_dict)
            print(sbs_dict)
            print(dbs_dict)
            
  
            sample_sequence_file = '/Users/DavidChen/Desktop/Sample/' + cancer_type + '_Generation_' + str(power * j) + '_Lineage_' + str(i) + '.fasta'
            with open(sample_sequence_file, 'w+') as fasta_file:
                fasta_file.write(sample_seq)
                








#%%
def boyer_moore_horspool(pattern, text):
    m = len(pattern)
    n = len(text)

    if m > n:
        return -1

    skip = defaultdict(lambda: m)
    found_indexes = []

    for k in range(m - 1):
        skip[ord(pattern[k])] = m - k - 1

    k = m - 1

    while k < n:
        j = m - 1
        i = k
        while j >= 0 and text[i] == pattern[j]:
            j -= 1
            i -= 1
        if j == -1:
            found_indexes.append(i + 1)

        k += skip[ord(text[k])]

    return found_indexes
#%%















#
#            func1 = functools.partial(mutation)
#            func2 = functools.partial(id_mutation)
#            pool = Pool(processes=4)
#            res = pool.map(foo, [func1, func2])
#            pool.close()
#            pool.join()
#            print(res)
#           
#
#
#
#
#
#
#          
#            sbs_sig_combination = random.choices(list(sbs_prop.index), list(sbs_prop.values))
#            sbs_mut_sig_df = pd.concat([sbs_index_mut_prob.loc[:,'Type':'SubType'], sbs_index_mut_prob[sbs_sig_combination[0]]], axis=1)
#            sbs_total_prob_df = sbs_mut_sig_df.iloc[:,:2]
#            sbs_total_prob_df['total_prob']= sbs_mut_sig_df.iloc[:, 2:].sum(axis=1)
#            sbs_sorted_df = sbs_total_prob_df.sort_values(['SubType'])
#            sbs_sorted_df.set_index(['SubType'], inplace=True)
#            
#            dbs_sig_combination = random.choices(list(dbs_prop.index), list(dbs_prop.values))    
#            dbs_mut_sig_df = pd.concat([dbs_index_mut_prob.loc[:,'Mutation Type':'2mer_index'], dbs_index_mut_prob[dbs_sig_combination[0]]], axis=1)
#            dbs_total_prob_df = dbs_mut_sig_df.copy()    
#            dbs_total_prob_df['total_prob']= dbs_mut_sig_df.iloc[:, 1:].sum(axis=1)  
#            dbs_sorted_df = dbs_total_prob_df.sort_values(['Mutation Type'])
#            dbs_sorted_df.set_index(['2mer_index'], inplace=True)
#            
#            index = list(range(len(sample)-2))
#            random.shuffle(index)
#            
#            
#            for l in index:
#                window = str(sample[i:i+3])
#                    
#                sbs_types = list(sbs_sorted_df.loc[window,'Type']) + [None]
#                sbs_prob = list(sbs_sorted_df.loc[window,'total_prob'])+ [1- sum(list(sbs_sorted_df.loc[window,'total_prob']))]
#                sbs = np.random.choice(sbs_types,p=sbs_prob) 
#                if sbs == None:
#                    pass
#                else:
#                    sbs_dict[str(i)] = sbs[2]
#                    
#                if window in dbs_sorted_df.index:
#                    dbs_types = list(dbs_sorted_df.loc[window[:2],'Mutation Type']) + [None]     
#                    dbs_prob = list(dbs_sorted_df.loc[window[:2],'total_prob'])+ [1- sum(list(dbs_sorted_df.loc[window[:2],'total_prob']))]
#                    dbs = np.random.choice(dbs_types,p=dbs_prob)
#                    if dbs == None:
#                        pass
#                    else:
#                        dbs_dict[str(i)] = dbs[3:]
#                else:
#                    pass
#             
#                
#                
#                
#                
#                
#            for index_dicts in [insertion_dict, deletion_dict, sbs_dict, dbs_dict]:
#                for n in list(index_dicts.keys()):
#                    sample_seq = sample_seq[:n] + sbs_dict[n] + sample_seq[n+1:]
#            
#            
#            
#                
#            print(insertion_dict)
#            print(deletion_dict)
#            print(sbs_dict)
#            print(dbs_dict)

    
    

  




  
"""
Run the following for mutation simulation.
"""




def somatic_sim():
    """
    Input sequence file path
    """
    input_file_path=input('Input file path of sequence: ')
    #/Users/DavidChen/Desktop/NC_000015.10[20150295..20200171].fasta
    
    """
    SBS input reference files
    """
    sbs_num_file_path = input('Input file path of PCAWG_sigProfiler_SBS_signatures_in_samples.csv: ')
    #/Users/DavidChen/Desktop/Project/Reference/PCAWG_sigProfiler_SBS_signatures_in_samples.csv
    sbs_prop_file_path = input('Input file path of sigProfiler_SBS_signatures.csv: ')
    #/Users/DavidChen/Desktop/Project/Reference/sigProfiler_SBS_signatures.csv
    
    sbs_num_data = pd.read_csv(sbs_num_file_path)
    sbs_num_data.sort_values(by='Cancer Types', axis=0, inplace=True)
    sbs_num_data.set_index(keys=['Cancer Types'], drop=False,inplace=True)
    
    sbs_cancer_type = sbs_num_data['Cancer Types'].unique().tolist()
    sbs_prop_data = pd.read_csv(sbs_prop_file_path)
      
    """
    DBS input reference files
    """
    dbs_num_file_path = input('Input file path of PCAWG_sigProfiler_DBS_signatures_in_samples.csv: ')
    #/Users/DavidChen/Desktop/Project/Reference/PCAWG_sigProfiler_DBS_signatures_in_samples.csv
    dbs_prop_file_path = input('Input file path of sigProfiler_DBS_signatures.csv: ')
    #/Users/DavidChen/Desktop/Project/Reference/sigProfiler_DBS_signatures.csv
    
    dbs_num_data = pd.read_csv(dbs_num_file_path)
    dbs_num_data.sort_values(by='Cancer Types', axis=0, inplace=True)
    dbs_num_data.set_index(keys=['Cancer Types'], drop=False,inplace=True)
    
    dbs_cancer_type = dbs_num_data['Cancer Types'].unique().tolist()
    dbs_prop_data = pd.read_csv(dbs_prop_file_path)
    
    """
    Insertion/Deletion input reference files
    """
    id_num_file_path = input('Input file path of PCAWG_sigProfiler_ID_signatures_in_samples.csv: ')
    #/Users/DavidChen/Desktop/Project/Reference/PCAWG_SigProfiler_ID_signatures_in_samples.csv
    id_prop_file_path = input('Input file path of sigProfiler_ID_signatures.csv: ')
    #/Users/DavidChen/Desktop/Project/Reference/sigProfiler_ID_signatures.csv
    
    id_num_data = pd.read_csv(id_num_file_path)
    id_num_data.sort_values(by='Cancer Types', axis=0, inplace=True)
    id_num_data.set_index(keys=['Cancer Types'], drop=False,inplace=True)
    
    id_cancer_type = id_num_data['Cancer Types'].unique().tolist()
    id_prop_data = pd.read_csv(id_prop_file_path)
                
    
    
    
    
    
    cancer_type = input('Cancer type: ')
    sbs_sig_proportion(cancer_type, sbs_num_data)
    dbs_sig_proportion(cancer_type, dbs_num_data)
    id_sig_proportion(cancer_type, id_num_data)
    
    sbs_freq_folder_path = input('Input absolute folder path for SBS_Expected_Frequency: ')
    #/Users/DavidChen/Desktop/Project/SBS_Expected_Frequency/
    dbs_freq_folder_path = input('Input absolute folder path for DBS_Expected_Frequency: ')
    #/Users/DavidChen/Desktop/Project/DBS_Expected_Frequency/
    id_freq_folder_path = input('Input absolute folder path for ID_Expected_Frequency: ')
    #/Users/DavidChen/Desktop/Project/ID_Expected_Frequency/ID_
    
    sbs_mean(cancer_type, sbs_freq_folder_path)
    dbs_mean(cancer_type, dbs_freq_folder_path)
    id_mean(cancer_type, id_freq_folder_path)
    
    sbs_pur_pyr_table(sbs_average_freq_df)
    dbs_pur_pyr_table(dbs_average_freq_df)
    id_pur_pyr_table(id_average_freq_df)
    
    kmer_count_folder_path = input('Input absolute folder path for kmer_ref_count: ')
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
    
    sbs_mut_prob()
    dbs_mut_prob()
    id_mut_prob_1_6()
    

    
    def foo(f):
        return f()
    func1 = functools.partial(mutation)
    func2 = functools.partial(id_mutation)
    pool = Pool(processes=4)
    res = pool.map(foo, [func1, func2])
    pool.close()
    pool.join()
    print(res)





