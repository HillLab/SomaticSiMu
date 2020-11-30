import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as patches
import numpy as np
import os
import sys
import argparse

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

def mut_catalog(cancer_type, simulation_type, gen_start, gen_end, mut_type, output_path):
    try:
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
            graph_data.fillna(0, inplace=True)
            
            if len(range(gen_start, gen_end+1)) > 1:
                
                upper_error = [graph_data.iloc[mut,:-2].quantile(0.95) for mut in range(96)]
                per_95 = []
                for i in upper_error:
                    per_95.append(i - graph_data.iloc[upper_error.index(i), :-2].mean(axis=0))
                per_95 = [100*x for x in per_95]
                
                lower_error = [graph_data.iloc[mut,:-2].quantile(0.05) for mut in range(96)]
                per_5 = []
                counter = 0
                for i in lower_error:
                    if i == 0:
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0))
                    else:  
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0) - i)
                    counter += 1
                per_5 = [100*x for x in per_5]
            
                fig = plt.figure(figsize=(50, 10))
                ax = fig.add_subplot(111)
                
                N = 96
                ind = np.arange(96)   
                
                color_list = [(0.416, 0.733, 0.918), (0,0,0), (0.765, 0.172, 0.157), (0.785, 0.785, 0.785), (0.678, 0.808, 0.412), (0.878, 0.773, 0.769)]
                
                p1 = ax.bar(range(0,16), graph_data.iloc[0:16, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[0], yerr = (per_5[0:16], per_95[0:16]),capsize=5)
                p2 = ax.bar(range(16,32), graph_data.iloc[16:32, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[1], yerr = (per_5[16:32], per_95[16:32]),capsize=5)
                p3 = ax.bar(range(32,48), graph_data.iloc[32:48, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[2], yerr = (per_5[32:48], per_95[32:48]),capsize=5)
                p4 = ax.bar(range(48,64), graph_data.iloc[48:64, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[3], yerr = (per_5[48:64], per_95[48:64]),capsize=5)
                p5 = ax.bar(range(64,80), graph_data.iloc[64:80, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[4], yerr = (per_5[64:80], per_95[64:80]),capsize=5)
                p6 = ax.bar(range(80,96), graph_data.iloc[80:96, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[5], yerr = (per_5[80:96], per_95[80:96]),capsize=5)
                  
            else:
            
                fig = plt.figure(figsize=(50, 10))
                ax = fig.add_subplot(111)
                
                N = 96
                ind = np.arange(96)   
                
                color_list = [(0.416, 0.733, 0.918), (0,0,0), (0.765, 0.172, 0.157), (0.785, 0.785, 0.785), (0.678, 0.808, 0.412), (0.878, 0.773, 0.769)]
                
                p1 = ax.bar(range(0,16), graph_data.iloc[0:16, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[0])
                p2 = ax.bar(range(16,32), graph_data.iloc[16:32, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[1])
                p3 = ax.bar(range(32,48), graph_data.iloc[32:48, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[2])
                p4 = ax.bar(range(48,64), graph_data.iloc[48:64, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[3])
                p5 = ax.bar(range(64,80), graph_data.iloc[64:80, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[4])
                p6 = ax.bar(range(80,96), graph_data.iloc[80:96, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[5])
                
                
                
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
            fig.savefig(output_path, bbox_inches='tight')
            
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
            graph_data.fillna(0, inplace=True)
            
            for gen in range(gen_start, gen_end+1):
                 graph_data[gen] = graph_data[gen].div(graph_data[gen].sum())
            
            graph_data.sort_values(by=["Mutation Type", "2mer_index"],inplace=True)
            graph_data.fillna(0, inplace=True)
            

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
                per_95 =  [0 if i < 0 else i*100 for i in per_95]
                       
                lower_error = [graph_data.iloc[mut,:-2].quantile(0.05)*100 for mut in range(78)]
                per_5 = []
                counter = 0
                for i in lower_error:
                    if i == 0:
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0)*100)
                    else:  
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0)*100 - i)
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
                
                p1 = ax.bar(range(0,9), graph_data.iloc[0:9, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[0], yerr = (per_5[0:9], per_95[0:9]),capsize=5)
                p2 = ax.bar(range(9,15), graph_data.iloc[9:15, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[1], yerr = (per_5[9:15], per_95[9:15]),capsize=5)
                p3 = ax.bar(range(15,24), graph_data.iloc[15:24, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[2], yerr = (per_5[15:24], per_95[15:24]),capsize=5)
                p4 = ax.bar(range(24,30), graph_data.iloc[24:30, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[3], yerr = (per_5[24:30], per_95[24:30]),capsize=5)
                p5 = ax.bar(range(30,39), graph_data.iloc[30:39, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[4], yerr = (per_5[30:39], per_95[30:39]),capsize=5)
                p6 = ax.bar(range(39,45), graph_data.iloc[39:45, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[5], yerr = (per_5[39:45], per_95[39:45]),capsize=5)
                p7 = ax.bar(range(45,51), graph_data.iloc[45:51, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[6], yerr = (per_5[45:51], per_95[45:51]),capsize=5)
                p8 = ax.bar(range(51,60), graph_data.iloc[51:60, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[7], yerr = (per_5[51:60], per_95[51:60]),capsize=5)
                p9 = ax.bar(range(60,69), graph_data.iloc[60:69, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[8], yerr = (per_5[60:69], per_95[60:69]),capsize=5)
                p10 = ax.bar(range(69,78), graph_data.iloc[69:78, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[9], yerr = (per_5[69:78], per_95[69:78]),capsize=5)
                
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
                
                p1 = ax.bar(range(0,9), graph_data.iloc[0:9, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[0])
                p2 = ax.bar(range(9,15), graph_data.iloc[9:15, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[1])
                p3 = ax.bar(range(15,24), graph_data.iloc[15:24, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[2])
                p4 = ax.bar(range(24,30), graph_data.iloc[24:30, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[3])
                p5 = ax.bar(range(30,39), graph_data.iloc[30:39, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[4])
                p6 = ax.bar(range(39,45), graph_data.iloc[39:45, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[5])
                p7 = ax.bar(range(45,51), graph_data.iloc[45:51, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[6])
                p8 = ax.bar(range(51,60), graph_data.iloc[51:60, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[7])
                p9 = ax.bar(range(60,69), graph_data.iloc[60:69, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[8])
                p10 = ax.bar(range(69,78), graph_data.iloc[69:78, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[9])

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
            fig.savefig(output_path, bbox_inches='tight')
            
        
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
            graph_data.fillna(0, inplace=True)
            

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
                per_95 =  [0 if i < 0 else i*100 for i in per_95]
                       
                lower_error = [graph_data.iloc[mut,:-2].quantile(0.05) for mut in range(12)]
                per_5 = []
                counter = 0
                for i in lower_error:
                    if i == 0:
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0)*100)
                    else:  
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0)*100 - i)
                    counter += 1

                fig = plt.figure(figsize=(50, 10))
                ax = fig.add_subplot(111)
                
                ind = np.arange(12)   
                
                color_list = [(0.949, 0.753, 0.478), 
                              (0.937, 0.512, 0.2)]
                
                p1 = ax.bar(range(0,6), graph_data.iloc[0:6, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[0],  yerr = (per_5[0:6], per_95[0:6]),capsize=5)
                p2 = ax.bar(range(6,12), graph_data.iloc[6:12, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[1],  yerr = (per_5[6:12], per_95[6:12]),capsize=5)
                
            else:
                fig = plt.figure(figsize=(50, 10))
                ax = fig.add_subplot(111)
                
                ind = np.arange(12)   
                
                color_list = [(0.949, 0.753, 0.478), 
                              (0.937, 0.512, 0.2)]
                
                p1 = ax.bar(range(0,6), graph_data.iloc[0:6, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[0])
                p2 = ax.bar(range(6,12), graph_data.iloc[6:12, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[1])
                
               
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
            fig.savefig(output_path, bbox_inches='tight')
            
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
            graph_data.fillna(0, inplace=True)
            

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
                per_95 =  [0 if i < 0 else i*100 for i in per_95]
                       
                lower_error = [graph_data.iloc[mut,:-2].quantile(0.05) for mut in range(12)]
                per_5 = []
                counter = 0
                for i in lower_error:
                    if i == 0:
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0)*100)
                    else:  
                        per_5.append(graph_data.iloc[counter, :-2].mean(axis=0)*100 - i)
                    counter += 1
            
                fig = plt.figure(figsize=(50, 10))
                ax = fig.add_subplot(111)
            
                ind = np.arange(12)   
            
                color_list = [(0.727, 0.855, 0.569), 
                          (0.345, 0.612, 0.247)]
            
                p1 = ax.bar(range(0,6), graph_data.iloc[0:6, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[0],  yerr = (per_5[0:6], per_95[0:6]),capsize=5)
                p2 = ax.bar(range(6,12), graph_data.iloc[6:12, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[1],  yerr = (per_5[6:12], per_95[6:12]),capsize=5)
                
            else:
                fig = plt.figure(figsize=(50, 10))
                ax = fig.add_subplot(111)
            
                ind = np.arange(12)   
            
                color_list = [(0.727, 0.855, 0.569), 
                          (0.345, 0.612, 0.247)]
            
                p1 = ax.bar(range(0,6), graph_data.iloc[0:6, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[0])
                p2 = ax.bar(range(6,12), graph_data.iloc[6:12, :-2].mean(axis=1).multiply(100), bottom=0, color = color_list[1])
                
               
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
            fig.savefig(output_path, bbox_inches='tight')


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
            fig, ax = plt.subplots(1, 2, figsize=(20, 10), sharey=True, gridspec_kw={'width_ratios': [10, 1]})
            fig.tight_layout()
            X = np.arange(gen_start, gen_end+1)
            
            ax[0].bar(X, mutation_burden['sbs'], color="#E74C3C")
            ax[0].bar(X, mutation_burden['dbs'], bottom=mutation_burden['sbs'], color = "#3498DB")
            ax[0].bar(X, mutation_burden['ins'], bottom=mutation_burden['dbs'] + mutation_burden['sbs'], color="#27AE60")
            ax[0].bar(X, mutation_burden['del'], bottom=mutation_burden['dbs'] + mutation_burden['sbs'] + mutation_burden['ins'], color="#9B59B6")
            ax[1].boxplot(mutation_burden.sum(axis=1))
            
            #ax1 = mutation_burden.plot.bar(stacked=True, color= ["#E74C3C", "#3498DB","#27AE60", "#9B59B6"], figsize=(10,7))
            ax[0].set_ylim(mutation_burden.min().sum() - 0.2*mutation_burden.max().sum() - (mutation_burden.min().sum() - 0.1*mutation_burden.max().sum())%10, mutation_burden.max().sum() + 0.15*mutation_burden.max().sum()+10 - (mutation_burden.max().sum() + 0.2*mutation_burden.max().sum())%10)
            ax[0].set_xlim(X[0]-0.5, X[-1]+0.5)
            
            ax[0].tick_params(axis='x', rotation=0)
            ax[0].set_xlabel("Sequence ID")
            ax[0].set_ylabel("Number of simulated mutations")
            ax[0].xaxis.labelpad= 10
            ax[0].yaxis.labelpad= 10
            ax[1].set_xlabel("Boxplot of Total Simulated \n Mutation Burden")
            
            colors = {'SBS':"#E74C3C", 'DBS':"#3498DB", 'INS':"#27AE60", 'DEL':"#9B59B6"}         
            labels = list(colors.keys())
            handles = [plt.Rectangle((0,0),1,1, color=colors[label]) for label in labels]
            ax[0].legend(handles, labels)
            ax[0].set_xticks(list(X))
            ax[1].set_xticklabels([])
            fig.savefig(output_path, bbox_inches='tight')
    except:
        print('Error')
        
    
        
        
def main(): 
    
    # Initiate the parser
    parser = argparse.ArgumentParser(description="SomaticSiMu Parameters")
    
    # Add long and short argument
    parser.add_argument("--cancer", "-c", help="cancer type")
    parser.add_argument("--simulation", "-s", help="[temporal stage, default is End]", default="End")
    parser.add_argument("--gen_start", "-a", help="start of range of sequence metadata")
    parser.add_argument("--gen_end", "-b", help="end of range of sequence metadata")
    parser.add_argument("--mut_type", "-m", help="mutation type: sbs, dbs, insertion, deletion, mutation_burden")
    parser.add_argument("--output", "-o", help="output file path", default="image.png")
    args = parser.parse_args()
    
    if str(args.cancer) not in cancer_type_list:
        print("cancer is not in known list")
        sys.exit()
    if isinstance(int(args.gen_start), int) == False:
        print("gen_start must be an integer")
        sys.exit()
    if isinstance(int(args.gen_end), int) == False or args.gen_end < args.gen_start:
        print("gen_end must be an integer and larger than gen_start")
        sys.exit()
    if str(args.mut_type).lower() not in ['sbs', 'dbs', 'insertion', 'deletion', 'mutation_burden'] :
        print("mut_type must be one of sbs, dbs, insertion, deletion, mutation_burden")
        sys.exit()
        
    mut_catalog(cancer_type=str(args.cancer), simulation_type=str(args.simulation), gen_start=int(args.gen_start), gen_end=int(args.gen_end), mut_type=str(args.mut_type), output_path=str(args.output))
    
if __name__ ==  '__main__':

    main()