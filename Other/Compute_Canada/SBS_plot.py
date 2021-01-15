import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as patches
import numpy as np
import os
import sys
import argparse

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

def mut_catalog(signature, gen_start, gen_end, username):
    try:
       
        graph_data = pd.DataFrame(index=(range(96)), columns=range(gen_start, gen_end+1))
        
        for gen in range(gen_start, gen_end+1):
            try:
                df_read_in = pd.read_csv("/home/" + str(username) + "/projects/def-khill22/" + str(username) + "/Metadata/" + signature + '_Sequence_' + str(gen) + '_sbs_freq_table.csv', 'File')
                graph_data[gen] = df_read_in['Frequency']
            except:
                print(signature + '_Sequence_' + str(gen) + '_sbs_freq_table.csv does not exist!')
                
        graph_data['SubType'] = df_read_in['SubType']
        graph_data['Type'] = df_read_in['Type']
        
        for gen in range(gen_start, gen_end+1):
             graph_data[gen] = graph_data[gen].div(graph_data[gen].sum())
        
        graph_data.sort_values(by=["Type","SubType"],inplace=True )
        graph_data.fillna(0, inplace=True)
        

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
        
        try:
            plot_directory = abs_path("Plot", "Directory")
        except:
            try:
                os.mkdir("Plot")
            except OSError:
                print ("Creation of the directory failed.")
                
                
        fig.savefig("/home/" + str(username) + "/projects/def-khill22/" + str(username) + "/Plot/" + str(signature) + ".png", bbox_inches='tight')

    except:
        print('Error')
        
    
        
        
def main(): 
    
    # Initiate the parser
    parser = argparse.ArgumentParser(description="SBS Plot Parameters")
    
    # Add long and short argument
    parser.add_argument("--signature", "-s", help="signature type")
    parser.add_argument("--gen_start", "-a", help="start of range of sequence metadata")
    parser.add_argument("--gen_end", "-b", help="end of range of sequence metadata")
    parser.add_argument("--username", "-u", help="Your Compute Canada Username")
    args = parser.parse_args()
    
    if str(args.signature) not in signatures:
        print("signature is not in known list")
        sys.exit()
    if isinstance(int(args.gen_start), int) == False:
        print("gen_start must be an integer")
        sys.exit()
    if isinstance(int(args.gen_end), int) == False or args.gen_end < args.gen_start:
        print("gen_end must be an integer and larger than gen_start")
        sys.exit()
        
    mut_catalog(signature = args.signature, gen_start = args.gen_start, gen_end = args.gen_end, username = args.username)
    
if __name__ ==  '__main__':

    main()
