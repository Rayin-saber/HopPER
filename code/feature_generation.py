# -*- coding: utf-8 -*-
"""
Created on Fri Dec 08 09:13:01 2017

@author: Rayin
"""

import os, sys
import pandas as pd
from PyBioMed import Pyprotein
from PyBioMed.PyProtein import CTD
from Bio import AlignIO
from Bio import SeqIO

os.chdir("C:/Users/Rayin/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data/Sequence")


def feature_generation(input_file):
    global feature_file, feature
#    feature_file = open(os.path.expanduser("C:/Users/Rayin/Google Drive/Tier_2_MOE2014/5_Journal/Plos_One/Data/feature/feature.csv"), 'w')
    input_file = pd.DataFrame(input_file)
    input_file = input_file['seq']
    count = 0
    for row in range(0, len(input_file)):
        seq = str(input_file.loc[row])
        feature = CTD.CalculateCTD(seq)
        if row == 0:
            write_header() 
        for key in sorted(feature.keys()):
            feature_value = feature[key]
            count = count + 1
            if count == 147:
                count = 0
                write_to_csv = str(feature_value) + '\n'
            else:
                write_to_csv = str(feature_value) + ','   
            #write_to_csv = str(feature_value) + ','     
            feature_file.write(write_to_csv)
            
    feature_file.close()                
    #return feature_value
        


def write_header():
#    feature_file = open(os.path.expanduser("C:/Users/Rayin/Google Drive/Tier_2_MOE2014/5_Journal/Plos_One/Data/feature/feature.csv"), 'w')
    column = 0
    for key in sorted(feature.keys()):
        feature_header = key
        column = column + 1
        if column == 147:
            column = 0
            write_header = str(feature_header) + '\n'
        else:
            write_header = str(feature_header) + ','
        feature_file.write(write_header)
    #return feature_header   
 
    
#this is feature generation for all length sequences of eight segments
if __name__ == '__main__':
    #feature generation of HA
    feature_file = open(os.path.expanduser("feature/all_length/HA_feature.csv"), 'w')
    HA_segment = pd.read_csv('csv/all_length/HA.csv', names=['accession', 'seq'])
    feature_generation(HA_segment)
    feature_file.close()

    #feature generation of M1
    feature_file = open(os.path.expanduser("feature/all_length/M1_feature.csv"), 'w')
    M1_segment = pd.read_csv('csv/all_length/M1.csv', names=['accession', 'seq'])
    feature_generation(M1_segment)
    feature_file.close()
    
    #feature generation of M2
    feature_file = open(os.path.expanduser("feature/all_length/M2_feature.csv"), 'w')
    M2_segment = pd.read_csv('csv/all_length/M2.csv', names=['accession', 'seq'])
    feature_generation(M2_segment)
    feature_file.close()
    
    #feature generation of NA
    feature_file = open(os.path.expanduser("feature/all_length/NA_feature.csv"), 'w')
    NA_segment = pd.read_csv('csv/all_length/NA.csv', names=['accession', 'seq'])
    feature_generation(NA_segment)
    feature_file.close()
    
    #feature generation of NP
    feature_file = open(os.path.expanduser("feature/all_length/NP_feature.csv"), 'w')
    NP_segment = pd.read_csv('csv/all_length/NP.csv', names=['accession', 'seq'])
    feature_generation(NP_segment)
    feature_file.close()
    
    #feature generation of NS1
    feature_file = open(os.path.expanduser("feature/all_length/NS1_feature.csv"), 'w')
    NS1_segment = pd.read_csv('csv/all_length/NS1.csv', names=['accession', 'seq'])
    feature_generation(NS1_segment)
    feature_file.close()

    #feature generation of NS2
    feature_file = open(os.path.expanduser("feature/all_length/NS2_feature.csv"), 'w')
    NS2_segment = pd.read_csv('csv/all_length/NS2.csv', names=['accession', 'seq'])
    feature_generation(NS2_segment)
    feature_file.close()

    #feature generation of PA
    feature_file = open(os.path.expanduser("feature/all_length/PA_feature.csv"), 'w')
    PA_segment = pd.read_csv('csv/all_length/PA.csv', names=['accession', 'seq'])
    feature_generation(PA_segment)
    feature_file.close()

    #feature generation of PB1
    feature_file = open(os.path.expanduser("feature/all_length/PB1_feature.csv"), 'w')
    PB1_segment = pd.read_csv('csv/all_length/PB1.csv', names=['accession', 'seq'])
    feature_generation(PB1_segment)
    feature_file.close() 
    
    #feature generation of PB1_F2
    feature_file = open(os.path.expanduser("feature/all_length/PB1_F2_feature.csv"), 'w')
    PB1_F2_segment = pd.read_csv('csv/all_length/PB1-F2.csv', names=['accession', 'seq'])
    feature_generation(PB1_F2_segment)
    feature_file.close() 
    
    #feature generation of PB2
    feature_file = open(os.path.expanduser("feature/all_length/PB2_feature.csv"), 'w')
    PB2_segment = pd.read_csv('csv/all_length/PB2.csv', names=['accession', 'seq'])
    feature_generation(PB2_segment)
    feature_file.close()  
    
    
    
    
######################################################################################################################    
#prepare the training data of full length for each segment, 0 is avian, 1 is human, 2 is swine   
os.chdir('C:/Users/Rayin/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data/Sequence')

HA_feature = pd.read_csv('feature/full_length/HA_feature.csv')
HA_year = pd.read_csv('feature/full_length/HA_year.csv', header=None, names=['year'])
#HA_data = HA_feature.join(HA_year, how='left')
HA_data = pd.concat([HA_feature, HA_year], axis=1)
HA_data['label'] = 0
HA_data['label'].loc[0:12248] = 0
HA_data['label'].loc[12248:25855] = 1
HA_data['label'].loc[25855:] = 2
#HA_data.to_csv('train/full_length/HA_data.csv')

M1_feature = pd.read_csv('feature/full_length/M1_feature.csv')
M1_year = pd.read_csv('feature/full_length/M1_year.csv', header=None, names=['year'])
M1_data = pd.concat([M1_feature, M1_year], axis=1)
M1_data['label'] = 0
M1_data['label'].loc[0:1682] = 0
M1_data['label'].loc[1682:2804] = 1
M1_data['label'].loc[2804:] = 2
#M1_data.to_csv('train/full_length/M1_data.csv')

M2_feature = pd.read_csv('feature/full_length/M2_feature.csv')
M2_year = pd.read_csv('feature/full_length/M2_year.csv', header=None, names=['year'])
M2_data = pd.concat([M2_feature, M2_year], axis=1)
M2_data['label'] = 0
M2_data['label'].loc[0:2237] = 0
M2_data['label'].loc[2237:3641] = 1
M2_data['label'].loc[3641:] = 2
#M2_data.to_csv('train/full_length/M2_data.csv')

NA_feature = pd.read_csv('feature/full_length/NA_feature.csv')
NA_year = pd.read_csv('feature/full_length/NA_year.csv', header=None, names=['year'])
NA_data = pd.concat([NA_feature, NA_year], axis=1)
NA_data['label'] = 0
NA_data['label'].loc[0:9452] = 0
NA_data['label'].loc[9452:19559] = 1
NA_data['label'].loc[19559:] = 2
#NA_data.to_csv('train/full_length/NA_data.csv')

NP_feature = pd.read_csv('feature/full_length/NP_feature.csv')
NP_year = pd.read_csv('feature/full_length/NP_year.csv', header=None, names=['year'])
NP_data = pd.concat([NP_feature, NP_year], axis=1)
NP_data['label'] = 0
NP_data['label'].loc[0:4841] = 0
NP_data['label'].loc[4841:7500] = 1
NP_data['label'].loc[7500:] = 2
#NP_data.to_csv('train/full_length/NP_data.csv')

NS1_feature = pd.read_csv('feature/full_length/NS1_feature.csv')
NS1_year = pd.read_csv('feature/full_length/NS1_year.csv', header=None, names=['year'])
NS1_data = pd.concat([NS1_feature, NS1_year], axis=1)
NS1_data['label'] = 0
NS1_data['label'].loc[0:6115] = 0
NS1_data['label'].loc[6115:10248] = 1
NS1_data['label'].loc[10248:] = 2
#NS1_data.to_csv('train/full_length/NS1_data.csv')

NS2_feature = pd.read_csv('feature/full_length/NS2_feature.csv')
NS2_year = pd.read_csv('feature/full_length/NS2_year.csv', header=None, names=['year'])
NS2_data = pd.concat([NS2_feature, NS2_year], axis=1)
NS2_data['label'] = 0
NS2_data['label'].loc[0:2336] = 0
NS2_data['label'].loc[2336:3463] = 1
NS2_data['label'].loc[3463:] = 2
#NS2_data.to_csv('train/full_length/NS2_data.csv')

PA_feature = pd.read_csv('feature/full_length/PA_feature.csv')
PA_year = pd.read_csv('feature/full_length/PA_year.csv', header=None, names=['year'])
PA_data = pd.concat([PA_feature, PA_year], axis=1)
PA_data['label'] = 0
PA_data['label'].loc[0:8428] = 0
PA_data['label'].loc[8428:13926] = 1
PA_data['label'].loc[13926:] = 2
#PA_data.to_csv('train/full_length/PA_data.csv')

PB1_feature = pd.read_csv('feature/full_length/PB1_feature.csv')
PB1_year = pd.read_csv('feature/full_length/PB1_year.csv', header=None, names=['year'])
PB1_data = pd.concat([PB1_feature, PB1_year], axis=1)
PB1_data['label'] = 0
PB1_data['label'].loc[0:7699] = 0
PB1_data['label'].loc[7699:12568] = 1
PB1_data['label'].loc[12568:] = 2
#PB1_data.to_csv('train/full_length/PB1_data.csv')

PB1_F2_feature = pd.read_csv('feature/full_length/PB1_F2_feature.csv')
PB1_F2_year = pd.read_csv('feature/full_length/PB1_F2_year.csv', header=None, names=['year'])
PB1_F2_data = pd.concat([PB1_F2_feature, PB1_F2_year], axis=1)
PB1_F2_data['label'] = 0
PB1_F2_data['label'].loc[0:4992] = 0
PB1_F2_data['label'].loc[4992:6694] = 1
PB1_F2_data['label'].loc[6694:] = 2
#PB1_F2_data.to_csv('train/full_length/PB1_F2_data.csv')

PB2_feature = pd.read_csv('feature/full_length/PB2_feature.csv')
PB2_year = pd.read_csv('feature/full_length/PB2_year.csv', header=None, names=['year'])
PB2_data = pd.concat([PB2_feature, PB2_year], axis=1)
PB2_data['label'] = 0
PB2_data['label'].loc[0:8106] = 0
PB2_data['label'].loc[8106:13596] = 1
PB2_data['label'].loc[13596:] = 2
#PB2_data.to_csv('train/full_length/PB2_data.csv')
#########################################################################################################################################


#########################################################################################################################################
#prepare the training data of all length for each segment, 0 is avian, 1 is human, 2 is swine   
os.chdir('C:/Users/Rayin/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data/Sequence')

HA_feature = pd.read_csv('feature/all_length/HA_feature.csv')
HA_year = pd.read_csv('feature/all_length/HA_year.csv', header=None, names=['year'])
#HA_data = HA_feature.join(HA_year, how='left')
HA_data = pd.concat([HA_feature, HA_year], axis=1)
HA_data['label'] = 0
HA_data['label'].loc[0:28713] = 0
HA_data['label'].loc[28713:89336] = 1
HA_data['label'].loc[89336:] = 2
#HA_data.to_csv('train/all_length/HA_data.csv')

M1_feature = pd.read_csv('feature/all_length/M1_feature.csv')
M1_year = pd.read_csv('feature/all_length/M1_year.csv', header=None, names=['year'])
M1_data = pd.concat([M1_feature, M1_year], axis=1)
M1_data['label'] = 0
M1_data['label'].loc[0:20549] = 0
M1_data['label'].loc[20549:52631] = 1
M1_data['label'].loc[52631:] = 2
#M1_data.to_csv('train/all_length/M1_data.csv')

M2_feature = pd.read_csv('feature/all_length/M2_feature.csv')
M2_year = pd.read_csv('feature/all_length/M2_year.csv', header=None, names=['year'])
M2_data = pd.concat([M2_feature, M2_year], axis=1)
M2_data['label'] = 0
M2_data['label'].loc[0:19623] = 0
M2_data['label'].loc[19623:51994] = 1
M2_data['label'].loc[51994:] = 2
#M2_data.to_csv('train/all_length/M2_data.csv')

NA_feature = pd.read_csv('feature/all_length/NA_feature.csv')
NA_year = pd.read_csv('feature/all_length/NA_year.csv', header=None, names=['year'])
NA_data = pd.concat([NA_feature, NA_year], axis=1)
NA_data['label'] = 0
NA_data['label'].loc[0:22336] = 0
NA_data['label'].loc[22336:63117] = 1
NA_data['label'].loc[63117:] = 2
#NA_data.to_csv('train/all_length/NA_data.csv')

NP_feature = pd.read_csv('feature/all_length/NP_feature.csv')
NP_year = pd.read_csv('feature/all_length/NP_year.csv', header=None, names=['year'])
NP_data = pd.concat([NP_feature, NP_year], axis=1)
NP_data['label'] = 0
NP_data['label'].loc[0:19375] = 0
NP_data['label'].loc[19375:44273] = 1
NP_data['label'].loc[44273:] = 2
#NP_data.to_csv('train/all_length/NP_data.csv')

NS1_feature = pd.read_csv('feature/all_length/NS1_feature.csv')
NS1_year = pd.read_csv('feature/all_length/NS1_year.csv', header=None, names=['year'])
NS1_data = pd.concat([NS1_feature, NS1_year], axis=1)
NS1_data['label'] = 0
NS1_data['label'].loc[0:19921] = 0
NS1_data['label'].loc[19921:45169] = 1
NS1_data['label'].loc[45169:] = 2
#NS1_data.to_csv('train/all_length/NS1_data.csv')

NS2_feature = pd.read_csv('feature/all_length/NS2_feature.csv')
NS2_year = pd.read_csv('feature/all_length/NS2_year.csv', header=None, names=['year'])
NS2_data = pd.concat([NS2_feature, NS2_year], axis=1)
NS2_data['label'] = 0
NS2_data['label'].loc[0:19437] = 0
NS2_data['label'].loc[19437:44585] = 1
NS2_data['label'].loc[44585:] = 2
#NS2_data.to_csv('train/all_length/NS2_data.csv')

PA_feature = pd.read_csv('feature/all_length/PA_feature.csv')
PA_year = pd.read_csv('feature/all_length/PA_year.csv', header=None, names=['year'])
PA_data = pd.concat([PA_feature, PA_year], axis=1)
PA_data['label'] = 0
PA_data['label'].loc[0:19670] = 0
PA_data['label'].loc[19670:43641] = 1
PA_data['label'].loc[43641:] = 2
#PA_data.to_csv('train/all_length/PA_data.csv')

PB1_feature = pd.read_csv('feature/all_length/PB1_feature.csv')
PB1_year = pd.read_csv('feature/all_length/PB1_year.csv', header=None, names=['year'])
PB1_data = pd.concat([PB1_feature, PB1_year], axis=1)
PB1_data['label'] = 0
PB1_data['label'].loc[0:19714] = 0
PB1_data['label'].loc[19714:43632] = 1
PB1_data['label'].loc[43632:] = 2
#PB1_data.to_csv('train/all_length/PB1_data.csv')

PB1_F2_feature = pd.read_csv('feature/all_length/PB1_F2_feature.csv')
PB1_F2_year = pd.read_csv('feature/all_length/PB1_F2_year.csv', header=None, names=['year'])
PB1_F2_data = pd.concat([PB1_F2_feature, PB1_F2_year], axis=1)
PB1_F2_data['label'] = 0
PB1_F2_data['label'].loc[0:14188] = 0
PB1_F2_data['label'].loc[14188:28930] = 1
PB1_F2_data['label'].loc[28930:] = 2
#PB1_F2_data.to_csv('train/all_length/PB1_F2_data.csv')

PB2_feature = pd.read_csv('feature/all_length/PB2_feature.csv')
PB2_year = pd.read_csv('feature/all_length/PB2_year.csv', header=None, names=['year'])
PB2_data = pd.concat([PB2_feature, PB2_year], axis=1)
PB2_data['label'] = 0
PB2_data['label'].loc[0:19730] = 0
PB2_data['label'].loc[19730:43861] = 1
PB2_data['label'].loc[43861:] = 2
#PB2_data.to_csv('train/all_length/PB2_data.csv')
########################################################################################################################################



 

        
        
        
        
  

    
    
    
    
