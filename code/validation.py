# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 14:44:58 2018

@author: Rayin
"""

import os, sys
import pandas as pd
import numpy as np
import difflib
import re

from sklearn.cross_validation import StratifiedKFold
from scipy import interp
from sklearn.multiclass import OneVsRestClassifier
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
from itertools import cycle
from Bio import SeqIO
from PyBioMed.PyProtein import CTD
from sklearn.externals import joblib

sys.path.append('C:/Users/Rayin/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/code')

#import method_random
from method_optimized import split_train_test
from method_optimized import randomforest_cross_validation
from method_optimized import generate_features
from method_optimized import prediction_genome_csv
from method_optimized import split_description
from method_optimized import calculate_reassort_prob
from method_optimized import random_training_selection
from method_optimized import training_data_preparation
from method_optimized import training_data_preparation_2
from method_optimized import genome_training_model
#from method_random import split_protein
#from method_random import fill_with_missing_protein


#this function is to fill with the missing protein description for each segment    
def fill_with_missing_protein(input_file):
    count = 0
    protein_len = []
    for i in range(0, input_file.shape[0]):
        max_similarity = 0
        if input_file.ix[i]['seq'] == '':
            count = count + 1
            input_file.ix[i]['seq'] = input_file.ix[i]['protein']
            input_file.ix[i]['protein'] = ''
            protein_len.append(len(input_file.ix[i]['seq']))
            for protein, seq in genome_sample.items():
                similarity = difflib.SequenceMatcher(None, input_file.ix[i]['seq'], seq).ratio()
                if similarity >= max_similarity:
                    max_similarity = similarity
                    input_file.ix[i]['protein'] = protein
    print count
    return input_file
  


#find the missing genome between our database and literature report
def find_missing_description(result_input, description_input):
    reassort_description = result_input['description']
    for i in range(result_input.shape[0]):
        reassort_description[i] = reassort_description[i].replace(" 1", ",1")
        reassort_description[i] = reassort_description[i].replace(" 2", ",2")
        a, b = reassort_description[i].split(",")
        reassort_description[i] = a
    missing_description = pd.concat([pd.DataFrame(description_input), pd.DataFrame(reassort_description)]).drop_duplicates(keep=False) 
    return missing_description
    
#this is the function to split the protein information from IRD database
def split_protein(input_file):
    df = input_file
    df_protein = df["protein"]
    for row in range(0, df.shape[0]):
        #df_protein[row] = re.search(r'/:(.*?)\s/', df_protein[row])
        df_protein[row] = re.findall(r'(NP|PB1|PB2|HA|NA|PA|NS1|NS2|M1|M2)\s.*', df_protein[row])
        df_protein[row] = ''.join(df_protein[row])

    df["protein"] = df_protein
    return df    

###########################################################################################################################################

################################################################################################################################################
os.chdir("C:/Users/Rayin/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data/")

##this is to deal with all the Hossein's genome data 
Hossein_data_all = pd.read_csv("Hossein/Hossein_data_all.csv", na_filter = False)
Hossein_data_all = Hossein_data_all.ix[:, 1:5]
#Hossein_data_all = pd.concat([Hossein_data, Hossein_missing_data])
#Hossein_data_all.to_csv("Hossein/Hossein_data_all.csv")
HA, NA, PA, NP, PB1, PB2, NS1, NS2, M1, M2 = genome_training_model(Hossein_data_all, 0.95)    
#prediction_genome_csv(Hossein_data_all, HA, NA, PA, NP, PB1, PB2, NS1, NS2, M1, M2)
prediction, pred_avian, pred_human, pred_swine = prediction_genome_csv(Hossein_data_all, HA, NA, PA, NP, PB1, PB2, NS1, NS2, M1, M2) 
pred_df = pd.DataFrame(prediction)
pred_avian_df = pd.DataFrame(pred_avian)
pred_avian_df = pd.DataFrame(pred_avian)
pred_human_df = pd.DataFrame(pred_human)
pred_swine_df = pd.DataFrame(pred_swine)
reassortment_result = calculate_reassort_prob(pred_avian_df, pred_human_df, pred_swine_df)
reassortment_result = pd.DataFrame(reassortment_result)  
reassortment_result.to_csv("Hossein/Hossein_all_reassortment_prob_0.95.csv")
pred_df.to_csv("Hossein/Hossein_all_prediction_0.95.csv")










