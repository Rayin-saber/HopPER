# -*- coding: utf-8 -*-
"""
Created on Sat Jul 14 20:07:21 2018

@author: Rayin
"""
from __future__ import division

import os, sys
import pandas as pd
import numpy as np
#import cPickle
import matplotlib.pyplot as plt
import warnings
import random
import gc
import re
import math
import difflib
    

from sklearn import datasets
from sklearn import preprocessing
from sklearn import neighbors
from sklearn import svm
from sklearn import model_selection
from sklearn import ensemble
from sklearn import metrics
from sklearn import linear_model
from sklearn.neural_network import MLPClassifier 
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest, chi2
from matplotlib import pyplot 
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import train_test_split
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import roc_curve, auc
from imblearn.metrics import geometric_mean_score

from scipy import interp
from sklearn.multiclass import OneVsRestClassifier
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
from itertools import cycle
#from Bio import SeqIO
from PyBioMed.PyProtein import CTD
from sklearn.externals import joblib

import warnings
warnings.filterwarnings("ignore")


os.chdir("C:/Users/yinr0002/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data")
#os.chdir("C:/Users/Rayin/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data")
#os.chdir("/Users/rayin/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data")
#########################################################################################################################################
#read in the sequence data of full length for each segment after feature transformation
#global HA_raw, HA_data, HA_feature, HA_label
HA_raw = pd.read_csv('sequence/csv_update/HA_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
HA_data = pd.read_csv('sequence/train/full_length/HA_data.csv')
HA_feature = HA_data.iloc[:, 1:148]
HA_label = HA_data.iloc[:, 149]
HA_num = HA_raw.shape[0]

#global NA_raw, NA_data, NA_feature, NA_label
NA_raw = pd.read_csv('sequence/csv_update/NA_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
NA_data = pd.read_csv('sequence/train/full_length/NA_data.csv')
NA_feature = NA_data.iloc[:, 1:148]
NA_label = NA_data.iloc[:, 149]
NA_num = NA_raw.shape[0]

#global NP_raw, NP_data, NP_feature, NP_label
NP_raw = pd.read_csv('sequence/csv_update/NP_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
NP_data = pd.read_csv('sequence/train/full_length/NP_data.csv')
NP_feature = NP_data.iloc[:, 1:148]
NP_label = NP_data.iloc[:, 149]
NP_num = NP_raw.shape[0]

#global PA_raw, PA_data, PA_feature, PA_label
PA_raw = pd.read_csv('sequence/csv_update/PA_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
PA_data = pd.read_csv('sequence/train/full_length/PA_data.csv')
PA_feature = PA_data.iloc[:, 1:148]
PA_label = PA_data.iloc[:, 149]
PA_num = PA_raw.shape[0]

#global PB1_raw, PB1_data, PB1_feature, PB1_label
PB1_raw = pd.read_csv('sequence/csv_update/PB1_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
PB1_data = pd.read_csv('sequence/train/full_length/PB1_data.csv')
PB1_feature = PB1_data.iloc[:, 1:148]
PB1_label = PB1_data.iloc[:, 149]
PB1_num = PB1_raw.shape[0]

#global PB2_raw, PB2_data, PB2_feature, PB2_label
PB2_raw = pd.read_csv('sequence/csv_update/PB2_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
PB2_data = pd.read_csv('sequence/train/full_length/PB2_data.csv')
PB2_feature = PB2_data.iloc[:, 1:148]
PB2_label = PB2_data.iloc[:, 149]
PB2_num = PB2_raw.shape[0]

#global M1_raw, M1_data, M1_feature, M1_label
#M1_raw = pd.read_csv('sequence/csv_update/M1_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
#M1_data = pd.read_csv('sequence/train/full_length/M1_data.csv')
#M1_feature = M1_data.iloc[:, 1:148]
#M1_label = M1_data.iloc[:, 149]
#M1_num = M1_raw.shape[0]

#global M2_raw, M2_data, M2_feature, M2_label
M2_raw = pd.read_csv('sequence/csv_update/M2_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
M2_data = pd.read_csv('sequence/train/full_length/M2_data.csv')
M2_feature = M2_data.iloc[:, 1:148]
M2_label = M2_data.iloc[:, 149]
M2_num = M2_raw.shape[0]

#global NS1_raw, NS1_data, NS1_feature, NS1_label
NS1_raw = pd.read_csv('sequence/csv_update/NS1_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
NS1_data = pd.read_csv('sequence/train/full_length/NS1_data.csv')
NS1_feature = NS1_data.iloc[:, 1:148]
NS1_label = NS1_data.iloc[:, 149]
NS1_num = NS1_raw.shape[0]

#global NS2_raw, NS2_data, NS2_feature, NS2_label
#NS2_raw = pd.read_csv('sequence/csv_update/NS2_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
#NS2_data = pd.read_csv('sequence/train/full_length/NS2_data.csv')
#NS2_feature = NS2_data.iloc[:, 1:148]
#NS2_label = NS2_data.iloc[:, 149]
#NS2_num = NS2_raw.shape[0]
######################################################################################################################################

def split_train_test(feature, label, train_ratio):
    np.random.seed(42)
    shuffled_indices = np.random.permutation(len(feature))
    train_set_size = int(len(feature)*train_ratio)
    train_indices = shuffled_indices[: train_set_size]
    train_feature = feature.iloc[train_indices]
    train_label = label.iloc[train_indices]
    return train_feature, train_label

    
    
def randomforest_cross_validation(train_x, train_y):
    np.random.seed(100)
    clf = ensemble.RandomForestClassifier()
    
    clf.fit(train_x, train_y)
    
    #calculate the accuracy
    accuracy = cross_val_score(clf, train_x, train_y, cv=10, scoring='accuracy') 
    #print "accuracy: %f" %accuracy.mean() + '\n'
    #print accuracy
    
    #calculate the precision
    precision = cross_val_score(clf, train_x, train_y, cv=10, scoring='precision_macro')
    #print "precision: %f" %precision.mean() + '\n'
    
    #calculate the recall score
    recall = cross_val_score(clf, train_x, train_y, cv=10, scoring='recall_macro')
    #print "recall: %f" %recall.mean() + '\n'
    
    #calculate the f_measure
    f_measure = cross_val_score(clf, train_x, train_y, cv=10, scoring='f1_macro')
    #print "f_measure: %f " %f_measure.mean() + '\n' 

    #generate classification report and MCC and G-mean value
    y_pred = cross_val_predict(clf, train_x, train_y, cv=10)
    G_mean = geometric_mean_score(train_y, y_pred)
    MCC = matthews_corrcoef(train_y, y_pred)
    #print "G_mean: %f" %G_mean.mean() + '\n'
    #print "MCC: %f" %np.mean(MCC) + '\n'
    
    #print "Classification_report:"
    #print(metrics.classification_report(train_y, y_pred))    
    
    return clf
 

def generate_features(sequence):
    seq = str(sequence)
    feature = CTD.CalculateCTD(seq)
    feature_all = []
    count = 0
    for key in sorted(feature.keys()):
        feature_all.append(feature[key])
        count = count + 1
        if count == 147:
            feature_numpy = np.array(feature_all)
            # print(feature_numpy)
            return feature_numpy

            
            
#this is the genome prediction function with respect with csv file as input 
def prediction_genome_csv(csv_file, HA, NA, PA, NP, PB1, PB2, NS1, M2):
    prediction = {'description':[], 'HA':[], 'NA':[], 'NP':[], 'PA':[], 'PB1':[],'PB2':[], 'M2':[], 'NS1':[]}
    avian_prob = {'description': [], 'HA_prob':[], 'NA_prob':[], 'NP_prob':[], 'PA_prob':[], 'PB1_prob':[],'PB2_prob':[], 'M2_prob':[], 'NS1_prob':[]}
    human_prob = {'description': [], 'HA_prob':[], 'NA_prob':[], 'NP_prob':[], 'PA_prob':[], 'PB1_prob':[],'PB2_prob':[], 'M2_prob':[], 'NS1_prob':[]}
    swine_prob = {'description': [], 'HA_prob':[], 'NA_prob':[], 'NP_prob':[], 'PA_prob':[], 'PB1_prob':[],'PB2_prob':[], 'M2_prob':[], 'NS1_prob':[]}
    #print(csv_file)
    #df = pd.read_csv(csv_file)
    df = csv_file
    df = df.fillna(0)
    
    description = []
    #probability = []
    for x in df['description']:
        description.append(x)
        
    for x in set(description):  
        print(x)
        df_new = df.ix[df['description'] == x]
        prediction['description'].append(x)
        avian_prob['description'].append(x)
        human_prob['description'].append(x)
        swine_prob['description'].append(x)
        HA_flag = 0
        NA_flag = 0
        NP_flag = 0
        PA_flag = 0
        PB1_flag = 0
        PB2_flag = 0
        #M1_flag = 0
        M2_flag = 0
        NS1_flag = 0
        #NS2_flag = 0
        
        for index, row in df_new.iterrows():
            if "HA" == row['protein']:
                if HA_flag == 0:
                    HA_flag = 1
                    HA_features = generate_features(row['seq'])
                    HA_features = np.reshape(HA_features, (1, -1))
                    HA_prediction = HA.predict(HA_features)
                
                    prediction['HA'].append(int(HA_prediction))
                    #HA_prob = HA.predict_proba(HA_features)*100
                    HA_prob = HA.predict_proba(HA_features)
                    avian_prob['HA_prob'].append(HA_prob[0][0])
                    human_prob['HA_prob'].append(HA_prob[0][1])
                    swine_prob['HA_prob'].append(HA_prob[0][2])
                    print("HA\t\t", HA_prediction, HA_prob)
            elif "NA" == row['protein']:
                if NA_flag == 0:
                    NA_flag = 1
                    NA_features = generate_features(row['seq'])
                    NA_features = np.reshape(NA_features, (1, -1))
                    NA_prediction = NA.predict(NA_features)
                
                    prediction['NA'].append(int(NA_prediction))
                    #NA_prob = NA.predict_proba(NA_features)*100
                    NA_prob = NA.predict_proba(NA_features)
                    avian_prob['NA_prob'].append(NA_prob[0][0])
                    human_prob['NA_prob'].append(NA_prob[0][1])
                    swine_prob['NA_prob'].append(NA_prob[0][2])
                    print("NA\t\t", NA_prediction, NA_prob)
            
            elif "NP" == row['protein']:
                if NP_flag == 0:
                    NP_flag = 1
                    NP_features = generate_features(row['seq'])
                    NP_features = np.reshape(NP_features, (1, -1))
                    NP_prediction = NP.predict(NP_features)
                
                    prediction['NP'].append(int(NP_prediction))
                    #NP_prob = NP.predict_proba(NP_features)*100
                    NP_prob = NP.predict_proba(NP_features)
                    avian_prob['NP_prob'].append(NP_prob[0][0])
                    human_prob['NP_prob'].append(NP_prob[0][1])
                    swine_prob['NP_prob'].append(NP_prob[0][2])
                    print("NP\t\t", NP_prediction, NP_prob)
            
            elif "PA" == row['protein']:
                if PA_flag == 0:
                    PA_flag = 1               
                    PA_features = generate_features(row['seq'])
                    PA_features = np.reshape(PA_features, (1, -1))
                    PA_prediction = PA.predict(PA_features)
                
                    prediction['PA'].append(int(PA_prediction))
                    #PA_prob = PA.predict_proba(PA_features)*100
                    PA_prob = PA.predict_proba(PA_features)
                    avian_prob['PA_prob'].append(PA_prob[0][0])
                    human_prob['PA_prob'].append(PA_prob[0][1])
                    swine_prob['PA_prob'].append(PA_prob[0][2])
                    print("PA\t\t", PA_prediction, PA_prob)
            
            elif "PB1" == row['protein']:
                if PB1_flag == 0:
                    PB1_flag = 1
                    PB1_features = generate_features(row['seq'])
                    PB1_features = np.reshape(PB1_features, (1, -1))
                    PB1_prediction = PB1.predict(PB1_features)
                
                    prediction['PB1'].append(int(PB1_prediction))
                    #PB1_prob = PB1.predict_proba(PB1_features)*100
                    PB1_prob = PB1.predict_proba(PB1_features)
                    avian_prob['PB1_prob'].append(PB1_prob[0][0])
                    human_prob['PB1_prob'].append(PB1_prob[0][1])
                    swine_prob['PB1_prob'].append(PB1_prob[0][2])
                    print("PB1\t\t", PB1_prediction, PB1_prob)
            
            elif "PB2" == row['protein']:
                if PB2_flag == 0:
                    PB2_flag = 1
                    PB2_features = generate_features(row['seq'])
                    #PB2_features = PB2_features.reshape(1, -1)
                    PB2_features = np.reshape(PB2_features, (1, -1))
                    PB2_prediction = PB2.predict(PB2_features)
                
                    prediction['PB2'].append(int(PB2_prediction))
                    #PB2_prob = PB2.predict_proba(PB2_features)*100
                    PB2_prob = PB2.predict_proba(PB2_features)
                    avian_prob['PB2_prob'].append(PB2_prob[0][0])
                    human_prob['PB2_prob'].append(PB2_prob[0][1])
                    swine_prob['PB2_prob'].append(PB2_prob[0][2])
                    print("PB2\t\t", PB2_prediction, PB2_prob)
            
#            elif "M1" == row['protein']:
#                if M1_flag == 0:
#                    M1_flag = 1
#                    M1_features = generate_features(row['seq'])
#                    M1_features = np.reshape(M1_features, (1, -1))
#                    M1_prediction = M1.predict(M1_features)
#                
#                    prediction['M1'].append(int(M1_prediction))
#                    #M1_prob = M1.predict_proba(M1_features)*100
#                    M1_prob = M1.predict_proba(M1_features)
#                    avian_prob['M1_prob'].append(M1_prob[0][0])
#                    human_prob['M1_prob'].append(M1_prob[0][1])
#                    swine_prob['M1_prob'].append(M1_prob[0][2])
#                    print "M1\t\t", M1_prediction, M1_prob 
                    
            elif "M2" == row['protein']:        
                if M2_flag == 0:
                    M2_flag = 1
                    M2_features = generate_features(row['seq'])
                    M2_features = np.reshape(M2_features, (1, -1))
                    M2_prediction = M2.predict(M2_features)
                
                    prediction['M2'].append(int(M2_prediction))
                    #M2_prob = M2.predict_proba(M2_features)*100
                    M2_prob = M2.predict_proba(M2_features)
                    avian_prob['M2_prob'].append(M2_prob[0][0])
                    human_prob['M2_prob'].append(M2_prob[0][1])
                    swine_prob['M2_prob'].append(M2_prob[0][2])
                    print("M2\t\t", M2_prediction, M2_prob)
             
            elif "NS1" == row['protein']:
                if NS1_flag == 0:
                    NS1_flag = 1
                    NS1_features = generate_features(row['seq'])
                    NS1_features = np.reshape(NS1_features, (1, -1))
                    NS1_prediction = NS1.predict(NS1_features)
                
                    prediction['NS1'].append(int(NS1_prediction))
                    #NS1_prob = NS1.predict_proba(NS1_features)*100
                    NS1_prob = NS1.predict_proba(NS1_features)
                    avian_prob['NS1_prob'].append(NS1_prob[0][0])
                    human_prob['NS1_prob'].append(NS1_prob[0][1])
                    swine_prob['NS1_prob'].append(NS1_prob[0][2])
                    print("NS1\t\t", NS1_prediction, NS1_prob)
                    
#            elif "NS2" == row['protein']:        
#                if NS2_flag == 0:
#                    NS2_flag = 1
#                    NS2_features = generate_features(row['seq'])
#                    NS2_features = np.reshape(NS2_features, (1, -1))
#                    NS2_prediction = NS2.predict(NS2_features)
#                
#                    prediction['NS2'].append(int(NS2_prediction))
#                    #NS2_prob = NS2.predict_proba(NS2_features)*100
#                    NS2_prob = NS2.predict_proba(NS2_features)
#                    avian_prob['NS2_prob'].append(NS2_prob[0][0])
#                    human_prob['NS2_prob'].append(NS2_prob[0][1])
#                    swine_prob['NS2_prob'].append(NS2_prob[0][2])
#                    print "NS2\t\t", NS2_prediction, NS2_prob 
        
        if HA_flag == 0:
            prediction['HA'].append("missing")
            avian_prob['HA_prob'].append("missing")
            human_prob['HA_prob'].append("missing")
            swine_prob['HA_prob'].append("missing")
        if NA_flag == 0:
            prediction['NA'].append("missing")
            avian_prob['NA_prob'].append("missing")
            human_prob['NA_prob'].append("missing")
            swine_prob['NA_prob'].append("missing")
        if NP_flag == 0:
            prediction['NP'].append("missing")
            avian_prob['NP_prob'].append("missing")
            human_prob['NP_prob'].append("missing")
            swine_prob['NP_prob'].append("missing")
        if PA_flag == 0:
            prediction['PA'].append("missing")
            avian_prob['PA_prob'].append("missing")
            human_prob['PA_prob'].append("missing")
            swine_prob['PA_prob'].append("missing")
        if PB1_flag == 0:
            prediction['PB1'].append("missing")
            avian_prob['PB1_prob'].append("missing")
            human_prob['PB1_prob'].append("missing")
            swine_prob['PB1_prob'].append("missing")
        if PB2_flag == 0:
            prediction['PB2'].append("missing")
            avian_prob['PB2_prob'].append("missing")
            human_prob['PB2_prob'].append("missing")
            swine_prob['PB2_prob'].append("missing")
#        if M1_flag == 0:
#            prediction['M1'].append("missing")
#            avian_prob['M1_prob'].append("missing")
#            human_prob['M1_prob'].append("missing")
#            swine_prob['M1_prob'].append("missing")
        if M2_flag == 0:
            prediction['M2'].append("missing")
            avian_prob['M2_prob'].append("missing")
            human_prob['M2_prob'].append("missing")
            swine_prob['M2_prob'].append("missing")
        if NS1_flag == 0:
            prediction['NS1'].append("missing")
            avian_prob['NS1_prob'].append("missing")
            human_prob['NS1_prob'].append("missing")
            swine_prob['NS1_prob'].append("missing")
#        if NS2_flag == 0:
#            prediction['NS2'].append("missing")
#            avian_prob['NS2_prob'].append("missing")
#            human_prob['NS2_prob'].append("missing")
#            swine_prob['NS2_prob'].append("missing")
            #if HA_flag == 1 and NA_flag == 1 and NP_flag == 1 and PA_flag == 1 and PB1_flag == 1 and PB2_flag == 1
    return prediction, avian_prob, human_prob, swine_prob




    
#this is the function that split the description
def split_description(input_file):
    df = input_file
    df_description = df["description"]
    for row in range(0, df.shape[0]):
        #df_description[row] = re.search(r'\bA/\S.*', df_description[row]).group()
        df_description[row] = re.search(r'\bA/\S.*', df_description[row]).group()
    df["description"] = df_description
    return df

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
    
#calculate the probability of reassortment of genome set
def calculate_reassort_prob(avian, human, swine):
    avian = avian
    human = human
    swine = swine
    row = avian.shape[0]
    #old version
    #columns = avian.columns[:-1]
    columns = avian.columns[-8:]
    #description = avian.columns[-1]
    reassortment_result = {'description':[], 'reassort_prob':[]}
    avian_prob = []
    human_prob = []
    swine_prob = []
    reassortment_prob = []
    for i in range(0, row):
        avian_initial = 1.0
        human_initial = 1.0
        swine_initial = 1.0
        reassort_prob = 0
        for column in columns:
            if avian[column].iloc[i] == 'missing':
                avian[column].iloc[i] = 1.0
            avian_initial = avian_initial * avian[column].iloc[i]
            if human[column].iloc[i] == 'missing':
                human[column].iloc[i] = 1.0
            human_initial = human_initial * human[column].iloc[i]
            if swine[column].iloc[i] == 'missing':
                swine[column].iloc[i] = 1.0
            swine_initial = swine_initial * swine[column].iloc[i]
        reassort_prob = 1 - avian_initial - human_initial - swine_initial
        avian_prob.append(avian[column].iloc[i])
        human_prob.append(human[column].iloc[i])
        swine_prob.append(swine[column].iloc[i])
        reassortment_prob.append(reassort_prob)    
        reassortment_result['description'].append(avian.iloc[i].description)
        reassortment_result['reassort_prob'].append(reassort_prob) 

    return reassortment_result    



def training_data_preparation(input_data, threshold):
    #threshold = 0.99
    HA_count = 0
    NA_count = 0
    PA_count = 0
    NP_count = 0
    PB1_count = 0
    PB2_count = 0
    NS1_count = 0
    #NS2_count = 0
    #M1_count = 0
    M2_count = 0
    HA_training = {'accession':[], 'description':[], 'protein':[], 'seq':[]}
    NA_training = {'accession':[], 'description':[], 'protein':[], 'seq':[]}
    PA_training = {'accession':[], 'description':[], 'protein':[], 'seq':[]}
    NP_training = {'accession':[], 'description':[], 'protein':[], 'seq':[]}
    PB1_training = {'accession':[], 'description':[], 'protein':[], 'seq':[]}
    PB2_training = {'accession':[], 'description':[], 'protein':[], 'seq':[]}
    NS1_training = {'accession':[], 'description':[], 'protein':[], 'seq':[]}
    #NS2_training = {'accession':[], 'description':[], 'protein':[], 'seq':[]}
    #M1_training = {'accession':[], 'description':[], 'protein':[], 'seq':[]}
    M2_training = {'accession':[], 'description':[], 'protein':[], 'seq':[]}
    if input_data.shape[1] == 4:
        validate_data = input_data
    else:
        validate_data = input_data.iloc[:,[1,8,9,14]]
        validate_data.columns = ['accession', 'description', 'protein', 'seq']
        validate_data = split_protein(validate_data)
    for i in range(0, validate_data.shape[0]):
        if validate_data.iloc[i]['protein'] == 'HA':
            for j in range(0, HA_raw.shape[0]):
                #if difflib.SequenceMatcher(None, HA_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                if HA_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    HA_count = HA_count + 1
                    HA_training['accession'].append(HA_raw.iloc[j]['accession'])
                    HA_training['description'].append(HA_raw.iloc[j]['description'])
                    HA_training['protein'].append(HA_raw.iloc[j]['protein'])
                    HA_training['seq'].append(HA_raw.iloc[j]['seq'])
                    #HA_raw = HA_raw.drop(HA_raw.index[j])
        elif validate_data.iloc[i]['protein'] == 'NA':
            for j in range(0, NA_raw.shape[0]):
                #if difflib.SequenceMatcher(None, NA_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                if NA_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    NA_count = NA_count + 1
                    NA_training['accession'].append(NA_raw.iloc[j]['accession'])
                    NA_training['description'].append(NA_raw.iloc[j]['description'])
                    NA_training['protein'].append(NA_raw.iloc[j]['protein'])
                    NA_training['seq'].append(NA_raw.iloc[j]['seq'])
                    #NA_raw = NA_raw.drop(NA_raw.index[j])

        elif validate_data.iloc[i]['protein'] == 'PA':
            for j in range(0, PA_raw.shape[0]):
                #if difflib.SequenceMatcher(None, PA_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                if PA_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    PA_count = PA_count + 1
                    PA_training['accession'].append(PA_raw.iloc[j]['accession'])
                    PA_training['description'].append(PA_raw.iloc[j]['description'])
                    PA_training['protein'].append(PA_raw.iloc[j]['protein'])
                    PA_training['seq'].append(PA_raw.iloc[j]['seq'])
                    #PA_raw = PA_raw.drop(PA_raw.index[j])   
        elif validate_data.iloc[i]['protein'] == 'NP':
            for j in range(0, NP_raw.shape[0]):
                #if difflib.SequenceMatcher(None, NP_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                if NP_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    NP_count = NP_count + 1
                    NP_training['accession'].append(NP_raw.iloc[j]['accession'])
                    NP_training['description'].append(NP_raw.iloc[j]['description'])
                    NP_training['protein'].append(NP_raw.iloc[j]['protein'])
                    NP_training['seq'].append(NP_raw.iloc[j]['seq'])
                    #NP_raw = NP_raw.drop(NP_raw.index[j])
        elif validate_data.iloc[i]['protein'] == 'PB1':
            for j in range(0, PB1_raw.shape[0]):
                #if difflib.SequenceMatcher(None, PB1_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                if PB1_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    PB1_count = PB1_count + 1
                    PB1_training['accession'].append(PB1_raw.iloc[j]['accession'])
                    PB1_training['description'].append(PB1_raw.iloc[j]['description'])
                    PB1_training['protein'].append(PB1_raw.iloc[j]['protein'])
                    PB1_training['seq'].append(PB1_raw.iloc[j]['seq'])
                    #PB1_raw = PB1_raw.drop(PB1_raw.index[j])
        elif validate_data.iloc[i]['protein'] == 'PB2':
            for j in range(0, PB2_raw.shape[0]):
                #if difflib.SequenceMatcher(None, PB2_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                if PB2_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    PB2_count = PB2_count + 1
                    PB2_training['accession'].append(PB2_raw.iloc[j]['accession'])
                    PB2_training['description'].append(PB2_raw.iloc[j]['description'])
                    PB2_training['protein'].append(PB2_raw.iloc[j]['protein'])
                    PB2_training['seq'].append(PB2_raw.iloc[j]['seq'])
                    #PB2_raw = PB2_raw.drop(PB2_raw.index[j])
        elif validate_data.iloc[i]['protein'] == 'NS1':
            for j in range(0, NS1_raw.shape[0]):
                #if difflib.SequenceMatcher(None, NS1_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                if NS1_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    NS1_count = NS1_count + 1
                    NS1_training['accession'].append(NS1_raw.iloc[j]['accession'])
                    NS1_training['description'].append(NS1_raw.iloc[j]['description'])
                    NS1_training['protein'].append(NS1_raw.iloc[j]['protein'])
                    NS1_training['seq'].append(NS1_raw.iloc[j]['seq'])
                    #NS1_raw = NS1_raw.drop(NS1_raw.index[j])
#        elif validate_data.iloc[i]['protein'] == 'NS2':
#            for j in range(0, NS2_raw.shape[0]):
#                if difflib.SequenceMatcher(None, NS2_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
#                #if NS2_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
#                    NS2_count = NS2_count + 1
#                    NS2_training['accession'].append(NS2_raw.iloc[j]['accession'])
#                    NS2_training['description'].append(NS2_raw.iloc[j]['description'])
#                    NS2_training['protein'].append(NS2_raw.iloc[j]['protein'])
#                    NS2_training['seq'].append(NS2_raw.iloc[j]['seq'])
#                    #NS2_raw = NS2_raw.drop(NS2_raw.index[j])
#        elif validate_data.iloc[i]['protein'] == 'M1':
#            for j in range(0, M1_raw.shape[0]):
#                if difflib.SequenceMatcher(None, M1_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
#                #if M1_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
#                    M1_count = M1_count + 1
#                    M1_training['accession'].append(M1_raw.iloc[j]['accession'])
#                    M1_training['description'].append(M1_raw.iloc[j]['description'])
#                    M1_training['protein'].append(M1_raw.iloc[j]['protein'])
#                    M1_training['seq'].append(M1_raw.iloc[j]['seq'])
#                    #M1_raw = M1_raw.drop(M1_raw.index[j])
        elif validate_data.iloc[i]['protein'] == 'M2':
            for j in range(0, M2_raw.shape[0]):
                #if difflib.SequenceMatcher(None, M2_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                if M2_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    M2_count = M2_count + 1
                    M2_training['accession'].append(M2_raw.iloc[j]['accession'])
                    M2_training['description'].append(M2_raw.iloc[j]['description'])
                    M2_training['protein'].append(M2_raw.iloc[j]['protein'])
                    M2_training['seq'].append(M2_raw.iloc[j]['seq'])
                    #M2_raw = M2_raw.drop(M2_raw.index[j])
                    
                    
                    
    print(HA_count, NA_count, PA_count, NP_count, PB1_count, PB2_count, NS1_count, M2_count)
    HA_training = pd.DataFrame(HA_training).drop_duplicates('accession', inplace=False)
    NA_training = pd.DataFrame(NA_training).drop_duplicates('accession', inplace=False)
    PA_training = pd.DataFrame(PA_training).drop_duplicates('accession', inplace=False)
    NP_training = pd.DataFrame(NP_training).drop_duplicates('accession', inplace=False)
    PB1_training = pd.DataFrame(PB1_training).drop_duplicates('accession', inplace=False)
    PB2_training = pd.DataFrame(PB2_training).drop_duplicates('accession', inplace=False)
    NS1_training = pd.DataFrame(NS1_training).drop_duplicates('accession', inplace=False)
    #NS2_training = pd.DataFrame(NS2_training).drop_duplicates('accession', inplace=False)
    #M1_training = pd.DataFrame(M1_training).drop_duplicates('accession', inplace=False)
    M2_training = pd.DataFrame(M2_training).drop_duplicates('accession', inplace=False)
    print(HA_training.shape[0], NA_training.shape[0], PA_training.shape[0], NP_training.shape[0], PB1_training.shape[0], PB2_training.shape[0], NS1_training.shape[0], M2_training.shape[0])
    return HA_training, NA_training, PA_training, NP_training, PB1_training, PB2_training, NS1_training, M2_training
    


def training_data_preparation_version_2(input_data, threshold):
    #threshold = 0.99
    HA_count = 0
    NA_count = 0
    PA_count = 0
    NP_count = 0
    PB1_count = 0
    PB2_count = 0
    NS1_count = 0
    #NS2_count = 0
    #M1_count = 0
    M2_count = 0

    HA_training_index = []
    NA_training_index = []
    PA_training_index = []
    NP_training_index = []
    PB1_training_index = []
    PB2_training_index = []
    NS1_training_index = []
    #NS2_training_index = []
    #M1_training_index = []
    M2_training_index = []
    if input_data.shape[1] == 3 or 4:
        validate_data = input_data
    else:
        validate_data = input_data.iloc[:,[1,8,9,14]]
        validate_data.columns = ['accession', 'description', 'protein', 'seq']
        validate_data = split_protein(validate_data)
    for i in range(0, validate_data.shape[0]):
        if validate_data.iloc[i]['protein'] == 'HA':
            for j in range(0, HA_raw.shape[0]):
                if difflib.SequenceMatcher(None, HA_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                #if HA_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    HA_count = HA_count + 1
                    HA_training_index.append(j)
        elif validate_data.iloc[i]['protein'] == 'NA':
            for j in range(0, NA_raw.shape[0]):
                if difflib.SequenceMatcher(None, NA_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                #if NA_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    NA_count = NA_count + 1
                    NA_training_index.append(j)
        elif validate_data.iloc[i]['protein'] == 'PA':
            for j in range(0, PA_raw.shape[0]):
                if difflib.SequenceMatcher(None, PA_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                #if PA_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    PA_count = PA_count + 1
                    PA_training_index.append(j)  
        elif validate_data.iloc[i]['protein'] == 'NP':
            for j in range(0, NP_raw.shape[0]):
                if difflib.SequenceMatcher(None, NP_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                #if NP_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    NP_count = NP_count + 1
                    NP_training_index.append(j)
        elif validate_data.iloc[i]['protein'] == 'PB1':
            for j in range(0, PB1_raw.shape[0]):
                if difflib.SequenceMatcher(None, PB1_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                #if PB1_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    PB1_count = PB1_count + 1
                    PB1_training_index.append(j)
        elif validate_data.iloc[i]['protein'] == 'PB2':
            for j in range(0, PB2_raw.shape[0]):
                if difflib.SequenceMatcher(None, PB2_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                #if PB2_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    PB2_count = PB2_count + 1
                    PB2_training_index.append(j)
        elif validate_data.iloc[i]['protein'] == 'NS1':
            for j in range(0, NS1_raw.shape[0]):
                if difflib.SequenceMatcher(None, NS1_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                #if NS1_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    NS1_count = NS1_count + 1
                    NS1_training_index.append(j)
#        elif validate_data.iloc[i]['protein'] == 'NS2':
#            for j in range(0, NS2_raw.shape[0]):
#                if difflib.SequenceMatcher(None, NS2_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
#                #if NS2_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
#                    NS2_count = NS2_count + 1
#                    NS2_training_index.append(j)
#        elif validate_data.iloc[i]['protein'] == 'M1':
#            for j in range(0, M1_raw.shape[0]):
#                if difflib.SequenceMatcher(None, M1_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
#                #if M1_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
#                    M1_count = M1_count + 1
#                    M1_training_index.append(j)
        elif validate_data.iloc[i]['protein'] == 'M2':
            for j in range(0, M2_raw.shape[0]):
                if difflib.SequenceMatcher(None, M2_raw.iloc[j]['seq'], validate_data.iloc[i]['seq']).quick_ratio() <= threshold:
                #if M2_raw.iloc[j]['seq'] == validate_data.iloc[i]['seq']:
                    M2_count = M2_count + 1
                    M2_training_index.append(j)
                    
                    
                    
    print(HA_count, NA_count, PA_count, NP_count, PB1_count, PB2_count, NS1_count, M2_count)
    HA_training_index = list(set(HA_training_index))
    NA_training_index = list(set(NA_training_index))
    PA_training_index = list(set(PA_training_index))
    NP_training_index = list(set(NP_training_index))
    PB1_training_index = list(set(PB1_training_index))
    PB2_training_index = list(set(PB2_training_index))
    NS1_training_index = list(set(NS1_training_index))
    #NS2_training_index = list(set(NS2_training_index))
    #M1_training_index = list(set(M1_training_index))
    M2_training_index = list(set(M2_training_index))
    print(len(HA_training_index), len(NA_training_index), len(PA_training_index), len(NP_training_index), len(PB1_training_index), len(PB2_training_index), len(NS1_training_index), len(M2_training_index))
    if len(HA_training_index) == 0:
        HA_training_index = list(range(0,HA_raw.shape[0]))
        
    if len(NA_training_index) == 0:
        NA_training_index = list(range(0,NA_raw.shape[0]))

    if len(PA_training_index) == 0:
        PA_training_index = list(range(0,PA_raw.shape[0]))
        
    if len(NP_training_index) == 0:
        NP_training_index = list(range(0,NP_raw.shape[0]))
        
    if len(PB1_training_index) == 0:
        PB1_training_index = list(range(0,PB1_raw.shape[0]))
        
    if len(PB2_training_index) == 0:
        PB2_training_index = list(range(0,PB2_raw.shape[0]))
        
    if len(NS1_training_index) == 0:
        NS1_training_index = list(range(0,NS1_raw.shape[0]))
        
#    if len(NS2_training_index) == 0:
#        NS2_training_index = list(range(0,NS2_raw.shape[0]))
#        
#    if len(M1_training_index) == 0:
#        M1_training_index = list(range(0,M1_raw.shape[0]))
        
    if len(M2_training_index) == 0:
        M2_training_index = list(range(0,M2_raw.shape[0]))
        
    return HA_training_index, NA_training_index, PA_training_index, NP_training_index, PB1_training_index, PB2_training_index, NS1_training_index, M2_training_index
    

#random select training dataset
def random_training_selection(input_index, segment_feature, segment_label):
    df_feature = []
    df_label = []
    #df_random_feature = []
    #df_random_label = []
    training_index = input_index
    seg_num = len(segment_label)

    for index in training_index:
        df_feature.append(segment_feature.iloc[index])
        df_label.append(segment_label.iloc[index])
    df_feature = pd.DataFrame(df_feature)
    df_feature.index = range(len(df_feature))
    df_label = pd.DataFrame(df_label)
    df_label.index = range(len(df_label))
    training_num = len(df_feature)
    train_ratio = round(float(seg_num)/(float(training_num*2)), 2)
    #print(len(df_feature))
    #print(len(training_index))
    if training_num <= int(seg_num/2):
        #model = 'error'
        model = randomforest_cross_validation(df_feature, df_label)
    else:
        df_random_feature, df_random_label = split_train_test(df_feature, df_label, train_ratio)
        print(df_random_feature.shape[0])
        model = randomforest_cross_validation(df_random_feature, df_random_label)
        #model = "error"
    return model
    

        
def genome_training_model(validate_data, threshold):
    HA_training_index, NA_training_index, PA_training_index, NP_training_index, PB1_training_index, PB2_training_index, NS1_training_index, M2_training_index =  training_data_preparation_version_2(validate_data, threshold)
    HA = random_training_selection(HA_training_index, HA_feature, HA_label)
    NA = random_training_selection(NA_training_index, NA_feature, NA_label)
    PA = random_training_selection(PA_training_index, PA_feature, PA_label)
    NP = random_training_selection(NP_training_index, NP_feature, NP_label)
    PB1 = random_training_selection(PB1_training_index, PB1_feature, PB1_label)
    PB2 = random_training_selection(PB2_training_index, PB2_feature, PB2_label)
    NS1 = random_training_selection(NS1_training_index, NS1_feature, NS1_label)
    #NS2 = random_training_selection(NS2_training_index, NS2_feature, NS2_label)
    #M1 = random_training_selection(M1_training_index, M1_feature, M1_label)
    M2 = random_training_selection(M2_training_index, M2_feature, M2_label)
    return HA, NA, PA, NP, PB1, PB2, NS1, M2
    
    






