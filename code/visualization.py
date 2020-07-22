from __future__ import division
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 14:34:44 2019

@author: yinr0002
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Jul 14 20:07:21 2018

@author: Rayin
"""


import os, sys
import pandas as pd
import numpy as np
import cPickle
import matplotlib.pyplot as plt
import warnings
import random
import gc
import re
import math
import difflib
import matplotlib.ticker as ticker
    
from sklearn import ensemble
from sklearn import metrics
from matplotlib import pyplot 
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import train_test_split
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import roc_curve, auc
from imblearn.metrics import geometric_mean_score
#from sklearn.model_selection import StratifiedKFold 
from scipy import interp
from sklearn.multiclass import OneVsRestClassifier
from Bio import SeqIO
from PyBioMed.PyProtein import CTD



genome_sample = { "PB2":"MERIKELRDLMSQSRTREILTKTTVDHMAIIKKYTSGRQEKNPALRMKWMMAMRYPITADKRIMDMIPERNEQGQTLWSKTNDAGSDRVMVSPLAVTWWNRNGPTTSTVHYPKVYKTYFEKVERLKHGTFGPVHFRNQVKIRRRVDTNPGHADLSAKEAQDVIMEVVFPNEVGARILTSESQLAITKEKKEELQDCKIAPLMVAYMLERELVRKTRFLPVAGGTGSVYIEVLHLTQGTCWEQMYTPGGEVRNDDVDQSLIIAARNIVRRAAVSADPLASLLEMCHSTQIGGVRMVDILRQNPTEEQAVDICKAAIGLRISSSFSFGGFTFKRTSGSSVKKEEEVLTGNLQTLKIRVHEGYEEFTMVGRRATAILRKATRRLIQLIVSGRDEQSIAEAIIVAMVFSQEDCMIKAVRGDLNFVNRANQRLNPMHQLLRHFQKDAKVLFQNWGIESIDNVMGMIGILPDMTPSTEMSLRGIRVSKMGVDEYSSTERVVVSIDRFLRVRDQRGNVLLSPEEVSETQGTEKLTITYSSSMMWEINGPESVLVNTYQWIIRNWEIVKIQWSQDPTMLYNKMEFEPFQSLVPKATRSRYSGFVRTLFQQMRDVLGTFDTVQIIKLLPFAAAPPEQSRMQFSSLTVNVRGSGLRILVRGNSPVFNYNKATKRLTVLGKDAGALTEDPDEGTSGVESAVLRGFLILGKEDKRYGPALSINELSNLAKGEKANVLIGQGDVVLVMKRKRDSSILTDSQTATKRIRMAIN",

"PB1":"MDVNPTLLFLKIPAQNAISTTFPYTGDPPYSHGTGTGYTMDTVNRTHQYSEKGKWTTNTETGAPQLNPIDGPLPEDNEPSGYAQTDCVLEAMAFLEESHPGIFENSCLETMEVVQQTRVDKLTQGRQTYDWTLNRNQPAATALANTIEVFRSNGLTANESGRLIDFLKDVMESMNKEEIEITTHFQRKRRVRDNMTKKMVTQRTIGKKKQRLNKRGYLIRALTLNTMTKDAERGKLKRRAIATPGMQIRGFVYFVETLARSICEKLEQSGLPVGGNEKKAKLANVVRKMMTNSQDTEISFTITGDNTKWNENQNPRMFLAMITYITRNQPEWFRNILSMAPIMFSNKMARLGKGYMFESKRMKIRTQIPAEMLASIDLKYFNESTKKKIEKIRPLLIDGTASLSPGMMMGMFNMLSTVLGVSVLNLGQKKYTKTIYWWDGLQSSDDFALIVNAPNHEGIQAGVDRFYRTCKLVGINMSKKKSYINKTGTFEFTSFFYRYGFVANFSMELPSFGVSGVNESADMSIGVTVIKNNMINNDLGPATAQMALQLFIKDYRYTYRCHRGDTQIQTRRSFELKKLWDQTQSKVGLLVSDGGPNLYNIRNLHIPEVCLKWELMDDDYRGRLCNPLNPFVSHKEIDSVNNAVVMPAHGPAKSMEYDAVATTHSWIPKRNRSILNTSQRGILEDEQMYQKCCNLFEKFFPSSSYRRPVGISSMVEAMVSRARIDARVDFESGRIKKEEFSEIMKICSTIEELRRQK",

"PA":"MEDFVRQCFNPMIVELAEKAMKEYGEDPKIETNKFAAICTHLEVCFMYSDFHFIDERGESIIVESGDPNALLKHRFEIIEGRDRIMAWTVVNSICNTTGVEKPKFLPDLYDYKENRFIEIGVTRREVHIYYLEKANKIKSEKTHIHIFSFTGEEMATKADYTLDEESRARIKTRLFTIRQEMASRSLWDSFRQSERGEETIEEKFEITGTMRKLADQSLPPNFSSLENFRAYVDGFEPNGCIEGKLSQMSKEVNAKIEPFLRTTPRPLRLPDGPLCHQRSKFLLMDALKLSIEDPSHEGEGIPLYDAIKCMKTFFGWKEPNIVKPHEKGINPNYLMAWKQVLAELQDIENEEKIPRTKNMKRTSQLKWALGENMAPEKVDFDDCKDVGDLKQYDSDEPEPRSLASWVQNEFNKACELTDSSWIELDEIGEDVAPIEHIASMRRNYFTAEVSHCRATEYIMKGVYINTALLNASCAAMDDFQLIPMISKCRTKEGRRKTNLYGFIIKGRSHLRNDTDVVNFVSMEFSLTDPRLEPHKWEKYCVLEIGDMLLRTAIGQVSRPMFLYVRTNGTSKIKMKWGMEMRRCLLQSLQQIESMIEAESSVKEKDMTKEFFENKSETWPIGESPRGVEEGSIGKVCRTLLAKSVFNSLYASPQLEGFSAESRKLLLIVQALRDNLEPGTFDLGGLYEAIEECLINDPWVLLNASWFNSFLTHALK",

"HA":"MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETSSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADAYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNVPSVQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI",

"NP":"MASQGTKRSYEQMETGGERQDATEIRASVGRMIGGIGRFYIQMCTELKLSDYDGRLIQNSITIERMVLSAFDERRNKYLEEHPSAGKDPKKTGGPIYRRIDGKWMRELILYDKEEIRRVWRQANNGEDATAGLTHIMIWHSNLNDATYQRTRALVRTGMDPRMCSLMQGSTLPRRSGAAGAAVKGVGTIAMELIRMIKRGINDRNFWRGENGRRTRVAYERMCNILKGKFQTAAQRAMMDQVRESRNPGNAEIEDLIFLARSALILRGSVAHKSCLPACVYGLAVASGHDFEREGYSLVGIDPFKLLQNSQVVSLMRPNENPAHKSQLVWMACHSAAFEDLRVSSFIRGKKVIPRGKLSTRGVQIASNENVETMDSNTLELRSRYWAIRTRSGGNTNQQKASAGQISVQPTFSVQRNLPFERATVMAAFSGNNEGRTSDMRTEVIRMMESAKPEDLSFQGRGVFELSDEKATNPIVPSFDMSNEGSYFFGDNAEEYDS",

"NA":"MNPNQKIITIGSVCMTIGMANLILQIGNIISIWISHSIQLGNQNQIETCNQSVITYENNTWVNQTYVNISNTNFAAGQSVVSVKLAGNSSLCPVSGWAIYSKDNSIRIGSKGDVFVIREPFISCSPLECRTFFLTQGALLNDKHSNGTIKDRSPYRTLMSCPIGEVPSPYNSRFESVAWSASACHDGINWLTIGISGPDNGAVAVIKYNGIITDTIKSWRNNILRTQESECACVNGSCFTVMTDGPSDGQASYKIFRIEKGKIVKSVEMNAPNYHYEECSCYPDSSEITCVCRDNWHGSNRPWVSFNQNLEYQIGYICSGIFGDNPRPNDKTGSCGPVSSNGANGVKGFSFKYGNGVWIGRTKSISSRNGFEMIWDPNGWTGTDNNFSIKQDIVGINEWSGYSGSFVQHPELTGLDCIRPCFWVELIRGRPKENTIWTSGSSISFCGVNSDTVGWSWPDGAELPFTIDK",    

"M1":"MSLLTEVETYVLSIIPSGPLKAEIAQRLESVFAGKNTDLEALMEWLKTRPILSPLTKGILGFVFTLTVPSERGLQRRRFVQNALNGNGDPNNMDRAVKLYKKLKREITFHGAKEVSLSYSTGALASCMGLIYNRMGTVTTEAAFGLVCATCEQIADSQHRSHRQMATTTNPLIRHENRMVLASTTAKAMEQMAGSSEQAAEAMEVANQTRQMVHAMRTIGTHPSSSAGLKDDLLENLQAYQKRMGVQMQRFK",  

"M2":"MSLLTEVETPTRSEWECRCSDSSDPLVIAANIIGILHLILWITDRLFFKCIYRRFKYGLKRGPSTEGVPESMREEYQQEQQSAVDVDDGHFVNIELE",

"NS1":"MDSNTMSSFQVDCFLWHIRKRFADNGLGDAPFLDRLRRDQKSLKGRGNTLGLDIETATLVGKQIVEWILKEESSETLRMTIASVPTSRYLSDMTLEEMSRDWFMLMPRQKIIGPLCVRLDQAIMEKNIVLKANFSVIFNRLETLILLRAFTEKGAIVGEISPLPSLPGHTYEDVKNAVGVLIGGLEWNGNTVRVSENIQRFAWRNCDENGRPSLPPEQK",

"NS2":"MDSNTMSSFQDILMRMSKMQLGSSSEDLNGMVTRFESLKIYRDSLGETVMRMGDLHYLQSRNEKWREQLGQKFEEIRWLIEEMRHRLKATENSFEQITFMQALQLLLEVEQEIRAFSFQLI"
}


def split_year(input_file):
    df = input_file
    df_year = []
    df_description = df["description"]
    df_year = df["year"]
    for row in range(0, df.shape[0]):
        #df_year[row] = re.findall(r'(\d{4})\s.*', df_description[row])
        df_year[row] = re.findall(r'\s(\d{4})', df_description[row])
        #if df_year[row] == '':
            #df_year[row] = re.findall(r'\s(\d{4})', df_description[row])
        df_year[row] = ''.join(df_year[row])   
    df["year"] = df_year
    return df

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
  
#this is the function that split the description
def split_description(input_file):
    df = input_file
    df_description = df["description"]
    for row in range(0, df.shape[0]):
        try:
            df_description[row] = re.search(r'\bA/\S.*', df_description[row]).group()
        except AttributeError:
            df_description[row] = re.search(r'\bA/\S.*', df_description[row])
            print row
        #df_description[row] = re.search(r'\bA/\S.*', df_description[row]).group()
        #df_description[row] = re.search(r'\bA/\S.*', df_description[row])
    df["description"] = df_description
    return df


    
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

os.chdir("C:/Users/yinr0002/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data")    
#os.chdir("C:/Users/Rayin/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data")
#os.chdir("/Users/rayin/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data")
#########################################################################################################################################
#read in the sequence data of full length for each segment after feature transformation
#global HA_raw, HA_data, HA_feature, HA_label
HA_raw = pd.read_csv('sequence/csv_update/HA_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
HA_data = pd.read_csv('sequence/train/full_length/HA_data.csv')
HA_feature = HA_data.ix[:, 1:148]
HA_label = HA_data.ix[:, 149]
HA_raw['year'] = 0
HA_raw = split_year(HA_raw)

#global NA_raw, NA_data, NA_feature, NA_label
NA_raw = pd.read_csv('sequence/csv_update/NA_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
NA_data = pd.read_csv('sequence/train/full_length/NA_data.csv')
NA_feature = NA_data.ix[:, 1:148]
NA_label = NA_data.ix[:, 149]
NA_raw['year'] = 0
NA_raw = split_year(NA_raw)

#global NP_raw, NP_data, NP_feature, NP_label
NP_raw = pd.read_csv('sequence/csv_update/NP_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
NP_data = pd.read_csv('sequence/train/full_length/NP_data.csv')
NP_feature = NP_data.ix[:, 1:148]
NP_label = NP_data.ix[:, 149]
NP_raw['year'] = 0
NP_raw = split_year(NP_raw)

#global PA_raw, PA_data, PA_feature, PA_label
PA_raw = pd.read_csv('sequence/csv_update/PA_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
PA_data = pd.read_csv('sequence/train/full_length/PA_data.csv')
PA_feature = PA_data.ix[:, 1:148]
PA_label = PA_data.ix[:, 149]
PA_raw['year'] = 0
PA_raw = split_year(PA_raw)

#global PB1_raw, PB1_data, PB1_feature, PB1_label
PB1_raw = pd.read_csv('sequence/csv_update/PB1_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
PB1_data = pd.read_csv('sequence/train/full_length/PB1_data.csv')
PB1_feature = PB1_data.ix[:, 1:148]
PB1_label = PB1_data.ix[:, 149]
PB1_raw['year'] = 0
PB1_raw = split_year(PB1_raw)

#global PB2_raw, PB2_data, PB2_feature, PB2_label
PB2_raw = pd.read_csv('sequence/csv_update/PB2_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
PB2_data = pd.read_csv('sequence/train/full_length/PB2_data.csv')
PB2_feature = PB2_data.ix[:, 1:148]
PB2_label = PB2_data.ix[:, 149]
PB2_raw['year'] = 0
PB2_raw = split_year(PB2_raw)


#global M2_raw, M2_data, M2_feature, M2_label
M2_raw = pd.read_csv('sequence/csv_update/M2_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
M2_data = pd.read_csv('sequence/train/full_length/M2_data.csv')
M2_feature = M2_data.ix[:, 1:148]
M2_label = M2_data.ix[:, 149]
M2_raw['year'] = 0
M2_raw = split_year(M2_raw)

#global NS1_raw, NS1_data, NS1_feature, NS1_label
NS1_raw = pd.read_csv('sequence/csv_update/NS1_raw.csv', header=None, na_filter = False, names=['accession', 'description', 'protein', 'seq'])
NS1_data = pd.read_csv('sequence/train/full_length/NS1_data.csv')
NS1_feature = NS1_data.ix[:, 1:148]
NS1_label = NS1_data.ix[:, 149]
NS1_raw['year'] = 0
NS1_raw = split_year(NS1_raw)




os.chdir("C:/Users/yinr0002/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data")
genome_avian = pd.read_csv("Genome/host_type/csv/avian.csv", names=['accession', 'description', 'protein', 'seq'], na_filter = False)
genome_avian = fill_with_missing_protein(genome_avian)
genome_avian = split_description(genome_avian)
genome_avian['year'] = 0
genome_avian = split_year(genome_avian) 

genome_human = pd.read_csv("Genome/host_type/csv/human.csv", names=['accession', 'description', 'protein', 'seq'], na_filter = False)
genome_human = fill_with_missing_protein(genome_human)
genome_human = split_description(genome_human)
genome_human['year'] = 0
genome_human = split_year(genome_human) 

genome_swine = pd.read_csv("Genome/host_type/csv/swine.csv", names=['accession', 'description', 'protein', 'seq'], na_filter = False)
genome_swine = fill_with_missing_protein(genome_swine)
genome_swine = split_description(genome_swine)
genome_swine['year'] = 0
genome_swine = split_year(genome_swine) 
           
        
def randomforest_cross_validation(train_x, train_y):
    np.random.seed(100)
    clf = ensemble.RandomForestClassifier()
    
    clf.fit(train_x, train_y)
    
    #calculate the accuracy
    accuracy = cross_val_score(clf, train_x, train_y, cv=10, scoring='accuracy') 
    print "accuracy: %f" %accuracy.mean() + '\n'
    #print accuracy
    
    #calculate the precision
    precision = cross_val_score(clf, train_x, train_y, cv=10, scoring='precision_macro')
    print "precision: %f" %precision.mean() + '\n'
    
    #calculate the recall score
    recall = cross_val_score(clf, train_x, train_y, cv=10, scoring='recall_macro')
    print "recall: %f" %recall.mean() + '\n'
    
    #calculate the f_measure
    f_measure = cross_val_score(clf, train_x, train_y, cv=10, scoring='f1_macro')
    print "f_measure: %f " %f_measure.mean() + '\n' 

    #generate classification report and MCC and G-mean value
    y_pred = cross_val_predict(clf, train_x, train_y, cv=10)
    G_mean = geometric_mean_score(train_y, y_pred)
    MCC = matthews_corrcoef(train_y, y_pred)
    print "G_mean: %f" %G_mean.mean() + '\n'
    print "MCC: %f" %np.mean(MCC) + '\n'
    
    print "Classification_report:"
    print(metrics.classification_report(train_y, y_pred))    
    
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
           


HA = randomforest_cross_validation(HA_feature, HA_label)
NA = randomforest_cross_validation(NA_feature, NA_label)
NP = randomforest_cross_validation(NP_feature, NP_label)
PA = randomforest_cross_validation(PA_feature, PA_label)
PB1 = randomforest_cross_validation(PB1_feature, PB1_label)
PB2 = randomforest_cross_validation(PB2_feature, PB2_label)
#M1 = randomforest_cross_validation(M1_feature, M1_label)
M2 = randomforest_cross_validation(M2_feature, M2_label)
NS1 = randomforest_cross_validation(NS1_feature, NS1_label)
#NS2 = randomforest_cross_validation(NS2_feature, NS2_label)  

         
#this is the genome prediction function with respect with csv file as input 
def prediction_genome_csv(csv_file, HA, NA, PA, NP, PB1, PB2, NS1, M2):
    prediction = {'description':[], 'HA':[], 'NA':[], 'NP':[], 'PA':[], 'PB1':[],'PB2':[],'M2':[], 'NS1':[]}
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
        print x
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
                    print "HA\t\t", HA_prediction, HA_prob
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
                    print "NA\t\t", NA_prediction, NA_prob
            
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
                    print "NP\t\t", NP_prediction, NP_prob
            
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
                    print "PA\t\t", PA_prediction, PA_prob  
            
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
                    print "PB1\t\t", PB1_prediction, PB1_prob 
            
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
                    print "PB2\t\t", PB2_prediction, PB2_prob 
                    
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
                    print "M2\t\t", M2_prediction, M2_prob 
             
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
                    print "NS1\t\t", NS1_prediction, NS1_prob 
        
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

    return prediction, avian_prob, human_prob, swine_prob


    
#calculate the probability of reassortment of genome set
def calculate_reassort_prob(avian, human, swine):
    avian = avian
    human = human
    swine = swine
    row = avian.shape[0]
    columns = avian.columns[:-1]
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
            if avian[column].ix[i] == 'missing':
                avian[column].ix[i] = 1.0
            avian_initial = avian_initial * avian[column].ix[i]
            if human[column].ix[i] == 'missing':
                human[column].ix[i] = 1.0
            human_initial = human_initial * human[column].ix[i]
            if swine[column].ix[i] == 'missing':
                swine[column].ix[i] = 1.0
            swine_initial = swine_initial * swine[column].ix[i]
        reassort_prob = 1 - avian_initial - human_initial - swine_initial
        avian_prob.append(avian[column].ix[i])
        human_prob.append(human[column].ix[i])
        swine_prob.append(swine[column].ix[i])
        reassortment_prob.append(reassort_prob)    
        reassortment_result['description'].append(avian.ix[i].description)
        reassortment_result['reassort_prob'].append(reassort_prob) 

    return reassortment_result    


os.chdir("C:/Users/yinr0002/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data/Genome/all_year/reassortment_prob_host")

#data preparation for the genome in different years
def genome_extract_by_year(year):
    df = genome_avian
    df_year = df[df['year'] == str(year)]
    return df_year

def calculate_reassort_prob_by_year(year):
    genome_year = genome_extract_by_year(year)
    if genome_year.empty:
        print 'Empty in genome year %d' %year
    else:     
        genome_year = genome_year.drop(['year'], axis=1)
        filepath_prob = str('reassortment_prob') + '_' + str(year)+ '.csv'
        filepath_pred = str('prediction') + '_' + str(year)+ '.csv'
        #HA, NA, PA, NP, PB1, PB2, NS1, NS2, M1, M2 = genome_training_model_with_year(year)
        prediction, pred_avian, pred_human, pred_swine = prediction_genome_csv(genome_year, HA, NA, PA, NP, PB1, PB2, NS1, M2) 
        pred_df = pd.DataFrame(prediction)
        pred_avian_df = pd.DataFrame(pred_avian)
        pred_human_df = pd.DataFrame(pred_human)
        pred_swine_df = pd.DataFrame(pred_swine)
        reassortment_result = calculate_reassort_prob(pred_avian_df, pred_human_df, pred_swine_df)
        reassortment_result = pd.DataFrame(reassortment_result)  
        reassortment_result.to_csv(filepath_prob)
        pred_df.to_csv(filepath_pred)
        print 'Finish reassortment probability calculation in year %d' %year
        
        count_reassort = 0
        count_nonreassort = 0
        for i in range(0, reassortment_result.shape[0]):
            if reassortment_result.ix[i]['reassort_prob'] >= 0.5:
                count_reassort = count_reassort + 1
        count_nonreassort = reassortment_result.shape[0] - count_reassort
        accuracy = count_reassort / reassortment_result.shape[0]
        return count_reassort, count_nonreassort, accuracy

        
reassort_result_by_year = {'year':[], 'count_reassort':[], 'count_nonreassort':[], 'accuracy':[]}

for year in range(1918, 2018):
    genome_year = genome_extract_by_year(year)
    if genome_year.empty:
        print "empty"
    else:
        count_reassorted, count_nonreassorted, acc = calculate_reassort_prob_by_year(year)
        reassort_result_by_year['year'].append(year)
        reassort_result_by_year['count_reassort'].append(count_reassorted)
        reassort_result_by_year['count_nonreassort'].append(count_nonreassorted)
        reassort_result_by_year['accuracy'].append(acc)

reassort_result_by_year = pd.DataFrame(reassort_result_by_year)
reassort_result_by_year.to_csv("avian/reassort_avian_by_year.csv")

################################################################################################################################
#draw the figure of reassortment occurence for avian, human and swine genomes
os.chdir('C:/Users/yinr0002/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data/Genome/all_year/reassortment_prob_host/')

df = pd.read_csv('avian/reassort_avian_by_year.csv', delimiter=',', na_filter = False)
df = df.drop(df.columns[0], axis=1)
avail_df = []

for i in range(0, df.shape[0]):
    if df['count_nonreassort'].ix[i] + df['count_reassort'].ix[i] > 10:
        df['year'].ix[i] = str(df['year'].ix[i])
        avail_df.append(df.ix[i])
avail_df = pd.DataFrame(avail_df)  

#plot the figure present ratio between non-reassorted and reassorted strains across years
y = avail_df['accuracy']
x = list(avail_df['year'])
plt.bar(x, y)

def plot_reassortment_rate():
    fig, ax = plt.subplots()
    ax.plot(x, y, color='red', marker='o', linestyle='dashed', linewidth=2, markersize=6)
    
    plt.xlabel('year', fontsize=10)
    plt.ylabel('ratio', fontsize=10)
    #plt.title('The ratio between non-reassorted and reassorted strains in' + '\n' + 'different years detected by our model', fontsize=10)
    plt.xticks([1957, 1968, 1977, 2009])
    
    plt.savefig('C:/Users/yinr0002/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/results/avian_reassortment_rate.eps', format='eps', dpi=1000)
    plt.show()


fig, ax = plt.subplots()
ax.plot(x, y, color='red', marker='o', linestyle='dashed', linewidth=2, markersize=6)
formatter = ticker.FormatStrFormatter('%1.2f')
ax.yaxis.set_major_formatter(formatter)

# 
for tick in ax.yaxis.get_major_ticks():
    tick.label1On = True  # label1On 
    #tick.label2On = True  # label2On 
    tick.label1.set_color('black')
    #tick.label2.set_color('green')

##
#for line in ax.yaxis.get_ticklines():
#    # line is a Line2D instance
#    line.set_color('green')
#    line.set_markersize(25)
#    line.set_markeredgewidth(3)

#
for label in ax.xaxis.get_ticklabels():
    # label is a Text instance
    label.set_color('black')
    label.set_rotation(45)
    label.set_fontsize(10)
  
plt.show()
plt.savefig('C:/Users/yinr0002/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/results/reassortment_rate.eps', format='eps', dpi=300)


