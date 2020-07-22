# -*- coding: utf-8 -*-
"""
Created on Wed Dec 06 15:05:27 2017

@author: Rayin
"""

import os, sys, csv
import pandas as pd
import numpy as np
from PyBioMed.PyProtein import CTD
from Bio import AlignIO
from Bio import SeqIO
#from pybiomed.PyProtein import CTD

os.chdir("C:/Users/Rayin/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data/Sequence")

#convert fasta file to csv file  
def fasta_to_csv(input_file):
    flu = open(os.path.expanduser("csv/seq.csv"), 'w')
    year = open(os.path.expanduser("csv/year.csv"), 'w')
    for seq_record in input_file:
        flu_seq = ''
        record_id = seq_record.id
        seq = seq_record.seq
        flu_seq = str(record_id) + ',' + str(seq) + "\n"
        flu.write(flu_seq)
        
        flu_description = ''
        flu_description = seq_record.description.split(' ')
        for i in range(0, len(flu_description)):
            if flu_description[i][0:2] == '19':
                seq_year = flu_description[i][0:4] + '\n'
            elif flu_description[i][0:2] == '20':
                seq_year = flu_description[i][0:4] + '\n'
            else:
                next
        year.write(seq_year)
    flu.close()
    year.close()
    
#read fasta file of all subtypes with full length sequence data
HA_segment = SeqIO.parse('raw/full_length/HA/HA_all.fa', 'fasta')
M1_segment = SeqIO.parse('raw/full_length/M1/M1_all.fa', 'fasta')
M2_segment = SeqIO.parse('raw/full_length/M2/M2_all.fa', 'fasta')
NA_segment = SeqIO.parse('raw/full_length/NA/NA_all.fa', 'fasta')
NP_segment = SeqIO.parse('raw/full_length/NP/NP_all.fa', 'fasta')
NS1_segment = SeqIO.parse('raw/full_length/NS1/NS1_all.fa', 'fasta')
NS2_segment = SeqIO.parse('raw/full_length/NS2/NS2_all.fa', 'fasta')
PA_segment = SeqIO.parse('raw/full_length/PA/PA_all.fa', 'fasta')
PB1_segment = SeqIO.parse('raw/full_length/PB1/PB1_all.fa', 'fasta')
PB1_F2_segment = SeqIO.parse('raw/full_length/PB1-F2/PB1_F2_all.fa', 'fasta')
PB2_segment = SeqIO.parse('raw/full_length/PB2/PB2_all.fa', 'fasta')


fasta_to_csv(HA_segment)
fasta_to_csv(M1_segment) 
fasta_to_csv(M2_segment)  
fasta_to_csv(NA_segment) 
fasta_to_csv(NP_segment) 
fasta_to_csv(NS1_segment)  
fasta_to_csv(NS2_segment)
fasta_to_csv(PA_segment)
fasta_to_csv(PB1_segment)
fasta_to_csv(PB1_F2_segment)
fasta_to_csv(PB2_segment)



#read fasta file of all subtypes with all length sequence data
os.chdir("C:/Users/Rayin/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data/Sequence")

HA_segment = SeqIO.parse('raw/all_length/HA/HA_all.fa', 'fasta')
M1_segment = SeqIO.parse('raw/all_length/M1/M1_all.fa', 'fasta')
M2_segment = SeqIO.parse('raw/all_length/M2/M2_all.fa', 'fasta')
NA_segment = SeqIO.parse('raw/all_length/NA/NA_all.fa', 'fasta')
NP_segment = SeqIO.parse('raw/all_length/NP/NP_all.fa', 'fasta')
NS1_segment = SeqIO.parse('raw/all_length/NS1/NS1_all.fa', 'fasta')
NS2_segment = SeqIO.parse('raw/all_length/NS2/NS2_all.fa', 'fasta')
PA_segment = SeqIO.parse('raw/all_length/PA/PA_all.fa', 'fasta')
PB1_segment = SeqIO.parse('raw/all_length/PB1/PB1_all.fa', 'fasta')
PB1_F2_segment = SeqIO.parse('raw/all_length/PB1-F2/PB1-F2_all.fa', 'fasta')
PB2_segment = SeqIO.parse('raw/all_length/PB2/PB2_all.fa', 'fasta')

fasta_to_csv(HA_segment)
fasta_to_csv(M1_segment) 
fasta_to_csv(M2_segment)  
fasta_to_csv(NA_segment) 
fasta_to_csv(NP_segment) 
fasta_to_csv(NS1_segment)  
fasta_to_csv(NS2_segment)
fasta_to_csv(PA_segment)
fasta_to_csv(PB1_segment)
fasta_to_csv(PB1_F2_segment)
fasta_to_csv(PB2_segment)


#read fasta file in a genome format for avian, human and swine
os.chdir("C:/Users/yinr0002/Google Drive/Tier_2_MOE2014/5_Journal/Bioinformatics/data/Genome")  
Genome_avian = SeqIO.parse("host_type/raw/avian.fa", 'fasta')
Genome_human = SeqIO.parse("host_type/raw/human.fa", 'fasta')
Genome_swine = SeqIO.parse("host_type/raw/swine.fa", 'fasta')

fasta_to_csv(Genome_avian)
fasta_to_csv(Genome_human) 
fasta_to_csv(Genome_swine) 


#for seq_record in SeqIO.parse(protein_file, 'fasta'):
#    record_id = seq_record.id
#    sequence = seq_record.seq
#    protein_csv = str(sequence) + '\n'
#    output_file.write(protein_csv)
#
#output_file.close()

#alignment = AlignIO.read('Subtype_H2.fasta', 'fasta')
#align_array = np.array([list(rec) for rec in alignment], np.character)  #alignments as array
#align_array = pd.DataFrame(align_array)
#align_array.to_csv('align_array.csv')
#seq_id = {}
#seq_mapping = {}
#for i, record in enumerate(alignment):
#    seq_id[i] = record.id
#    seq_mapping[i] = record.seq
#    seq = str(seq_id[i]) + ',' + str(seq_mapping[i]) + '\n'
#    output_file.write(seq)
#
#output_file.close()

        
            
#myfile = open(os.path.expanduser("C:/Users/Rayin/Google Drive/Tier_2_MOE2014/5_Journal/Plos_One/Data/mafft/H2_csv.csv"), 'w')
#
#H2_alignment = SeqIO.parse('Subtype_H2.fasta', 'fasta')
#for seq_record in H2_alignment:
#    H2_csv = ''
#    record_id = seq_record.id
#    seq =seq_record.seq
#    H2_csv = str(record_id) + "," + str(seq) + "\n"
#    myfile.write(H2_csv)
#myfile.close()      
    
        
       
        
## open file and iterate through the lines, composing each single line as we go
#out_lines = []
#temp_line = ''
#with open('Subtype_H2.fasta','r') as fp:
#     for line in fp:
#         if line.startswith('>'):
#             out_lines.append(temp_line)
#             temp_line = line.strip() + '\t'
#         else:
#             temp_line += line.strip()
#
#with open('Subtype_H2.csv', 'w') as fp_out:
#    fp_out.write('\n'.join(out_lines))    
    
    



 