#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Maria Tsontaki

"""
# Import Libraries
import re
import collections
import numpy as np
import pandas as pd

# open file and save it in a variable
fasta=""
with open("blast_homologs_O76756_edited_alignment.fst") as f:
    line=" ".join(line.strip() for line in f)
    fasta=fasta+line

# get all sequences
allSequencesList=[]
sequences=re.findall(r' .{1,250} ',fasta)
for seq in sequences:
    seqN=seq.replace(' ','')
    allSequencesList.append(seqN)
    
print("# of aligned sequences: {}".format(len(allSequencesList)))
    
# get a list of the a.a in every column 
positionsAAlist=[]

for pos in range(0,180):
    posAAList=[]
    for seq in allSequencesList:
        posAAList.append(seq[pos])
    positionsAAlist.append(posAAList)
    
            
# getting a.a frequencies and calculate Shannon Frequency
allEntropy=[]
for position in positionsAAlist:
    H=0
    elements_count = collections.Counter(position)
    for key, value in elements_count.items():
        keyFreq=value/29
        H+=-keyFreq*np.log2(keyFreq)
    allEntropy.append(H)
       
# have a look on the results
print(allEntropy)

# check indicated conserved areas 
pos=0 
for entro in allEntropy:
    pos+=1
    if entro>0:
        print("position: {} \nH: {} \n \n".format(pos,entro))
    
    
    
#-----------Save output------------#
df=pd.DataFrame(columns=["position","Shannon Entropy"])
pos=0 
for entro in allEntropy:
    pos+=1
    df=df.append({"position":pos,"Shannon Entropy":entro},ignore_index=True)
    
# save them in a .csv file
df.to_csv(r'./conservation_results.csv',index=False)
#---------------------------------#
  
    
    
    
    
    
    
    
    
        


