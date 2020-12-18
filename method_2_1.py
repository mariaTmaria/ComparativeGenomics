#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Maria Tsontaki

"""

# Import Libraries 
import re
import pandas as pd

# Extract best pair hits and their score for a blast output file
def findMatchedPairs(file,selfDB):
    pairs=[]
    
    with open(file) as fp:
       line = fp.readline()
       foundSubject=0
       
       # Extract aligned pairs
       while line:
           line = fp.readline()
           
           # find protein 1 from query:
           query=re.findall(r'Query= [a-z]{2}\|[a-zA-Z]+.*\|',line)
           if len(query):
               foundSubject=1
               match1=re.findall(r'\|[a-zA-Z]+.*\|',str(query))
               if len(match1):
                   protein1=re.sub("[\[\]\'\|]",'',str(match1))
    
           # find protein 2 from subject DB and take the first match:
           dbFind=re.findall(r'^[a-z]{2}\|[a-zA-Z]+.*\|.*\s{2}[0-9]*\.{0,1}[0-9]*\s*\.{0,1}[0-9]*',line)
           if len(dbFind) and foundSubject!=0:
               if selfDB==False or (selfDB==True and foundSubject==2):
                   foundSubject=0
                   scoreFind=re.findall(r'  [0-9]+\.{0,1}[0-9]* ',dbFind[0])
                   score=float(re.sub("[\[\]\']",'',str(scoreFind[0])))
                   match2=re.findall(r'\|[a-zA-Z]+.*\|',dbFind[0])
        
                   if len(match2):
                       protein2=re.sub("[\[\]\'\|]",'',str(match2))
                       pairs.append([protein1,protein2,score])
               elif selfDB==True:
                   foundSubject=2
    
    return pairs

# Extract best pair hits and their score for every blast output file
pairs_sp1_sp2=findMatchedPairs('results_beeDB_waspQUERY.txt',False)
pairs_sp2_sp1=findMatchedPairs('results_waspDB_beeQUERY.txt',False)
pairs_sp1_sp1=findMatchedPairs('results_beeDB_beeQUERY.txt',True)
pairs_sp2_sp2=findMatchedPairs('results_waspDB_waspQUERY.txt',True)


# find bidirectional hits and calculate new score == initial orthologs
bbhInitial=[]
for pair1 in pairs_sp1_sp2:
    for pair2 in pairs_sp2_sp1:
        if pair1[0]==pair2[1] and pair1[1]==pair2[0]:
            bbhInitial.append([pair1[0],pair1[1],(pair1[2]+pair2[2])/2]) #0: wasp, 1:bee, 2:score

print("initial BBH: {}".format(len(bbhInitial)))
            
# find if same species genes are more similar than bidirectional hits => paralogs
bbhOrtho_only=bbhInitial

def get_paralogs_and_their_ortho(selfPairsList):
    paralogs_and_their_ortho=[]
    for pairsSelf in selfPairsList:
        for pairsBBH in bbhInitial:
            if pairsSelf[0]==pairsBBH[0] or pairsSelf[1]==pairsBBH[0]:
                if pairsSelf[2]>pairsBBH[2]:
                    paralogs_and_their_ortho.append([pairsSelf[0],pairsSelf[1],pairsBBH[1],pairsSelf[2],pairsBBH[2]])
                    bbhOrtho_only.remove([pairsBBH[0],pairsBBH[1],pairsBBH[2]])
            if pairsSelf[0]==pairsBBH[1] or pairsSelf[1]==pairsBBH[1]:
                if pairsSelf[2]>pairsBBH[2]:
                    paralogs_and_their_ortho.append([pairsSelf[0],pairsSelf[1],pairsBBH[0],pairsSelf[2],pairsBBH[2]])
                    bbhOrtho_only.remove([pairsBBH[0],pairsBBH[1],pairsBBH[2]])
    
    print("complete orthologs,paralogs pairs: {}".format(len(paralogs_and_their_ortho)))
    
    return paralogs_and_their_ortho

#get paralogs and their ortholog pairs for bee as DN and then wasp 
paralogs_and_their_ortho_BEE=get_paralogs_and_their_ortho(pairs_sp1_sp1) #245 triplets of 2 paralogs and their orholog found in bee
paralogs_and_their_ortho_WASP=get_paralogs_and_their_ortho(pairs_sp2_sp2) #359 triplets of 2 paralogs and their orholog found in wasp


print("only ortholog pairs: {}".format(len(bbhOrtho_only)))



#-----------Save output------------#
df=pd.DataFrame(columns=["ortholog/paralog1","paralog2","ortholog","score"])
for pair in bbhOrtho_only:
    df=df.append({"ortholog/paralog1":pair[0],"paralog2":"no","ortholog":pair[1],"score":pair[2]},ignore_index=True)
    
# for triplets of the 2 paralogs and 1 ortholog, as score will be indicated as the difference of the score of the BBH with the score of the self match
for pair in paralogs_and_their_ortho_BEE:
    df=df.append({"ortholog/paralog1":pair[0],"paralog2":pair[1],"ortholog":pair[2],"score":(pair[3]-pair[4])},ignore_index=True)

for pair in paralogs_and_their_ortho_WASP:
    df=df.append({"ortholog/paralog1":pair[0],"paralog2":pair[1],"ortholog":pair[2],"score":(pair[3]-pair[4])},ignore_index=True)

# save them in a .csv file
df.to_csv(r'./orthologs_found_bbh.csv',index=False)
#---------------------------------#































    
    
  