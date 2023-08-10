#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import math
import sys
from scipy.special import xlogy


# In[2]:


def d_chao(A, lev, q): 
    tot = np.sum(A)
    if tot == 0:
        return 0 
    
    eA = A / tot

    cA = np.sum(A, axis=0)
    N = A.shape[0]

    ecA = cA / tot
  
    if lev == 'alpha':
        if q != 1:
            Da = (1 / N) * np.sum(eA**q)**(1 / (1 - q))
            D_value = Da
        else:
            Da = np.exp(-np.sum(xlogy(eA, eA)) - np.log(N))
            D_value = Da
    
    if lev == 'beta':
        D_value = d_chao(A, lev='gamma', q=q) / d_chao(A, lev='alpha', q=q)
    
    if lev == 'gamma':
        if q != 1:
            Dg = np.sum(ecA**q)**(1 / (1 - q))
            D_value = Dg
        else:
            Dg = np.exp(-np.sum(xlogy(ecA, ecA)))
            D_value = Dg
    
    return D_value


# In[3]:


def Cq(A, beta, q, N):
    if q != 1:
        cq = (((1 / beta) ** (q - 1)) - ((1 / N) ** (q - 1))) / (1 - (1 / N) ** (q - 1))
    else:
        cq = 1 - (np.log(1 / beta)/np.log(1 / N))
    return cq


# In[4]:


def Uq(A, beta, q, N):
    if q != 1:
        uq = (((1 / beta) ** (1 - q)) - ((1 / N) ** (1 - q))) / (1 - (1 / N) ** (1 - q))
    else:
        uq = 1 - (np.log(1 / beta)/np.log(1 / N))    
    return uq


# In[5]:


def Sq(A, beta, q, N):
    sq = ((1 / beta) - (1 / N)) / (1 - 1 / N)
    return sq

def Vq(A, beta, q, N):
    vq = 1 - ((beta - 1) / (N - 1))
    return vq


# In[6]:


# Main process 
args = ["blood_pcg_child_overall.txt","blood_pcg_mother_overall.txt"]
df1 = pd.read_csv(args[0], sep="\t")
df2 = pd.read_csv(args[1], sep="\t")


# In[7]:


Res = pd.DataFrame(columns=["ID1", "ID2", "beta.D", "cq"])
visited_pairs = set()
grouped_data1 = df1.groupby('ID')
grouped_data2 = df2.groupby('ID')
for group_name1, group_data1 in grouped_data1:
    for group_name2, group_data2 in grouped_data2:
        #if group_name1 == group_name2:
            #continue
        if (group_name1, group_name2) in visited_pairs or (group_name2, group_name1) in visited_pairs:
            continue
        visited_pairs.add((group_name1, group_name2))            
        group_array1 = np.nan_to_num(group_data1.iloc[:, 1:].astype(float).values)
        group_array2 = np.nan_to_num(group_data2.iloc[:, 1:].astype(float).values)
        sda = np.vstack((group_array1, group_array2))
        beta = d_chao(sda, "beta", 1)
        cqn = Cq(sda, beta, 1, sda.shape[0]) 
        #uqn = Uq(sda, beta, 2, sda.shape[0]) 
        #sqn = Sq(sda, beta, 2, sda.shape[0]) 
        #vqn = Vq(sda, beta, 2, sda.shape[0]) 
        res = pd.DataFrame([[group_name1, group_name2, beta, cqn]], columns=["ID1", "ID2", "beta.D", "cq"])
        Res = pd.concat([Res, res], ignore_index=True)
Res.to_csv("Beta-diversity-and-similarity_blood_pcg_child_mother111.txt", sep="\t", index=False)


# In[ ]:




