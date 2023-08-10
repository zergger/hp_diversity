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


# Main process 
args = "mother_mito_cheek.txt"
df = pd.read_csv(args, sep="\t")


# In[4]:


# Calculating
Alpha = []
grouped_data = df.groupby('ID')
for group_name, group_data in grouped_data:
    group_array = np.nan_to_num(group_data.iloc[:, 1:].astype(float).values)
    result = d_chao(group_array, "alpha", 1)
    Alpha.append([group_name] + [result])
Alpha = np.array(Alpha)


# In[5]:


#colnames = ["ID", "q=0", "q=1", "q=2", "q=3", "q=4"]
colnames = ["ID", "q=1"]
Alpha = pd.DataFrame(Alpha, columns=colnames)
Alpha.to_csv("mother_mito_cheek_alpha111.txt", sep="\t", index=False)


# In[ ]:




