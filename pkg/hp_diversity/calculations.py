import numpy as np
from scipy.special import xlogy

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

def Cq(A, beta, q, N):
    if q != 1:
        cq = (((1 / beta) ** (q - 1)) - ((1 / N) ** (q - 1))) / (1 - (1 / N) ** (q - 1))
    else:
        cq = 1 - (np.log(1 / beta)/np.log(1 / N))
    return cq

def Uq(A, beta, q, N):
    if q != 1:
        uq = (((1 / beta) ** (1 - q)) - ((1 / N) ** (1 - q))) / (1 - (1 / N) ** (1 - q))
    else:
        uq = 1 - (np.log(1 / beta)/np.log(1 / N))    
    return uq

def Sq(A, beta, q, N):
    sq = ((1 / beta) - (1 / N)) / (1 - 1 / N)
    return sq

def Vq(A, beta, q, N):
    vq = 1 - ((beta - 1) / (N - 1))
    return vq
