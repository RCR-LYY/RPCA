# -*- codeing = utf-8 -*-
# @Time : 2023/3/11 16:20
# @Author : Chuanru Ren
# @File : GIPK.py
# @Software: PyCharm

import copy
import numpy as np
import math
import ipdb

def GIP_sm(A):
    w1 = np.linalg.norm(A,axis=1)
    widthSum = np.sum(np.square(w1))
    Y_sm = A.shape[0]/widthSum
    G = np.zeros((A.shape[0],A.shape[0]))
    for i in range(G.shape[0]):
        for j in range(i, G.shape[1]):
            G[i, j] = math.exp((-Y_sm)*np.square(np.linalg.norm(A[i] - A[j])))
            G[j, i] = G[i, j]
    return G

def GIP_m(A):
    w1 = np.linalg.norm(A,axis=0)
    widthSum = np.sum(np.square(w1))
    Y_m = A.shape[1]/widthSum
    G = np.zeros((A.shape[1],A.shape[1]))
    for i in range(G.shape[0]):
        for j in range(i, G.shape[1]):
            G[i, j] = math.exp((-Y_m)*np.square(np.linalg.norm(A[:,i] - A[:,j])))
            G[j, i] = G[i, j]
    return G


def InSm(sm1, sm2, w):
    simm = w * sm1 + (1-w) * sm2
    return simm

def CaT_GIP():
    w = 0.2
    # Load Dataset 1
    M = np.loadtxt(r"SM-miRNA association matrix.txt", dtype=int)
    SM = np.loadtxt(r"SM similarity matrix.txt",dtype=float)
    miRNA = np.loadtxt(r"miRNA similarity matrix.txt",dtype=float)
    # Calculate the GIPK similarity
    GIPK_sm = GIP_sm(A)
    GIPK_m = GIP_m(A)
    # Similarity integration
    SM_new = InSm(copy.deepcopy(SM),copy.deepcopy(GIPK_sm),w)
    miRNA_new = InSm(copy.deepcopy(miRNA),copy.deepcopy(GIPK_m),w)
    # Construct the adjacency matrix of the heterogeneous SM-miRNA network
    hs1 = np.hstack((SM,A))
    hs2 = np.hstack((np.transpose(A),miRNA))
    vs1 = np.vstack((hs1,hs2))
    # Save data
    np.savetxt(r'T-GIP.txt', vs1, fmt='%e')


