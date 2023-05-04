import copy
import numpy as np
import math
import ipdb

def GIP_sm(A):
    # 求权重Y_sm-行
    w1 = np.linalg.norm(A,axis=1)
    widthSum = np.sum(np.square(w1))#该值为664
    Y_sm = A.shape[0]/widthSum#该值为831/664
    #构造GIP核相似性矩阵
    G = np.zeros((A.shape[0],A.shape[0]))
    for i in range(G.shape[0]):
        for j in range(i, G.shape[1]):
            G[i, j] = math.exp((-Y_sm)*np.square(np.linalg.norm(A[i] - A[j])))
            G[j, i] = G[i, j]
    return G

def GIP_m(A):
    # 求权重Y_m-列
    w1 = np.linalg.norm(A,axis=0)
    widthSum = np.sum(np.square(w1))
    Y_m = A.shape[1]/widthSum
    #构造GIP核相似性矩阵
    G = np.zeros((A.shape[1],A.shape[1]))
    for i in range(G.shape[0]):
        for j in range(i, G.shape[1]):
            G[i, j] = math.exp((-Y_m)*np.square(np.linalg.norm(A[:,i] - A[:,j])))
            G[j, i] = G[i, j]
    return G


def InSm(sm1, sm2, w):
    simm = w * sm1 + (1-w) * sm2
    return simm

def CaT_GIP(w):
    # A = np.loadtxt(r"D:\A论文编辑\数据集及代码\dataset1\SM-miRNA关联矩阵.txt", dtype=int)
    # SM = np.loadtxt(r"D:\A论文编辑\数据集及代码\dataset1\SM相似性矩阵.txt",dtype=float)
    # miRNA = np.loadtxt(r"D:\A论文编辑\数据集及代码\dataset1\miRNA相似性矩阵.txt",dtype=float)
    A = np.loadtxt(r"D:\A论文编辑\数据集及代码\dataset2\SM-miRNA association matrix 2.txt", dtype=int)
    SM = np.loadtxt(r"D:\A论文编辑\数据集及代码\dataset2\SM相似性矩阵2.txt",dtype=float)
    miRNA = np.loadtxt(r"D:\A论文编辑\数据集及代码\dataset2\miRNA相似性矩阵2.txt",dtype=float)
    B = GIP_sm(A)
    C = GIP_m(A)
    #相似性融合
    SM = InSm(copy.deepcopy(SM),copy.deepcopy(B),w)
    miRNA = InSm(copy.deepcopy(miRNA),copy.deepcopy(C),w)
    #相似性保存
    # np.savetxt(r'D:\A论文编辑\数据集及代码\dataset1\SM-GIP相似性矩阵.txt',SM,fmt='%e')
    # np.savetxt(r'D:\A论文编辑\数据集及代码\dataset1\miRNA-GIP相似性矩阵.txt',miRNA,fmt='%e')
    np.savetxt(r'D:\A论文编辑\数据集及代码\dataset2\SM-GIP相似性矩阵2.txt',SM,fmt='%e')
    np.savetxt(r'D:\A论文编辑\数据集及代码\dataset2\miRNA-GIP相似性矩阵2.txt',miRNA,fmt='%e')
    #矩阵拼接
    hs1 = np.hstack((SM,A))
    hs2 = np.hstack((np.transpose(A),miRNA))
    vs1 = np.vstack((hs1,hs2))
    #拼接矩阵保存
    # np.savetxt(r'D:\A论文编辑\数据集及代码\dataset1\T-GIP1.txt', vs1, fmt='%e')
    np.savetxt(r'D:\A论文编辑\数据集及代码\dataset2\T-GIP2.txt',vs1,fmt='%e')


