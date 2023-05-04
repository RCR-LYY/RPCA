import numpy as np
import copy
import ipdb
import time

def DC(D,mu,T0,g):
    U,S,V = np.linalg.svd(D)
    T1 = np.zeros(np.size(T0))
    for i in range(1,100):
        T1 = DCInner(S,mu,T0,g)
        err = np.sum(np.square(T1-T0))#向量2范数的平方
        if err < 1e-6:
            break
        T0 = T1
    l_1 = np.dot(U, np.diag(T1))
    l = np.dot(l_1, V)
    return l,T1


def DCInner(S,mu,T_k,gam):
    lamb = 1/mu
    grad = (1+gam)*gam/(np.square(gam+T_k))
    T_k1 = S-lamb*grad
    T_k1[T_k1<0]=0
    return T_k1


def errorsol(Y_1,H_1,L_1,lamb,mu,type):
    if type ==1:
        D=-Y_1/mu+(H_1-L_1)
        E=np.zeros(np.shape(D))
        eps = lamb/mu
        DD = np.abs(D)-eps
        DD2= DD*np.sign(D)
        ID = np.abs(D)>eps
        E[ID]=DD2[ID]
    else:
        alpha = lamb /mu
        G = H_1 - L_1 - Y_1 / mu
        G1=np.sqrt(np.sum(np.square(G),1))
        G1[G1 == 0] = alpha#该行存在的意义是保证下一行做除法时不会因除数G1==0而出错，不会对最终结果造成影响。
        G2 = (G1 - alpha)/ G1
        E = np.dot(G ,np.diag((G1 > alpha)* G2))
    return E


def GAMA(H): # H = L+S
    # 可变参数
    muzero = 18
    mu = muzero
    # 固定参数
    gamma = 0.01
    type = 21
    lamb = 1e-3
    rho = 1.1
    tol = 1e-3
    # 迭代初始化
    m, n = np.shape(H)
    L = copy.deepcopy(H)
    S = np.zeros((m,n))
    Y = np.zeros((m,n))
    # 矩阵迭代更新
    for i in range(0, 500):
        D = H-S-Y/mu
        sig = np.zeros(min(m, n))
        L, sig = DC(copy.deepcopy(D),mu,copy.deepcopy(sig),gamma)
        S = errorsol(copy.deepcopy(Y),copy.deepcopy(H),copy.deepcopy(L),lamb,mu,type)
        Y= Y+mu*(L-H+S)
        mu = mu*rho
        sigma = np.linalg.norm(H-S-L,'fro')
        RRE = sigma/np.linalg.norm(H,'fro')
        # r = np.linalg.matrix_rank(L)
        # print(['Relative err', RRE, ' Rank-L ', r])
        if RRE < tol:
            break
    # ipdb.set_trace()
    #np.savetxt(r"C:\Users\16578\Desktop\FRPCASMMA投稿\S.txt",S,fmt="%e")
    return L
#
# if __name__ == "__main__":
#     starttime = time.time()
#     main()
#     endtime = time.time()
#     print('总共的时间为:', round(endtime - starttime, 4), '秒')
