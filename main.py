#矩阵处理库
import numpy as np
#自建库
import GAMA
import GIP
import BNNR
import 实验的main集合
#机器学习库
from sklearn import metrics
from sklearn.model_selection import KFold
from sklearn.metrics import roc_curve, auc
import scipy
import scipy.io as sio
#基础库
import copy
import math
import ipdb
import time
#画图的库
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import ConnectionPatch

def main():
    # w = 0.2
    # GIP.CaT_GIP(w)
    starttime = time.time()
    for i in range(1):
        fpr,tpr,AUC,roc_sum=实验的main集合.CV(5)
        # print(len(fpr[0]))
    # 实验的main集合.main_w()
    #实验的main集合.CV_miRNA()
    #实验的main集合.CV_SM()
    #实验的main集合.main_5fold_Globalfold()
    # 实验的main集合.cv_1()
    #实验的main集合.CV(5)
    # for i in range(100):
    #     实验的main集合.new()
    #实验的main集合.onetimes()
    endtime = time.time()
    print('该实验花费时间为{0}',endtime-starttime)

if __name__ == '__main__':
    main()
