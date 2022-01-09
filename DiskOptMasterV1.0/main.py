"""

    main.m

    Description:
        This program is used to optimize the disk in aero engine.


    All the codes had been written by Zhizhen Dong in 2022

"""

import numpy as np
import pandas as pd
import numpy as py
import xlwt as xlwt

mn = 4  # 网格的中间节点数（输入！！！）

# 导入数据
# couplednodes = pd.read_table(r"couplednodes.txt", sep=",", header=None) # 按行读文本
couplednodes = np.loadtxt("couplednodes.txt", dtype=int, delimiter=' ')
left_1 = np.array([[2956, 2947, 2948]])
left_2 = np.array([couplednodes[:, 0]])
left_3 = np.array([[2611, 2612, 2623, 2624]])
left = np.vstack((left_1.T, left_2.T, left_3.T))

right_1 = np.array([[2730, 2734, 2739]])
right_2 = np.array([couplednodes[:, mn + 1]])
right_3 = np.array([[2669, 2672, 2681, 2684]])
right = np.vstack((right_1.T, right_2.T, right_3.T))

# 计算chebyshev点
a = 32.5  # 范围
b1 = 95
b2 = 98.8
n = 8  # 阶数-1

