"""

    main.m

    Description:
        This program is used to optimize the disk in aero engine.


    All the codes had been written by Zhizhen Dong in 2022

"""

import numpy as np
import math
import pandas as pd
from functools import reduce
from numpy import exp, abs, angle
import numpy as py
import xlwt as xlwt

mn = 4  # 网格的中间节点数（输入！！！）
Size = 50        # 粒子组数

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
def chebyshev(a_t, b_t, n_t):  # k坐标变换，多项式阶数
    k1 = (a_t + b_t) / 2
    k2 = (-a_t + b_t) / 2
    c = np.zeros(n_t)

    for i in range(n_t):
        c_temp = (k1 + k2 * math.cos((2 * i + 1) * math.pi / (2 * (n_t + 1))))
        c[i] = float('%.2f' % c_temp)

    return c


a = 32.5  # 范围
b1 = 95
b2 = 98.8
n = 8  # 阶数-1

left_cheby_r = chebyshev(a, b1, n)
right_cheby_r = chebyshev(a, b2, n)

''' 导入全部节点信息（编号、坐标）极坐标换算 '''
nodeini = pd.read_excel(r"nodedata.xlsx", sheet_name='datapaper', header=None)
nodeinitial = np.array(nodeini)
nn = np.array([nodeinitial[:, 0]]).astype(int)
xini = np.array([nodeinitial[:, 1]])
yini = np.array([nodeinitial[:, 2]])
zini = np.array([nodeinitial[:, 3]])


def cart2pol(x_t, y_t, z_t=None):
    rho_t = np.sqrt(x_t ** 2 + y_t ** 2)
    phi_t = np.arctan2(y_t, x_t)
    if z_t is None:
        return (phi_t, rho_t)
    else:
        return (phi_t, rho_t, z_t)


def pol2cart(rho_t, phi_t, z_t=None):
    x_t = rho_t * np.cos(phi_t)
    y_t = rho_t * np.sin(phi_t)
    if z_t is None:
        return (x_t, y_t)
    else:
        return (x_t, y_t, z_t)


th, r, z = cart2pol(xini, zini, yini)
nodetrans = np.hstack((r.T, z.T, nn.T, th.T))  # 节点、单元全部信息整合

# left、right网格索引同节点编号


''' 周向各层角坐标提取 '''
th_t = np.around(th, 4)


def del_repeatnum(s):
    s1 = []
    for i in s:
        if i not in s1:
            s1.append(i)
    return s1


th_unique = np.array(del_repeatnum(th_t.T)).T  # 周向各层角坐标提取
zn = th_unique.shape[1] - 1  # 周向层数
delta_th = (np.max(th_unique) - np.min(th_unique)) / zn  # 周向夹角

''' 2D端面内节点关联 '''
# 2D端面内节点期望坐标计算
facenodes = np.array(pd.read_excel(r"nodedata.xlsx", sheet_name='dataface', header=None))  # 导入左边界节点、坐标
copnds1 = np.array(pd.read_excel(r"nodedata.xlsx", sheet_name='core', header=None))[:, :2]
copnds2 = np.array(pd.read_excel(r"nodedata.xlsx", sheet_name='expanded', header=None))[:, :2]
copnds = np.vstack((copnds1, copnds2))

deltadis_r = np.zeros((copnds.shape[0], 1))
deltadis_z = np.zeros((copnds.shape[0], 1))

for i in range(copnds.shape[0]):
    discopnds_r = nodetrans[copnds[i, 1] - 1, 0] - nodetrans[copnds[i, 0] - 1, 0]
    discopnds_z = nodetrans[copnds[i, 1] - 1, 1] - nodetrans[copnds[i, 0] - 1, 1]
    deltadis_r[i] = discopnds_r / (mn + 1)
    deltadis_z[i] = discopnds_z / (mn + 1)

coupled_2D_temp = np.zeros((copnds.shape[0], mn + 2))  # 初始化全部端面耦合节点在APDL网格中的索引
coupled_2D_ndr = np.zeros((copnds.shape[0], mn + 2))  # 初始化全部端面耦合节点在APDL网格中的r坐标
coupled_2D_ndz = np.zeros((copnds.shape[0], mn + 2))  # 初始化全部端面耦合节点在APDL网格中的z坐标

for i in range(copnds.shape[0]):
    coupled_2D_ndr[i, 0] = nodetrans[copnds[i, 0] - 1, 0]
    coupled_2D_ndr[i, mn + 1] = nodetrans[copnds[i, 1] - 1, 0]
    coupled_2D_ndz[i, 0] = nodetrans[copnds[i, 0] - 1, 1]
    coupled_2D_ndz[i, mn + 1] = nodetrans[copnds[i, 1] - 1, 1]
    for j in range(1, mn + 1):
        coupled_2D_ndr[i, j] = nodetrans[copnds[i, 0] - 1, 0] + deltadis_r[i] * j
        coupled_2D_ndz[i, j] = nodetrans[copnds[i, 0] - 1, 1] + deltadis_z[i] * j

# 2D端面内节点索引匹配
deltadis_t = np.sqrt(deltadis_r ** 2 + deltadis_z ** 2)
trshd = np.min(deltadis_t) / 2  # 比较关键

temp_rz_face = np.zeros((facenodes.shape[0], 3))  # 端面nrz坐标索引
for i in range(facenodes.shape[0]):
    temp_rz_face[i, 0] = facenodes[i, 0]
    temp_rz_face[i, 1] = nodetrans[int(facenodes[i, 0] - 1), 0]
    temp_rz_face[i, 2] = nodetrans[int(facenodes[i, 0] - 1), 1]

for i in range(copnds.shape[0]):
    for j in range(mn + 2):
        temp_r_index = temp_rz_face[reduce(np.intersect1d, [np.where(temp_rz_face[:, 1] < coupled_2D_ndr[i, j] + trshd),
                                                            np.where(temp_rz_face[:, 1] > coupled_2D_ndr[
                                                                i, j] - trshd)]), 0]  # 多数组求交
        temp_z_index = temp_rz_face[np.intersect1d(np.where(temp_rz_face[:, 2] < coupled_2D_ndz[i, j] + trshd),
                                                   np.where(temp_rz_face[:, 2] > coupled_2D_ndz[i, j] - trshd)), 0]

        for k in temp_r_index:  # 提高通用性
            if np.intersect1d(k, temp_z_index) is not None:
                coupled_2D_temp[i, j] = nodetrans[int(k), 2]

tt = 0
coupled_2D_index = np.zeros((copnds.shape[0] * (mn + 2), 1))
for j in range(mn + 2):           # 数据整理
    for i in range(copnds.shape[0]):
        coupled_2D_index[tt, 0] = coupled_2D_temp[i, j]
        tt += 1

coupled_2D_ntr = np.zeros((coupled_2D_index.shape[0], 1))             # 后续使用
coupled_2D_ntz = np.zeros((coupled_2D_index.shape[0], Size))           # 后续使用
for i in range(coupled_2D_index.shape[0]):
    coupled_2D_ntr[i, 0] = nodetrans[int(coupled_2D_index[i, 0]) - 1, 0]

# 3D待优化区内节点索引匹配
trshd = 0.2    # 比较关键
coupled_3D_index = np.zeros((coupled_2D_index.shape[0], mn+2))
coupled_3D_index[:, 5] = coupled_2D_index
for i in range(coupled_2D_index.shape[0]):
    for j in range(zn):
