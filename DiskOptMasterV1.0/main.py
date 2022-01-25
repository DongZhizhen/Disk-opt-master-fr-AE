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
from scipy import interpolate
import pylab as pl

mn = 4  # 网格的中间节点数（输入！！！）
bn = 1  # 叶片网格层数
Size = 50  # 粒子组数

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
facenodes = np.array(pd.read_excel(r"nodedata.xlsx", sheet_name='dataface', header=None))  # 导入左边界节点、坐标
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

sita_2 = np.sort(th_unique, axis=1)  # 正序，非际序
sita_3 = abs(np.sort(-th_unique, axis=1))  # 倒序，实际序

# 极坐标按角坐标分类，按实际序
sita_part = np.zeros((int((nodetrans.shape[0] - facenodes.shape[0] * sita_2.shape[1]) / (bn + 1) + facenodes.shape[0]),
                      facenodes.shape[1], sita_2.shape[1]))  # 注意叶片网格层数
for i in range(sita_2.shape[1]):
    tt = 0
    for j in range(nodetrans.shape[0]):
        if (nodetrans[j, 3] <= sita_3[0, i] + 0.001).all() and (nodetrans[j, 3] >= sita_3[0, i] - 0.001).all():
            for k in range(facenodes.shape[1]):
                sita_part[tt, k, i] = nodetrans[j, k]
            tt += 1

''' 2D端面内节点关联 '''
# 2D端面内节点期望坐标计算
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
for j in range(mn + 2):  # 数据整理
    for i in range(copnds.shape[0]):
        coupled_2D_index[tt, 0] = coupled_2D_temp[i, j]
        tt += 1

coupled_2D_ntr = np.zeros((coupled_2D_index.shape[0], 1))  # 后续使用
coupled_2D_ntz = np.zeros((coupled_2D_index.shape[0], Size))  # 后续使用
for i in range(coupled_2D_index.shape[0]):
    coupled_2D_ntr[i, 0] = nodetrans[int(coupled_2D_index[i, 0]) - 1, 0]

# 3D待优化区内节点索引匹配
trshd = 0.2  # 比较关键
coupled_3D_index = np.zeros((coupled_2D_index.shape[0], mn + 2))
coupled_3D_index[:, 5] = coupled_2D_index[:, 0]
for i in range(coupled_2D_index.shape[0]):
    for j in range(zn):
        sita_temp = sita_part[:, :, j]
        temp2_r_index = sita_temp[
            np.intersect1d(np.where(sita_temp[:, 0] < nodetrans[int(coupled_2D_index[i]), 0] + trshd),
                           np.where(nodetrans[int(coupled_2D_index[i]), 0] - trshd < sita_temp[:, 0])), 0]
        temp2_z_index = sita_temp[
            np.intersect1d(np.where(sita_temp[:, 1] < nodetrans[int(coupled_2D_index[i]), 1] + trshd * 2),
                           np.where(nodetrans[int(coupled_2D_index[i]), 1] - trshd * 2 < sita_temp[:, 1])), 1]

        for k in temp2_r_index:  # 提高通用性
            if np.intersect1d(k, temp2_z_index) is not None:
                coupled_3D_index[i, j] = sita_temp[int(k), 2]

# 3D补充节(过渡区)节点索引匹配
trshd = 0.001  # 比较关键
suply = np.array(pd.read_excel(r"nodedata.xlsx", sheet_name='suply', header=None))[:32,
        :]  # 从mapping文件的suply_det变量获取，输入excel
coupled_2D_index = np.vstack([coupled_2D_index, suply])  # 2D补充

suply_r = np.zeros((suply.shape[0], 1))
suply_z = np.zeros((suply.shape[0], 1))
for i in range(suply.shape[0]):
    suply_r[i] = nodetrans[int(suply[i]) - 1, 0]
    suply_z[i] = nodetrans[int(suply[i]) - 1, 1]

coupled_3D_temp = np.zeros((1, mn + 2))
for i in range(suply.shape[0]):
    tt = 0
    temp3_r_index = np.intersect1d(np.where(nodetrans[:, 0] < suply_r[i] + trshd),
                                   np.where(suply_r[i] - trshd < nodetrans[:, 0]))
    temp3_z_index = np.intersect1d(np.where(nodetrans[:, 1] < suply_z[i] + trshd),
                                   np.where(suply_z[i] - trshd < nodetrans[:, 1]))

    for j in temp3_r_index:  # 通用性提高
        if np.intersect1d(j, temp3_z_index) is not None:  # & & find(temp_rz_index(j) == temp_z_index(:))
            coupled_3D_temp[0, tt] = j
            tt += 1

    coupled_3D_index = np.concatenate((coupled_3D_index, coupled_3D_temp), axis=0)  # 拼接行

''' 筛选除待更新节点外的固定节点 '''
delete_3D_index = coupled_3D_index.reshape(1, coupled_3D_index.shape[0] * coupled_3D_index.shape[1]).astype(int)
nodefixed = np.delete(nodeinitial, delete_3D_index - 1, 0)  # 删除待更新节点坐标信息

''' cubic插值（求解chebyshev点对应的Z坐标值） '''
left_R = np.zeros((left.shape[0], 1))
left_Z = np.zeros((left.shape[0], 1))
right_R = np.zeros((right.shape[0], 1))
right_Z = np.zeros((right.shape[0], 1))
for i in range(left.shape[0])
    left_R[i] = nodetrans[left[i] - 1, 0]  # left(:,3)
    left_Z[i] = nodetrans[left[i] - 1, 1]  # left(:,2)
    
# for kind in ["cubic"]:  # 插值方式["nearest", "zero", "slinear", "quadratic", "cubic"] , "nearest","zero"为阶梯插值, slinear 线性插值, "quadratic","cubic" 为2阶、3阶B样条曲线插值
    f = interpolate.interp1d(left_R, left_Z, kind="cubic")
    left_cheby_z = f(left_cheby_r)        # R的值不能重复   作为PSO初始边界

for i in range(right.shape[0])
    right_R[i] = nodetrans[right[i] - 1, 0]
    right_Z[i] = nodetrans[right[i] - 1, 1]

f = interpolate.interp1d(right_R, right_Z, kind="cubic")
right_cheby_z = f(right_cheby_r)      # R的值不能重复   作为PSO初始边界

''' PSOA初始化 '''


