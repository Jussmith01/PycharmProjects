__author__ = 'jujuman'

import numpy as np
import statsmodels.api as sm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
import random
import graphtools as gt


# -----------------------
cmap = mpl.cm.brg
# ------------5412.mordor
# AM1 vs Act
# ------------
user = os.environ['USER']

dir = '/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/train_06AA/'

file = 'pp_01_test.dat_graph'

data1 = gt.getfltsfromfile('/home/' + user + dir + file, [0])
data2 = gt.getfltsfromfile('/home/' + user + dir + file, [1])
data3 = gt.getfltsfromfile('/home/' + user + dir + file, [2])


data2 = gt.calculateelementdiff(data2)
data3 = gt.calculateelementdiff(data3)

rmse = gt.calculaterootmeansqrerror(data2[:,1],data3[:,1]) / 64.0


print('Datasize: ' + str(data1.shape[0]))

font = {'family' : 'Bitstream Vera Sans',
            'weight' : 'normal',
            'size'   : 14}

plt.rc('font', **font)


plt.plot(data2[:,0], data2[:,1], color='black', label='wB97X/6-31G*',linewidth=2)
plt.scatter(data2[:,0], data2[:,1], color='black',linewidth=4)
plt.plot(data3[:,0], data3[:,1],'r--', color='red', label='ANN - GDB-5 RMSE: ' + "{:.6f}".format(rmse) + "Ha/atom",linewidth=2)
plt.scatter(data3[:,0], data3[:,1], color='red',linewidth=4)

plt.title("Energy Differences Between 50 Random Structures\nPolypeptide Chain: H-Gly-Pro-Hyp-Gly-Ala-Gly-OH")
#plt.xlabel('Conformation Pair (Count 49)')
plt.xlabel('Theoretical dE (Hartree)')
plt.ylabel('Calculated dE (Hartree)')
plt.legend(bbox_to_anchor=(0.025, 0.975), loc=2, borderaxespad=0.)


# -----
# PLOT
# -----
plt.show()
