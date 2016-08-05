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
dir2 = '/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/train_05/'
dir3 = '/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/train_06/'

file = 'pp_01_test.dat_graph'

data1 = gt.getfltsfromfile('/home/' + user + dir + file, [0])
data2 = gt.getfltsfromfile('/home/' + user + dir + file, [1])
data3 = gt.getfltsfromfile('/home/' + user + dir + file, [2])
data4 = gt.getfltsfromfile('/home/' + user + dir2 + file, [2])
data5 = gt.getfltsfromfile('/home/' + user + dir3 + file, [2])

rmse1 = gt.calculaterootmeansqrerror(data2,data3)
rmse2 = gt.calculaterootmeansqrerror(data2,data4)
rmse3 = gt.calculaterootmeansqrerror(data2,data5)

print('Datasize: ' + str(data1.shape[0]))

font = {'family' : 'Bitstream Vera Sans',
            'weight' : 'normal',
            'size'   : 14}

plt.rc('font', **font)


plt.plot(data1, data2, color='black', label='wB97X/6-31G*',linewidth=2)
plt.scatter(data1, data2, color='black',linewidth=4)
plt.scatter(data1, data3, color='red', label='ANN - GDB-6AA RMSE: ' + "{:.6f}".format(rmse1) + "Ha",linewidth=4)

plt.scatter(data1, data4, color='green', label='ANN - GDB-62 RMSE: ' + "{:.6f}".format(rmse2) + "Ha",linewidth=4)
plt.scatter(data1, data5, color='blue', label='ANN - GDB-61 RMSE: ' + "{:.6f}".format(rmse3) + "Ha",linewidth=4)

plt.title("Energy Differences Between 50 Random Structures\nPolypeptide Chain: H-Gly-Pro-Hyp-Gly-Ala-Gly-OH")
#plt.xlabel('Conformation Pair (Count 49)')
plt.xlabel('Theoretical dE (Hartree)')
plt.ylabel('Calculated dE (Hartree)')
plt.legend(bbox_to_anchor=(0.025, 0.975), loc=2, borderaxespad=0.)


# -----
# PLOT
# -----
plt.show()
