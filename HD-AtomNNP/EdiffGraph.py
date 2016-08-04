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

dir1 = '/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/train_04/'
dir2 = '/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/train_05/'
dir3 = '/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/train_06_2/'
dir4 = '/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/train_06AA/'

file = 'polypep_test.dat_graph'

data1 = gt.getfltsfromfile('/home/' + user + dir1 + file, [0])
data2 = gt.getfltsfromfile('/home/' + user + dir1 + file, [1])
data3 = gt.getfltsfromfile('/home/' + user + dir1 + file, [2])
data4 = gt.getfltsfromfile('/home/' + user + dir2 + file, [2])
data5 = gt.getfltsfromfile('/home/' + user + dir3 + file, [2])
data6 = gt.getfltsfromfile('/home/' + user + dir4 + file, [2])

data2 = gt.calculateelementdiff(data2)
data3 = gt.calculateelementdiff(data3)
data4 = gt.calculateelementdiff(data4)
data5 = gt.calculateelementdiff(data5)
data6 = gt.calculateelementdiff(data6)

rmse4 = gt.calculaterootmeansqrerror(data2[:,1],data3[:,1]) / 63.0
rmse5 = gt.calculaterootmeansqrerror(data2[:,1],data4[:,1]) / 63.0
rmse6 = gt.calculaterootmeansqrerror(data2[:,1],data5[:,1]) / 63.0
rmse6AA = gt.calculaterootmeansqrerror(data2[:,1],data6[:,1]) / 63.0

print('Datasize: ' + str(data1.shape[0]))

font = {'family' : 'Bitstream Vera Sans',
            'weight' : 'normal',
            'size'   : 14}

plt.rc('font', **font)


plt.plot(data2[:,0], data2[:,1], color='black', label='wB97X/6-31G*',linewidth=2)
plt.scatter(data2[:,0], data2[:,1], color='black',linewidth=4)
plt.plot(data3[:,0], data3[:,1],'r--', color='blue', label='ANN - GDB-4 RMSE: ' + "{:.6f}".format(rmse4) + "Ha/atom",linewidth=2)
plt.scatter(data3[:,0], data3[:,1], color='blue',linewidth=4)#
plt.plot(data4[:,0], data4[:,1],'r--', color='red', label='ANN - GDB-5 RMSE: ' + "{:.6f}".format(rmse5) + "Ha/atom",linewidth=2)
plt.scatter(data4[:,0], data4[:,1], color='red',linewidth=4)
plt.plot(data5[:,0], data5[:,1],'r--', color='green', label='ANN - GDB-6 RMSE: ' + "{:.6f}".format(rmse6) + "Ha/atom",linewidth=2)
plt.scatter(data5[:,0], data5[:,1], color='green',linewidth=4)
plt.plot(data6[:,0], data6[:,1],'r--', color='orange', label='ANN - GDB-6AA RMSE: ' + "{:.6f}".format(rmse6AA) + "Ha/atom",linewidth=2)
plt.scatter(data6[:,0], data6[:,1], color='orange',linewidth=4)
'''

plt.plot(data2[:,1], data2[:,1], color='black', label='wB97X/6-31G*',linewidth=2)
plt.scatter(data2[:,1], data3[:,1], color='blue', label='ANN - GDB-4 RMSE: ' + "{:.6f}".format(rmse4) + "Ha/atom",linewidth=2)
plt.scatter(data2[:,1], data4[:,1], color='red', label='ANN - GDB-5 RMSE: ' + "{:.6f}".format(rmse5) + "Ha/atom",linewidth=2)
plt.scatter(data2[:,1], data5[:,1], color='green', label='ANN - GDB-6 RMSE: ' + "{:.6f}".format(rmse6) + "Ha/atom",linewidth=2)
plt.scatter(data2[:,1], data6[:,1], color='orange', label='ANN - GDB-6AA RMSE: ' + "{:.6f}".format(rmse6AA) + "Ha/atom",linewidth=2)
'''

plt.title("Energy Differences Between 50 Random Structures\nPolypeptide Chain: H-Gly-Pro-Hyp-Gly-Ala-Gly-OH")
#plt.xlabel('Conformation Pair (Count 49)')
plt.xlabel('Theoretical dE (Hartree)')
plt.ylabel('Calculated dE (Hartree)')
plt.legend(bbox_to_anchor=(0.025, 0.975), loc=2, borderaxespad=0.)


# -----
# PLOT
# -----
plt.show()
