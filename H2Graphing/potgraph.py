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
# ------------
# AM1 vs Act
# ------------
user = os.environ['USER']
#dir = '/Research/ANN-Test-Data/GDB-11/train5/'
#dir = '/Research/ANN-Test-Data/GDB-11/train3/'
'''
N = 8
for i in range(0,N):
    file = 'gdb11_s02-' + str(i) + '_train.dat_graph.dat'

    data1 = getfltsfromfile('/home/' + user + dir + file, [0])
    data2 = getfltsfromfile('/home/' + user + dir + file, [1])
    data3 = getfltsfromfile('/home/' + user + dir + file, [2])

    data2 = (data3 - data2)*(data3 - data2)

    print('Datasize: ' + str(data2.shape[0]))

    font = {'family' : 'Bitstream Vera Sans',
            'weight' : 'normal',
            'size'   : 14}

    plt.rc('font', **font)

    s = 288

    plt.plot(data1, data2, color=cmap((i+1)/float(N)), label=str(i),linewidth=1)
    plt.scatter(data1, data2, color=cmap((i+1)/float(N)), label=str(i),linewidth=2)
    #plt.scatter(data2, data3, color=cmap((i+1)/float(N)), label=str(i),linewidth=1)
'''

dir = '/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/train_05/'
dir2 = '/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/train_06/'
#file = 'RMSEperATOM.dat'
file = 'gdb11_s05-99_test.dat_graph'

data1 = gt.getfltsfromfile('/home/' + user + dir + file, [0])
data2 = gt.getfltsfromfile('/home/' + user + dir + file, [1])
data3 = gt.getfltsfromfile('/home/' + user + dir + file, [2])
data4 = gt.getfltsfromfile('/home/' + user + dir2 + file, [2])

data2 = gt.calculateelementdiff(data2)
data3 = gt.calculateelementdiff(data3)
data4 = gt.calculateelementdiff(data4)

rmse5 = gt.calculaterootmeansqrerror(data2[:,1],data3[:,1]) / 63.0
rmse6 = gt.calculaterootmeansqrerror(data2[:,1],data4[:,1]) / 63.0

print('Datasize: ' + str(data1.shape[0]))

font = {'family' : 'Bitstream Vera Sans',
            'weight' : 'normal',
            'size'   : 14}

plt.rc('font', **font)

#print(data2)

plt.plot(data2[:,0], data2[:,1], color='black', label='wB97X/6-31G*',linewidth=2)
plt.scatter(data2[:,0], data2[:,1], color='black',linewidth=4)
plt.plot(data3[:,0], data3[:,1],'r--', color='blue', label='ANN - GDB-5 RMSE: ' + "{:.6f}".format(rmse5) + "eV/atom",linewidth=2)
plt.scatter(data3[:,0], data3[:,1], color='blue',linewidth=4)
plt.plot(data4[:,0], data4[:,1],'r--', color='red', label='ANN - GDB-6 RMSE: ' + "{:.6f}".format(rmse6) + "eV/atom",linewidth=2)
plt.scatter(data4[:,0], data4[:,1], color='red',linewidth=4)
#plt.scatter(data1, data2, color='blue', label='ANN - GDB3',linewidth=2)
#plt.plot(data1, data2, color='blue', label='ANN - GDB3',linewidth=2)
#plt.scatter(data1, data3, color='green', label='ANN - GDB3',linewidth=2)


#plt.title("H-Gly-Pro-Hyp-Gly-Ala-Gly-OH")
plt.title("Energy Differences Between 50 Random Structures\nPolypeptide Chain: H-Gly-Pro-Hyp-Gly-Ala-Gly-OH")
plt.xlabel('Conformation Pair (Count 49)')
plt.ylabel('Calculated dE (Hartree)')
plt.legend(bbox_to_anchor=(0.025, 0.975), loc=2, borderaxespad=0.)


# -----
# PLOT
# -----
plt.show()
