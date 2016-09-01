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

dir1 = '/Research/trainingcases/wB97X-631gd-train-comet/train_07_a2.9A_r5.2A/'

file = 'pp_02_test.dat_graph'
#file = 'polypep_test.dat_graph'
#file = 'aminoacid_00-12_test.dat_graph'
#file = 'benzamide_conformers-0_test.dat_graph'
#file = 'pentadecane_test.dat_graph'
#file = 'retinolconformer_test.dat_graph'

#data1 = gt.getfltsfromfile('/home/' + user + dir1 + file, ' ', [0])
data0 = gt.getfltsfromfile('/home/' + user + dir1 + file, ' ', [1])
data1 = gt.getfltsfromfile('/home/' + user + dir1 + file, ' ', [2])

data0 = gt.calculateelementdiff(data0)
data1 = gt.calculateelementdiff(data1)

rmse1 = 27.2113825435 * gt.calculaterootmeansqrerror(data0[:,1],data1[:,1]) / 24.0

print('Datasize: ' + str(data1.shape[0]))

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 8}

plt.rc('font', **font)

data = data0
plt.plot(data[:,0], data[:,1], color='black', label='wB97X/6-31G*',linewidth=2)
plt.scatter(data[:,0], data[:,1], color='black',linewidth=4)

data = data1
plt.plot(data[:,0], data[:,1],'r--', color='blue', label='ANN - GDB-7 (6.1$\AA$) RMSE: ' + "{:.6f}".format(rmse1) + "eV/atom",linewidth=2)
plt.scatter(data[:,0], data[:,1], color='blue',linewidth=4)

plt.title("Energy Differences Between 8 Random Conformations of a Peptide")
#plt.xlabel('Conformation Pair (Count 49)')
plt.xlabel('Theoretical dE (Hartree)')
plt.ylabel('Calculated dE (Hartree)')
plt.legend(bbox_to_anchor=(0.3, 0.975), loc=2, borderaxespad=0.)


# -----
# PLOT
# -----
plt.show()
