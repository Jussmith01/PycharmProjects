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
# ------------
# AM1 vs Act
# ------------
user = os.environ['USER']

dir = '/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/train_05/'
dir = '/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/train_06/'
file = 'RMSEperATOM.dat'

data1 = gt.getfltsfromfile('/home/' + user + dir + file, [0])
data2 = gt.getfltsfromfile('/home/' + user + dir + file, [1])
data3 = gt.getfltsfromfile('/home/' + user + dir2 + file, [1])

print('Datasize: ' + str(data1.shape[0]))

font = {'family' : 'Bitstream Vera Sans',
            'weight' : 'normal',
            'size'   : 14}

plt.rc('font', **font)

plt.scatter(data1, data2, color='blue', label='ANN - GDB3',linewidth=2)
plt.plot(data1, data2, color='blue', label='ANN - GDB3',linewidth=2)

plt.title("RMSE per Atom for each dataset")
plt.xlabel('File index')
plt.ylabel('RMSE (eV/atom)')
plt.legend(bbox_to_anchor=(0.025, 0.975), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()
