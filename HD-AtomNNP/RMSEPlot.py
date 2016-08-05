__author__ = 'jujuman'

import numpy as np
import statsmodels.api as sm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
import random
import graphtools as gt

def graphRMSEperATOM(dir, label, color='black'):
    user = os.environ['USER']
    file = 'RMSEperATOM.dat'
    data1 = gt.getfltsfromfile('/home/' + user + dir + file, [0])
    data2 = gt.getfltsfromfile('/home/' + user + dir + file, [1])

    plt.scatter(data1, data2, color=color, label=label,linewidth=2)
    plt.plot(data1, data2, color=color,linewidth=2)



# -----------------------
cmap = mpl.cm.brg
# ------------
# AM1 vs Act
# ------------

font = {'family' : 'Bitstream Vera Sans',
            'weight' : 'normal',
            'size'   : 14}

plt.rc('font', **font)

graphRMSEperATOM('/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/train_06AA/','ANN - GDB6AA','blue')
graphRMSEperATOM('/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/train_06_2/','ANN - GDB6','red')

plt.title("RMSE per Atom for each dataset")
plt.xlabel('File index')
plt.ylabel('RMSE (eV/atom)')
plt.legend(bbox_to_anchor=(0.025, 0.975), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()
