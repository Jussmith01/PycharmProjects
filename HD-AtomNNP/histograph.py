__author__ = 'jujuman'

import numpy as np
import statsmodels.api as sm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
import random
import graphtools as gt

def histograph(dir, file, bins, label, color='black'):
    user = os.environ['USER']
    data = gt.getfltsfromfile('/home/' + user + dir + file,' ', [0])

    plt.hist(data, bins, color=color,normed=0, label=label,linewidth=2)


# -----------------------
cmap = mpl.cm.brg
# ------------5412.mordor
# AM1 vs Act
# ------------
user = os.environ['USER']

font = {'family' : 'Bitstream Vera Sans',
            'weight' : 'normal',
            'size'   : 14}

plt.rc('font', **font)

histograph('/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/dnntsgdb11_06/','output.dat',50,'','blue')
histograph('/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/dnntsgdb11_05/','output.dat',50,'','red')
histograph('/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/dnntsgdb11_04/','output.dat',50,'','green')

plt.title("Energy Differences Between 50 Random Structures\nPolypeptide Chain: H-Gly-Pro-Hyp-Gly-Ala-Gly-OH")
#plt.xlabel('Conformation Pair (Count 49)')
plt.xlabel('Theoretical dE (Hartree)')
plt.ylabel('Calculated dE (Hartree)')
plt.legend(bbox_to_anchor=(0.025, 0.975), loc=2, borderaxespad=0.)


# -----
# PLOT
# -----
plt.show()