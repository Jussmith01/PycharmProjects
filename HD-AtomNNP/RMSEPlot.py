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
    data1 = gt.getfltsfromfile('/home/' + user + dir + file, ' ', [0])
    data2 = 23.0609 * gt.getfltsfromfile('/home/' + user + dir + file, ' ', [1])

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

graphRMSEperATOM('/Research/trainingcases/wB97X-631gd-train-highgarden/nw-64-64-64-32/train_07-a3.1A_r4.5/','ANN - c07b (Rca: 3.1$\AA$;Rcr: 4.5$\AA$)','blue')
#graphRMSEperATOM('/Research/trainingdata/wB97X-631gd-train-mordor/train_07-3.4A/','ANN - GDB7 3.4A cutoff','red')
#graphRMSEperATOM('/Research/GDB-11-wB98X-6-31gd/train_07/tests1/','ANN - GDB6','red')

plt.title("Error per atom per data file of c07b data set\nNetwork tag: c07b-64-64-64-32-r32-a-8x8-ac3.1-rc4.5")
plt.xlabel('File index')
plt.ylabel('RMSE (kcal/mol/atom)')
plt.legend(bbox_to_anchor=(0.025, 0.975), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()
