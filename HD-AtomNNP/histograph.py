__author__ = 'jujuman'

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
import random
import graphtools as gt

def histograph(ax, dir, file, bins, norm, label, color='black',alpha=1.0):
    user = os.environ['USER']
    data = gt.getfltsfromfile('/home/' + user + dir + file,' ', [0])

    ax.set_title("Data: " + file)
    ax.set_ylabel('Normalized distant count')
    ax.set_xlabel('Distance ($\AA$)')

    ax.hist(data, bins, color=color,normed=norm, label=label,linewidth=2,alpha=alpha)

# -----------------------
cmap = mpl.cm.brg
# ------------5412.mordor
# AM1 vs Act
# ------------
user = os.environ['USER']

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 9}

plt.rc('font', **font)

fig, axes = plt.subplots(nrows=2, ncols=5)

N = 50

GDBstr = "08"

print ('Plotting All')
histograph(axes.flat[0],'/Research/distance_data/','distance_dist_HH_GDB' + GDBstr + '.dat',N,1,'ANI-1 Training Set','blue',1.0)
histograph(axes.flat[1],'/Research/distance_data/','distance_dist_HO_GDB' + GDBstr + '.dat',N,1,'ANI-1 Training Set','blue',1.0)
histograph(axes.flat[2],'/Research/distance_data/','distance_dist_HC_GDB' + GDBstr + '.dat',N,1,'ANI-1 Training Set','blue',1.0)
histograph(axes.flat[3],'/Research/distance_data/','distance_dist_HN_GDB' + GDBstr + '.dat',N,1,'ANI-1 Training Set','blue',1.0)
histograph(axes.flat[4],'/Research/distance_data/','distance_dist_OO_GDB' + GDBstr + '.dat',N,1,'ANI-1 Training Set','blue',1.0)
histograph(axes.flat[5],'/Research/distance_data/','distance_dist_OC_GDB' + GDBstr + '.dat',N,1,'ANI-1 Training Set','blue',1.0)
histograph(axes.flat[6],'/Research/distance_data/','distance_dist_ON_GDB' + GDBstr + '.dat',N,1,'ANI-1 Training Set','blue',1.0)
histograph(axes.flat[7],'/Research/distance_data/','distance_dist_CC_GDB' + GDBstr + '.dat',N,1,'ANI-1 Training Set','blue',1.0)
histograph(axes.flat[8],'/Research/distance_data/','distance_dist_CN_GDB' + GDBstr + '.dat',N,1,'ANI-1 Training Set','blue',1.0)
histograph(axes.flat[9],'/Research/distance_data/','distance_dist_NN_GDB' + GDBstr + '.dat',N,1,'ANI-1 Training Set','blue',1.0)

# -----
# PLOT
# -----
plt.show()