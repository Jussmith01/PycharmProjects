__author__ = 'jujuman'

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
import random
import graphtools as gt

def histograph(dir, file, bins, norm, label, color='black',alpha=1.0):
    user = os.environ['USER']
    data = gt.getfltsfromfile('/home/' + user + dir + file,' ', [0])

    plt.hist(data, bins, color=color,normed=norm, label=label,linewidth=2,alpha=alpha)


# -----------------------
cmap = mpl.cm.brg
# ------------5412.mordor
# AM1 vs Act
# ------------
user = os.environ['USER']

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 24}

plt.rc('font', **font)

print ('Plotting All')
histograph('/Research/GDB-11-wB97X-6-31gd/','dist_all.chk',150,1,'ANI-1 Training Set','blue',0.5)

#plt.title("Heavy atom  atomic distribution")
plt.ylabel('Normalized distant count')
plt.xlabel('Distance ($\AA$)')
plt.legend(bbox_to_anchor=(0.6, 0.975), loc=2, borderaxespad=0.)


# -----
# PLOT
# -----
plt.show()
