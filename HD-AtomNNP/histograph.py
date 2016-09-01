__author__ = 'jujuman'

import numpy as np
import statsmodels.api as sm
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
        'size'   : 14}

plt.rc('font', **font)

print ('Plotting 05')
histograph('/Research/GDB-11-wB97X-6-31gd/','dist_all.chk',150,1,'All upto GDB-7','red')
print ('Plotting 07 All')
histograph('/Research/GDB-11-wB97X-6-31gd/dnnts_testdata/PeptideCases/testdata/','pp_01_dist.chk',150,1,'pp01','orange',0.75)
print ('Plotting 06')
histograph('/Research/GDB-11-wB97X-6-31gd/dnnts_testdata/PeptideCases/testdata/','pp_02_dist.chk',150,1,'pp02','blue',0.5)
#print ('Plotting 05')
#histograph('/Research/GDB-11-wB98X-6-31gd/dnntsgdb11_05/','distchk.dat',125,0,'GDB-5','red')
#print ('Plotting 04')
#histograph('/Research/GDB-11-wB98X-6-31gd/dnntsgdb11_04/','distchk.dat',100,0,'GDB-4','green')

plt.title("Atomic Distance Distribution (OH only)")
plt.ylabel('Distance Count')
plt.xlabel('Distance ($\AA$)')
plt.legend(bbox_to_anchor=(0.025, 0.975), loc=2, borderaxespad=0.)


# -----
# PLOT
# -----
plt.show()