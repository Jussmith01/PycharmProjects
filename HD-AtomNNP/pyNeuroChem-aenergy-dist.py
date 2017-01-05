__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import numpy as np
import graphtools as gt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import matplotlib as mpl

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/smallAEV_testing/train_384-256-128-64-1_c08e/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

#dtdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_01/data/'
#xyz,typ,Eact,tmp    = gt.readncdat(dtdir + 'gdb11_s01-1_test.dat',np.float32)
file = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_02/testdata/gdb11_s02-0_test.dat'
xyz,typ,Eact,tmp = gt.readncdat(file,np.float32)

file2 = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_02/testdata/gdb11_s02-3_test.dat'
xyz2,typ2,Eact2,tmp2 = gt.readncdat(file2,np.float32)

file3 = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_02/testdata/gdb11_s02-7_test.dat'
xyz3,typ3,Eact3,tmp3 = gt.readncdat(file3,np.float32)

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile, saefile, nnfdir, 0)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ)

# Compute Energies of Conformations
E1 = nc.energy()

# Atomic energy return
AE1 = np.copy(nc.aenergies(sae=False))

print (AE1.shape)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz2,types=typ2)

# Compute Energies of Conformations
E2 = nc.energy()

# Atomic energy return
AE2 = np.copy(nc.aenergies(sae=False))

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz3,types=typ3)

# Compute Energies of Conformations
E3 = nc.energy()

# Atomic energy return
AE3 = np.copy(nc.aenergies(sae=False))

print (typ)
print (typ2)
print (typ3)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 9}

plt.rc('font', **font)

fig, axes = plt.subplots(nrows=1, ncols=1)

axes.set_title("Comparison of atomic energy from ANI-c08e(384) for hydrogen 1 in Ethane and Ethene conformers")
axes.set_ylabel('Energy count')
axes.set_xlabel('Energy (Ha)')

axes.hist(AE1[:,2], 50, color='red',normed=True, label='Ethane hydrogen 1',linewidth=2,alpha=0.6)
axes.hist(AE2[:,2], 50, color='blue',normed=True, label='Ethene hydrogen 1',linewidth=2,alpha=0.6)
axes.hist(AE3[:,2], 50, color='orange',normed=True, label='Ethyne hydrogen 1',linewidth=2,alpha=0.6)

#f = open(file + '.aedata', 'w')

#for i in range(0,len(typ)):
#    f.write( typ[i] + ' ' + "{:.7f}".format(AE[i]) + '\n' )

#f.close()

plt.legend(bbox_to_anchor=(0.6, 0.98), loc=2, borderaxespad=0., fontsize=14)

# -----
# PLOT
# -----
plt.show()