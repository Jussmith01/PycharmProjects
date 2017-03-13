import pygau09tools as g09
import hdnntools as hdn
import numpy as np
import matplotlib.pyplot as plt
import pyNeuroChem as pync

# Set required files for pyNeuroChem
anipath  = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk'
cnstfile = anipath + '/rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = anipath + '/sae_6-31gd.dat'
nnfdir   = anipath + '/networks/'

# Construct pyNeuroChem class
nc = pync.conformers(cnstfile, saefile, nnfdir, 0)

spc = ['C','C']

xyz = np.vstack([np.array([[0.0,0.0,0.0],[0.0,0.0,float(i*0.01)]],dtype=np.float32) for i in range(600)]).reshape(600,2,3)
print(xyz)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=spc)

ani1 = hdn.hatokcal * nc.energy()

xv = xyz[:,0,:]-xyz[:,1,:]
x = np.linalg.norm(xv, axis=1)

ani1 = ani1 - ani1.min()

plt.plot(x,ani1,'r-.', color='blue', label= 'ANI-c08e',linewidth=4)

plt.title("C-C dissociation curve")
#plt.xlabel('Conformation Pair (Count 49)')
plt.xlabel('Atomic distance ($\AA$)')
plt.ylabel('Calculated dE ($kcal*{mol}^-1$)')
plt.legend(bbox_to_anchor=(0.65, 0.975), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()