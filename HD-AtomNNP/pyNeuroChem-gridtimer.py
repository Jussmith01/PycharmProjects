__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
import numpy as np
import matplotlib.pyplot as plt
import time as tm
from scipy import stats as st

def generatemolgrid(typs,Nx,Ny):

    Na = len(typs) * Nx * Ny

    xyz = [[]]
    typ = []

    x = 0.0
    y = 0.0
    z = 0.0

    for i in typs:
        for j in range(Ny):
            for k in range(Nx):
                typ.append(i)
                xyz[0].append(x)
                xyz[0].append(y)
                xyz[0].append(z)
                x += 1.0
            y += 1.0
        z += 1.0

    return xyz,typ,Na

# Set required files for pyNeuroChem
wkdir1    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/train_08_1/'

#Network 1 Files
cnstfile1 = wkdir1 + 'rHCNO-4.5A_32-3.1A_a8-8.params'
saefile1  = wkdir1 + 'sae_6-31gd.dat'
nnfdir1   = wkdir1 + 'networks/'

# Construct pyNeuroChem classes
nc1 = pync.pyNeuroChem(cnstfile1,saefile1,nnfdir1,0)

idx = []
tme = []

for i in range(1,20):
    xyz,typ,Na = generatemolgrid(['H','C','N','O'],10,10)

    # Set the conformers in NeuroChem
    nc1.setConformers(confs=xyz,types=typ)

    # Print some data from the NeuroChem
    print( '1) Number of Atoms Loaded: ' + str(nc1.getNumAtoms()) )
    print( '1) Number of Confs Loaded: ' + str(nc1.getNumConfs()) )

    # Compute Forces of Conformations
    print('Computing energy 1...')
    _t1b = tm.time()
    Ecmp1 = np.array( nc1.computeEnergies() )
    _t1e = tm.time()
    print('Energy computation complete. Time: ' + "{:.4f}".format((_t1e - _t1b) * 1000.0)  + 'ms')

    print('Computing forces 1...')
    _t2b = tm.time()
    F = np.array( nc1.computeAnalyticalForces() )
    _t2e = tm.time()
    print('Force computation complete. Time: ' + "{:.4f}".format((_t2e - _t2b) * 1000.0)  + 'ms')

    print('Force to Energy: ' + "{:.4f}".format((_t2e - _t2b)/(_t1e - _t1b)))

    if i > 10:
        idx.append(nc1.getNumAtoms())
        tme.append((_t1e - _t1b + _t2e - _t2b ) * 1000.0)

print( st.linregress(np.log10(idx),np.log10(tme)) )

plt.plot   (np.log10(idx), np.log10(tme), color='black',linewidth=3)
plt.scatter(np.log10(idx), np.log10(tme), marker='o' , color='black',  linewidth=4)

#plt.title("300K NMS structures of\nNME-Gly-Pro-Hyp-Gly-Ala-Gly-ACE")
plt.title("ANN Compute time scaling by number of atoms")


plt.xlabel('Log10(Number of Atoms)')
plt.ylabel('Log10(Time(ms))')
#plt.legend(bbox_to_anchor=(0.1, 0.95), loc=2, borderaxespad=0.,fontsize=14)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

plt.show()
