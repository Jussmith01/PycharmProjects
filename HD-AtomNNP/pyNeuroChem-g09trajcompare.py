__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
import numpy as np
import matplotlib.pyplot as plt
import time as tm
from scipy import stats as st

def sortbyother(Y, X):
    xy = zip(X, Y)
    xy = sorted(xy, key=lambda x: x[0])
    X, Y = zip(*xy)
    return np.array(Y)

# Set required files for pyNeuroChem
wkdir1 = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-1-ntwk/'
wkdir1 = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/train_08_9/'

#Network 1 Files
cnstfile1 = wkdir1 + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile1  = wkdir1 + 'sae_6-31gd.dat'
nnfdir1   = wkdir1 + 'networks/'

dtdir = '/home/jujuman/Scratch/Research/ANN-Test-Data/MDinG09testing/'
dtdir = '/home/jujuman/Dropbox/Research/MDTests/'

#xyz,typ,Eact = gt.readncdat('../data_irc.dat')
xyz,typ,Eact = gt.readg09trajdat(dtdir + 'ranolazine_md.log',np.float32)
#xyz,typ,Eact = gt.readg09trajdat(dtdir + 'md.log')
#xyz2,typ2,Eact2,tmp = gt.readncdat(dtdir + 'ranolazine_dftb.dat')
typ = typ[0]

#gt.writexyzfile('traj.xyz',xyz,typ)

print(typ)
print(len(xyz))
print(Eact)

Eact = np.array(Eact)
#Eact2 = np.array(Eact2)

# Construct pyNeuroChem classes
nc1 = pync.pyNeuroChem(cnstfile1,saefile1,nnfdir1,0)

# Set the conformers in NeuroChem
nc1.setConformers(confs=xyz,types=typ)

# Print some data from the NeuroChem
print( '1) Number of Atoms Loaded: ' + str(nc1.getNumAtoms()) )
print( '1) Number of Confs Loaded: ' + str(nc1.getNumConfs()) )

# Compute Forces of Conformations
print('Computing energies 1...')
_t1b = tm.time()
Ecmp1 = np.array( nc1.energy() )
print('Computation complete 1. Time: ' + "{:.4f}".format((tm.time() - _t1b) * 1000.0)  + 'ms')

#n = 0
#m = 9
#Ecmp1 = gt.hatokcal * Ecmp1[n:m]
#Eact  = gt.hatokcal * Eact[n:m]
Ecmp1 = gt.hatokcal * Ecmp1
Eact  = gt.hatokcal * Eact
#Eact2  = gt.hatokcal * Eact2

IDX = np.arange(0,Eact.shape[0],1,dtype=float) + 1

#Ecmp1 = Ecmp1 - Ecmp1.min()
#Eact  = Eact  - Eact.min()
#Eact2  = Eact2  - Eact2.min()

rmse1 = gt.calculaterootmeansqrerror(Eact,Ecmp1)
#rmse2 = gt.calculaterootmeansqrerror(Eact,Eact2)

print ( "Spearman corr. 1: " + "{:.3f}".format(st.spearmanr(Ecmp1,Eact)[0]) )
#print ( "Spearman corr. 2: " + "{:.3f}".format(st.spearmanr(Eact2,Eact)[0]) )

#plt.plot   (IDX, Eact, '-', marker=r'o', color='black',label='DFT',  linewidth=2, markersize=10)
#plt.plot   (IDX, Ecmp1, ':', marker=r'D', color='red',  label='ANI-1 RMSE: ' + "{:.7f}".format(rmse1) + ' kcal/mol',  linewidth=4, markersize=8)
plt.plot   (IDX, Eact, '-', color='black',label='DFT',  linewidth=3)
plt.plot   (IDX, Ecmp1, '--', color='red',  label='ANI-1 RMSE: ' + "{:.3f}".format(rmse1) + ' kcal/mol',  linewidth=3)
#plt.plot   (IDX, Eact2, '--', color='blue',  label='DFTB   RMSE: ' + "{:.3f}".format(rmse2) + ' kcal/mol',  linewidth=2)

plt.title("BO MD Trajectory - Ranolazine")

plt.ylabel('$\Delta$E calculated (kcal/mol)')
plt.xlabel('Step number')
#plt.legend(bbox_to_anchor=(0.35, 0.98), loc=2, borderaxespad=0.,fontsize=16)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

plt.show()
