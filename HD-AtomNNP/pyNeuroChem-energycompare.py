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
wkdir1    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/smallAEV_testing/train_384-256-128-64-1_c08e/'

#Network 1 Files
cnstfile1 = wkdir1 + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile1  = wkdir1 + 'sae_6-31gd.dat'
nnfdir1   = wkdir1 + 'networks/'

dtdir = '/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/C10H20Isomers/'

#xyz,typ,Eact = gt.readncdat('../data_irc.dat')
xyz,typ,Eact,tmp    = gt.readncdat(dtdir + 'isomer_structures_DFT.dat')
#xyz2,typ2,Eact2,tmp = gt.readncdat(dtdir + 'isomer_structures_DFTB.dat')
#xyz3,typ3,Eact3,tmp = gt.readncdat(dtdir + 'isomer_structures_PM6.dat')

#xyz = [xyz[0],xyz[1]]
#xyz2 = [xyz2[0],xyz2[1]]
#xyz3 = [xyz3[0],xyz3[1]]

xyz = np.asarray(xyz,dtype=np.float32)
xyz = xyz.reshape((xyz.shape[0],len(typ),3))

Eact = np.array(Eact)
#Eact2 = np.array(Eact2)
#Eact3 = np.array(Eact3)

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
Ecmp1 = nc1.energy()
print('Computation complete 1. Time: ' + "{:.4f}".format((tm.time() - _t1b) * 1000.0)  + 'ms')


print(Ecmp1)

'''
# Construct pyNeuroChem classes
nc2 = pync.pyNeuroChem(cnstfile2,saefile2,nnfdir2,0)

# Set the conformers in NeuroChem
nc2.setConformers(confs=xyz,types=typ)

# Print some data from the NeuroChem
print( '2) Number of Atoms Loaded: ' + str(nc2.getNumAtoms()) )
print( '2) Number of Confs Loaded: ' + str(nc2.getNumConfs()) )

# Compute Forces of Conformations
print('Computing energies 1...')
_t2b = tm.time()
Ecmp2 = np.array( nc2.computeEnergies() )
print('Computation complete 1. Time: ' + "{:.4f}".format((tm.time() - _t2b) * 1000.0)  + 'ms')

print(Ecmp2)
'''

n = 0
m = 13
Ecmp1 = gt.hatokcal * Ecmp1
#Ecmp2 = gt.hatokcal * Ecmp2[n:m]
Eact  = gt.hatokcal * Eact
#Eact2 = gt.hatokcal * Eact2[n:m]
#Eact3 = gt.hatokcal * Eact3[n:m]

dE = abs(Eact - Ecmp1)

for i in range (0,dE.shape[0]):
    print ( str(i) + ' ' + str(dE[i]) )

IDX = np.arange(0,Eact.shape[0],1,dtype=float) + 1
IDX2 = sortbyother(IDX, Eact)
#print IDX2

Ecmp1 = sortbyother(Ecmp1, Eact)
#Ecmp2 = sortbyother(Ecmp2, Eact)
#Eact2 = sortbyother(Eact2, Eact)
#Eact3 = sortbyother(Eact3, Eact)
Eact  = np.sort( Eact )

#print (Eact2.min())

#Ecmp1 = Ecmp1 - Ecmp1.min()
#Ecmp2 = Ecmp2 - Ecmp2.min()
#Eact  = Eact  - Eact.min()
#Eact2 = Eact2 - Eact2.min()
#Eact3 = Eact3 - Eact3.min()

rmse1 = gt.calculaterootmeansqrerror(Eact,Ecmp1)
#rmse2 = gt.calculaterootmeansqrerror(Eact,Ecmp2)
#rmse3 = gt.calculaterootmeansqrerror(Eact,Eact2)
#rmse4 = gt.calculaterootmeansqrerror(Eact,Eact3)

print ( "Spearman corr. 1: " + "{:.3f}".format(st.spearmanr(Ecmp1,Eact)[0]) )
#print ( "Spearman corr. 2: " + "{:.3f}".format(st.spearmanr(Ecmp2,Eact)[0]) )
#print ( "Spearman corr. 3: " + "{:.3f}".format(st.spearmanr(Eact2,Eact)[0]) )
#print ( "Spearman corr. 4: " + "{:.3f}".format(st.spearmanr(Eact3,Eact)[0]) )

plt.plot   (IDX, Eact, color='black',label='DFT',  linewidth=2)

plt.scatter   (IDX, Ecmp1, color='red',  label='ANI-X RMSE: ' + "{:.3f}".format(rmse1) + ' kcal/mol',  linewidth=2)
#plt.plot   (IDX, Ecmp2, ':', marker=r'D', color='orange',  label='ANN - c08c RMSE: ' + "{:.3f}".format(rmse2) + ' kcal/mol',  linewidth=2, markersize=5)
#plt.plot   (IDX, Eact2, ':', marker=r'v', color='blue', label='DFTB  RMSE: ' + "{:.2f}".format(rmse3) + ' kcal/mol',  linewidth=2, markersize=8)
#plt.plot   (IDX, Eact3, ':', marker=r'*', color='orange',label='PM6   RMSE: ' + "{:.2f}".format(rmse4) + ' kcal/mol',  linewidth=2, markersize=9)
#plt.plot   (IDX, Ecmp1, color='red',  label='ANN - c08b RMSE: ' + "{:.3f}".format(rmse1) + ' kcal/mol',  linewidth=5)
#plt.plot   (IDX, Eact2, color='blue', label='DFTB  RMSE: ' + "{:.2f}".format(rmse3) + ' kcal/mol',  linewidth=5)
#plt.plot   (IDX, Eact3, color='orange',label='PM6   RMSE: ' + "{:.2f}".format(rmse4) + ' kcal/mol',  linewidth=5)


#plt.plot   (IDX, Eact, color='black',  label='DFT',linewidth=3)
#plt.scatter(IDX, Eact, marker='o' , color='black',  linewidth=4)

#plt.title("300K NMS structures of\nNME-Gly-Pro-Hyp-Gly-Ala-Gly-ACE")
plt.title("C10H20 - ANI vs DFT")

plt.ylabel('E cmp (kcal/mol)')
plt.xlabel('E act (kcal/mol)')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.,fontsize=16)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

plt.show()
