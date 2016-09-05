__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
import numpy as np
import matplotlib.pyplot as plt

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Research/trainingcases/wB97X-631gd-train-highgarden/'

#Network 1 Files
cnstfile1 = wkdir + 'train_08-a3.1A_r4.5_dn1/rHCNO-4.5A_32-3.1A_a8-8.params'
saefile1  = wkdir + 'sae_6-31gd.dat'
nnfdir1   = wkdir + 'train_08-a3.1A_r4.5_dn1/networks/'

# Network 2 Files
cnstfile2 = wkdir + 'nw-64-64-64-32/train_07-a3.1A_r4.5/rHCNO-4.5A_32-3.1A_a8-8.params'
saefile2  = wkdir + 'sae_6-31gd.dat'
nnfdir2   = wkdir + 'nw-64-64-64-32/train_07-a3.1A_r4.5/networks/'

# Construct pyNeuroChem classes
nc1 = pync.pyNeuroChem(cnstfile1,saefile1,nnfdir1,0)
nc2 = pync.pyNeuroChem(cnstfile2,saefile2,nnfdir2,0)

xyz,typ,Eact = gt.readncdat('/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_testdata/formamide_dhl_test.dat')

Eact = np.array(Eact)

# Set the conformers in NeuroChem
nc1.setConformers(confs=xyz,types=typ)
nc2.setConformers(confs=xyz,types=typ)

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(nc1.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(nc1.getNumConfs()) )

# Compute Forces of Conformations
print('Computing energies 1...')
Ecmp1 = np.array( nc1.computeEnergies() )
print('Computation complete 1.')

print('Computing energies 2...')
Ecmp2 = np.array( nc2.computeEnergies() )
print('Computation complete 2.')

Ecmp1 = gt.hatokcal * Ecmp1
Ecmp2 = gt.hatokcal * Ecmp2
Eact = gt.hatokcal * Eact

#Ecmp1 = Ecmp1 - gt.calculatemean(Ecmp1)
#Ecmp2 = Ecmp2 - gt.calculatemean(Ecmp2)
#Eact  = Eact  - gt.calculatemean(Eact)

IDX = np.arange(0,Eact.shape[0],1,dtype=float)

rmse1 = gt.calculaterootmeansqrerror(Eact,Ecmp1)
rmse2 = gt.calculaterootmeansqrerror(Eact,Ecmp2)

plt.plot   (IDX, Eact, color='black',      label='DFT',linewidth=2)
plt.scatter(IDX, Eact, color='black',  linewidth=3)

plt.plot   (IDX, Ecmp1, '--', color='red', label='ANN - c08b RMSE: ' + "{:.6f}".format(rmse1) + "kcal/mol",  linewidth=2)
plt.scatter(IDX, Ecmp1, color='red',  linewidth=3)

plt.plot   (IDX, Ecmp2, '--', color='green',label='ANN - c07b RMSE: ' + "{:.6f}".format(rmse2) + "kcal/mol",  linewidth=2)
plt.scatter(IDX, Ecmp2, color='green',linewidth=3)

plt.title("Formamide Dihedral Scan")

plt.ylabel('E (kcal/mol)')
plt.xlabel('Scan Coordinate ($\AA$)')
plt.legend(bbox_to_anchor=(0.3, 0.975), loc=2, borderaxespad=0.,fontsize=12)

plt.show()