__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import numpy as np
import graphtools as gt

import matplotlib.pyplot as plt

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-1-ntwk/'
cnstfile = wkdir + 'rHCNO-4.6A_32-3.1A_a8-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/smallAEV_testing/train_384-256-128-64-1_c08e/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

dtdir = '/home/jujuman/Research/'
xyz,typ,Eact,tmp    = gt.readncdat(dtdir + 'gdb11_s02-1_test.dat',np.float32)
#xyz,typ,Na = gt.readxyz(file)

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile, saefile, nnfdir, 0)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ)
#nc.setMolecule(coords=xyz,types=typ)

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(nc.getNumConfs()) )

#O = nc.optimize(conv=0.000001)

#print(O)

# Compute Energies of Conformations
E = nc.energy()

print('\nE:')
print(E)

#F = nc.force()

#print('\nF:')
#print(F)

AE = nc.aenergies()

print('\nAE:')
print(AE)
print(typ)

AV1 = nc.activations(atom_idx=0, layer_idx=2, molec_idx=0)

print('\nAV1:')
print(AV1)

AV2 = nc.activations(atom_idx=0, layer_idx=2, molec_idx=1)

print('\nAV2:')
print(AV2)

AV3 = AV1 - AV2

print('\nAV3:')
print(AV3)

IDX = np.arange(0,AV3.shape[0],1,dtype=float) + 1
plt.plot (IDX,AV1,marker='o', color='red',  label='ANI-X',  linewidth=2)
plt.plot (IDX,AV2,marker='o', color='blue',  label='ANI-X',  linewidth=2)

plt.title("C10H20 - ANI vs DFT")

plt.ylabel('E cmp (kcal/mol)')
plt.xlabel('E act (kcal/mol)')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.,fontsize=16)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

plt.show()