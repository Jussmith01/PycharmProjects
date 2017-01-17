__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import numpy as np
import graphtools as gt
import time

import matplotlib.pyplot as plt

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

dtdir = '/home/jujuman/Dropbox/ChemSciencePaper.AER/Poster-GTC-May-2017/Timings/'
#xyz,typ,Eact,tmp    = gt.readncdat(dtdir + 'gdb11_s01-1_test.dat',np.float32)
xyz,typ,Na = gt.readxyz(dtdir + 'water.xyz')

print (xyz)

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile, saefile, nnfdir, 0)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ)
#nc.setMolecule(coords=xyz,types=typ)

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(nc.getNumConfs()) )

#print ('Optimizing...')
#O = nc.optimize(conv=0.000001,max_iter=250)
#print(O)

#print(len(typ))
#for i in range(0,len(typ)):
#    print (typ[i] + ' ' + "{:.10f}".format(O[i,0]) + ' ' + "{:.10f}".format(O[i,1]) + ' ' + "{:.10f}".format(O[i,2]))

# Compute Energies of Conformations
print('Calculating energies and forces...')
start_time = time.time()
E = nc.energy()
#F = nc.force()
print('[ANI Total time:', time.time() - start_time, 'seconds]')

print('\nE:')
print(E)

#print('\nF:')
#print(F)
'''
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
'''