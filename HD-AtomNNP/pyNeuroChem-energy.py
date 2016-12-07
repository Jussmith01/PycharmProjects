__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import numpy as np
import graphtools as gt

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-1-ntwk/'
cnstfile = wkdir + 'rHCNO-4.6A_32-3.1A_a8-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'



#for i in range(0,3):

print (xyz2)

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile, saefile, nnfdir, 0)

# Set the conformers in NeuroChem
#nc.setConformers(confs=xyz,types=typ)
nc.setMolecule(coords=xyz2,types=typ)

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(nc.getNumConfs()) )

# Compute Energies of Conformations
E = nc.energy()

print('\nE:')
print(E)

F = nc.force()

print('\nF:')
print(F)

O1 = nc.optimize(conv=0.00001)

print (O1.dtype)
nc.setMolecule(coords=O1,types=typ)

E1 = nc.energy()
#print('\nE')
print(E1)

