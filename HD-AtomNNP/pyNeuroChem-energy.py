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

#dtdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_01/data/'
#xyz,typ,Eact,tmp    = gt.readncdat(dtdir + 'gdb11_s01-1_test.dat',np.float32)

xyz,typ,Na = gt.readxyz('/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/test_for_dipoles/nh2bz.xyz')

print(xyz)

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile, saefile, nnfdir, 0)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ)
#nc.setMolecule(coords=xyz,types=typ)

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

AE = nc.aenergies()

print('\nAE:')
print(AE)
print(typ)
