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

xyz = np.array([[[0.00000000,  0.07775414,  0.14191920],  [0.00000000,  0.87511754, -0.40200517],  [0.00000000, -0.64261788, -0.50009561]]
	  		   ,[[0.00000000,  0.07782031,  0.14223997],  [0.00000000,  0.87485045, -0.40213168],  [0.00000000, -0.64241695, -0.50028986]]
	  		   ,[[0.00000000,  0.07782238,  0.14223947],  [0.00000000,  0.87484127, -0.40213102],  [0.00000000, -0.64240986, -0.50029004]]
	  		   ,[[0.00000000,  0.07782238,  0.14223947],  [0.00000000,  0.87484127, -0.40213102],  [0.00000000, -0.64240986, -0.50029004]]],dtype=np.float32)

xyz2 = np.array([[0.00000000,  0.2775414,  0.14191920],  [0.00000000,  0.87511754, -0.40200517],  [0.00000000, -0.64261788, -0.50009561]],dtype=np.float32)


typ = ['O','H','H']
Na  = [3]

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

