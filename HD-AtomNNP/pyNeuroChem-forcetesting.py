__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
import numpy as np

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Research/NeuroChemForceTesting/train_01/'
cnstfile = wkdir + 'rH-3.0A_4-2.5A_a2-2.params'
saefile  = wkdir + '../sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

xyz = [[0.00,0.00,0.50,0.00,0.00,-0.50]
	  ,[0.00,0.52,0.00,0.00,-0.52,0.00]
	  ,[0.00,0.00,0.49,0.00,0.00,-0.49]]

typ = ['H','H']

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile, saefile, nnfdir, 0)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ)

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(nc.getNumConfs()) )

# Compute Energies of Conformations
E = np.array(nc.computeEnergies())


print ('-----------------DATA---------------')
#E = gt.hatokcal*E
print (E)