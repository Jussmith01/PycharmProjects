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

xyz = [[0.00000,0.00000,0.37124,0.00000, 0.00000,-0.37124]
	  ,[0.00000,0.372635,0.00000,0.00000,-0.372635, 0.00000]]
	  #,[0.00000,0.00000,0.41000,0.00000, 0.00000,-0.41000]]

typ = ['H','H']

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile, saefile, nnfdir, 1)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ)

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(nc.getNumConfs()) )

# Compute Energies of Conformations
E = np.array(nc.computeEnergies())
F = np.array(nc.computeAnalyticalForces())

print ('-----------------DATA---------------')
#E = gt.hatokcal*E
print (E)
print (F)