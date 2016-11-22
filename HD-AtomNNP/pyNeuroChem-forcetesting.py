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

xyz = [[0.00000,0.00000,0.37124,0.00000, 0.00000,-0.37124]
	  ,[0.00000,0.372635,0.00000,0.00000,-0.372635, 0.00000]]
	  #,[0.00000,0.00000,0.41000,0.00000, 0.00000,-0.41000]]

xyz = [[0.0, 0.0, 0.118604, 0.0, 0.67007291, -0.40437055, 0.0, -0.75981897, -0.474415], [0.0, 0.0, 0.118604, 0.0, 0.67007291, -0.4043705, 0.0, -0.75981897, -0.474415]]

typ = ['O','H','H']

xyz,typ,Na = gt.readxyz('/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/Retinol/optimization/opt_test_numeric.xyz')

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile, saefile, nnfdir, 0)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ[0])

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