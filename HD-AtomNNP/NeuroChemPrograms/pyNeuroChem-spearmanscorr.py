__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
import numpy as np

from scipy import stats as st

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/train_07/'

#Network Parameter Files
cnstfile = wkdir + 'rHCNO-4.5A_32-3.1A_a8-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

# Construct pyNeuroChem classes
nc = pync.pyNeuroChem(cnstfile,saefile,nnfdir,0)

xyz,typ,Eact = gt.readncdat('/home/jujuman/Dropbox/Research/ChemSciencePaper/TestCases/Atomoxetine/Atomoxetine_conformersC_test.dat')

Eact = np.array(Eact)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ)

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(nc.getNumConfs()) )

# Compute Forces of Conformations
print('Computing energies 1...')
Ecmp = np.array( nc.computeEnergies() )
print('Computation complete 1.')

Ecmp = gt.hatokcal * Ecmp
Eact = gt.hatokcal * Eact

print(st.spearmanr(Ecmp,Eact))
