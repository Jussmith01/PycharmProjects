__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
import numpy as np

from scipy import stats as st

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Research/trainingcases/wB97X-631gd-train-highgarden/'

#Network 1 Files
cnstfile = wkdir + 'train_08-a3.1A_r4.5_dn1/rHCNO-4.5A_32-3.1A_a8-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'train_08-a3.1A_r4.5_dn1/networks/'

# Construct pyNeuroChem classes
nc = pync.pyNeuroChem(cnstfile,saefile,nnfdir,0)

xyz,typ,Eact = gt.readncdat('/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_testdata/')

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
