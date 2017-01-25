__author__ = 'jujuman'

# Import pyNeuroChem
import sys
sys.path.append('/home/jujuman/Gits/NeuroChem/src-python')
from ase_interface import NeuroChem2ASE
import pyNeuroChem as pync

import numpy as np
import  ase
import time
from ase import units
from ase.io import read, write
from ase.optimize import BFGS, LBFGS

# Set required files for pyNeuroChem
anipath  = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk'
cnstfile = anipath + '/rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = anipath + '/sae_6-31gd.dat'
nnfdir   = anipath + '/networks/'

# Construct pyNeuroChem class
nc = pync.molecule(cnstfile, saefile, nnfdir, 0)

mol = read('/home/jujuman/Research/ChiralTest/mol1.xyz')

mol.set_calculator(NeuroChem2ASE(nc))

ei = mol.get_potential_energy()*23.061
print("Initial Energy: ",ei)

# O of Conformations
print("Optimizing...")
start_time = time.time()
dyn = LBFGS(mol)
dyn.run(fmax=0.00001)
print('[ANI Optimization - Total time:', time.time() - start_time, 'seconds]')

ef = mol.get_potential_energy()*23.061
print("Final Energy: ",ef)


# Write visualization of molecule
f = open("optmol.xyz",'w')
f.write('\n' + str(len(mol)) + '\n')
for i in mol:
    f.write(str(i.symbol) + ' ' + str(i.x) + ' ' + str(i.y) + ' ' + str(i.z) + '\n')
f.close()