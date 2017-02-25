__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import numpy as np
import hdnntools as gt
import os

# Set required files for pyNeuroChem
anipath  = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk'
cnstfile = anipath + '/rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = anipath + '/sae_6-31gd.dat'
nnfdir   = anipath + '/networks/'

# Construct pyNeuroChem class
nc = pync.conformers(cnstfile, saefile, nnfdir, 0)

dir = "/home/jujuman/Research/AEChargeTesting"
files = os.listdir(dir + '/xyz')

of = open(dir+'/localatomenergies.txt','w')

for i in files:
    print(i)
    of.write(i+'\n')
    xyz,typ,Na = gt.readxyz2(dir + "/xyz/" + i)

    # Set the conformers in NeuroChem
    nc.setConformers(confs=xyz,types=list(typ))

    # Compute Energies of Conformations
    E = nc.energy()

    # Atomic energy return
    AE = np.copy(nc.aenergies(sae=False))
    of.write(str(Na) + '\n')
    of.write(np.array_str(np.array(typ)) + '\n')
    of.write(np.array_str(AE)+'\n')

of.close()