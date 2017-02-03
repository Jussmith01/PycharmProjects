import sys
sys.path.append('/home/jujuman/Gits/NeuroChem/src-python')
from ase_interface import NeuroChem2ASE
import pyNeuroChem as pync

import numpy as np
import time

# ASE
from ase.io import read, write
from ase.optimize import BFGS, LBFGS
from ase.vibrations import Vibrations

# Set required files for pyNeuroChem
anipath  = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk'
cnstfile = anipath + '/rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = anipath + '/sae_6-31gd.dat'
nnfdir   = anipath + '/networks/'

# Construct pyNeuroChem class
nc = pync.molecule(cnstfile, saefile, nnfdir, 0)

#Read geometry from xyz file
geometry = read('/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_begdb/begdb-h2oclusters/xyz/4197_water6CB2.xyz')

# Setup ANI and calculate single point energy
geometry.set_calculator(NeuroChem2ASE(nc))
e = geometry.get_potential_energy()
print('Total energy', e, 'eV')

# Get the forces
print('Forces:')
print(geometry.get_forces())

f = open("pre-normoptmol.xyz",'w')
f.write('\n' + str(len(geometry)) + '\n')
for i in geometry:
    f.write(str(i.symbol) + ' ' + str(i.x) + ' ' + str(i.y) + ' ' + str(i.z) + '\n')
f.close()

# Geometry optimization with BFGS
start_time = time.time()
dyn = LBFGS(geometry)
dyn.run(fmax=0.00001)
print('[ANI Total time:', time.time() - start_time, 'seconds]')

f = open("post-normoptmol.xyz",'w')
f.write('\n' + str(len(geometry)) + '\n')
for i in geometry:
    f.write(str(i.symbol) + ' ' + str(i.x) + ' ' + str(i.y) + ' ' + str(i.z) + '\n')
f.close()

# Calc minimized energies
e = geometry.get_potential_energy()
print('Total energy', e, 'eV')

# Get the forces again
print('Forces:')
print(geometry.get_forces())

# Run the vibrational analysis
vib = Vibrations(geometry, nfree=2)
vib.run()
print( vib.summary() )
print( vib.get_zero_point_energy() )

# Get freqs
f = vib.get_frequencies()

# Print modes + freq
for i in range(0,len(f)):
    print('Mode(' + str(i) +') Freq: ' + "{:.7e}".format(f[i]) )
    print(vib.get_mode(i))

'''
# We will alse need to use these functions to get thermo-chem data
thermo = IdealGasThermo(vib_energies=vib_energies,
                       potentialenergy=potentialenergy,
                       atoms=atoms,
                       geometry='linear',
                       symmetrynumber=2, spin=0)
G = thermo.get_gibbs_energy(temperature=298.15, pressure=101325.)
'''

# remove mode pckl files
vib.clean()
