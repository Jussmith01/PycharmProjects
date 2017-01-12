import sys
sys.path.append('/home/jujuman/Gits/NeuroChem/pyase_interface')
from neurochemToASE_iface import NeuroChem2ASE
import pyNeuroChem as pync

import numpy as np
import time

# ASE
from ase.io import read, write
from ase.optimize import BFGS, LBFGS
from ase.vibrations import Vibrations

# Set required files for pyNeuroChem
anipath  = '/home/jujuman/Research/wB97X-631gd-train-highgarden/train_08-a3.1A_r4.6_AEV384_1'
cnstfile = anipath + '/rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = anipath + '/../sae_6-31gd.dat'
nnfdir   = anipath + '/networks/'

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile, saefile, nnfdir, 0)

#Read geometry from xyz file
geometry = read('/home/jujuman/Gits/ASE_ANI/examples/data/bz.xyz')

# Setup ANI and calculate single point energy
geometry.set_calculator(NeuroChem2ASE(nc))
e = geometry.get_potential_energy()
print('Total energy', e, 'eV')

# Get the forces
print('Forces:')
print(geometry.get_forces())

# Geometry optimization with BFGS
start_time = time.time()
dyn = LBFGS(geometry)
dyn.run(fmax=0.0001)
print('[ANI Total time:', time.time() - start_time, 'seconds]')

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
