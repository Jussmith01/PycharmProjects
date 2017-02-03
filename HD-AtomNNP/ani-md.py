import sys
sys.path.append('/home/jujuman/Gits/NeuroChem/src-python')
from ase_interface import NeuroChem2ASE
import pyNeuroChem as pync

import numpy as np
import  ase
import time
#from ase.build import molecule
#from ase.neb import NEB
#from ase.calculators.mopac import MOPAC
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
from ase import units

from ase.optimize.fire import FIRE as QuasiNewton

from ase.md.nvtberendsen import NVTBerendsen
from ase.md import MDLogger

#from ase.neb import NEBtools
from ase.io import read, write
from ase.optimize import BFGS, LBFGS

#import matplotlib
#import matplotlib as mpl

#import matplotlib.pyplot as plt

#import seaborn as sns
#%matplotlib inline

# Set required files for pyNeuroChem
anipath  = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk'
cnstfile = anipath + '/rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = anipath + '/sae_6-31gd.dat'
nnfdir   = anipath + '/networks/'

# Construct pyNeuroChem class
nc = pync.molecule(cnstfile, saefile, nnfdir, 0)

#bz = read('/home/jujuman/Dropbox/ChemSciencePaper.AER/Poster-GTC-May-2017/Timings/min1.xyz')
bz = read('/home/jujuman/Downloads/1le1.pdb')

#L = 16.7
#bz.set_cell(([[L,0,0],[0,L,0],[0,0,L]]))
#bz.set_pbc((True, True, True))

bz.set_calculator(NeuroChem2ASE(nc))

start_time = time.time()
dyn = LBFGS(bz)
dyn.run(fmax=0.001)
print('[ANI Total time:', time.time() - start_time, 'seconds]')

# Write visualization of molecule
f = open("optmol.xyz",'w')
f.write('\n' + str(len(bz)) + '\n')
for i in bz:
    f.write(str(i.symbol) + ' ' + str(i.x) + ' ' + str(i.y) + ' ' + str(i.z) + '\n')
f.close()
# Temperature
T = 300.0

# We want to run MD with constant energy using the Langevin algorithm
# with a time step of 5 fs, the temperature T and the friction
# coefficient to 0.02 atomic units.
dyn = Langevin(bz, 0.25 * units.fs, T * units.kB, 0.05)
#dyn = NVTBerendsen(bz, 0.5 * units.fs, 300.0, taut=3.0*1000*units.fs)

start_loop_time = time.time()

mdcrd = open("mdcrd.xyz",'w')
temp = open("temp.dat",'w')
def printenergy(a=bz,b=mdcrd,d=dyn,t=temp,loop_tm=start_loop_time):  # store a reference to atoms in the
    """Function to print the potential, kinetic and total energy."""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Step %i - Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV Time: %.3fs' % (d.get_number_of_steps(),epot, ekin, ekin / (1.5 * units.kB), epot + ekin, time.time() - loop_tm))

    loop_tm = time.time()

    t.write(str(d.get_number_of_steps()) + ' ' + str(ekin / (1.5 * units.kB)) + ' ' + str(epot) + ' ' +  str(ekin) + ' ' + str(epot + ekin) + '\n')
    b.write('\n' + str(len(a)) + '\n')
    for i in a:
        b.write(str(i.symbol) + ' ' + str(i.x) + ' ' + str(i.y) + ' ' + str(i.z) + '\n')

dyn.attach(printenergy, interval=5)
#dyn.attach(MDLogger(dyn, bz, 'bz_md_NVT_10ps_1fs.log', header=True, stress=False,
#           peratom=False, mode="w"), interval=50)

printenergy()
print("test")

start_time = time.time()
dyn.run(4*1000*1000) # Do 1ns of MD
print('[ANI Total time:', time.time() - start_time, 'seconds]')

mdcrd.close()
temp.close()
