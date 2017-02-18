import sys
import time

# Numpy
import numpy as np

# Neuro Chem
from ase_interface import ANI
import pyNeuroChem as pync

import  ase
#from ase.build import molecule
#from ase.neb import NEB
#from ase.calculators.mopac import MOPAC
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
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

#bz = read('C_100.xyz')
bz = read('/home/jujuman/Dropbox/ChemSciencePaper.AER/JustinsDocuments/Poster-GTC-May-2017/Timings/2naz_neutralized_manual2.pdb')


#L = 22.
#bz.set_cell(([[L,0,0],[0,L,0],[0,0,L]]))
#bz.set_pbc((True, True, True))

bz.set_calculator(ANI(False))
bz.calc.setnc(nc)

#start_time = time.time()
#dyn = LBFGS(bz)
#dyn.run(fmax=0.001)
#dyn = BFGS(bz)
#dyn.run(fmax=0.1)
#print('[ANI Total time:', time.time() - start_time, 'seconds]')

# Write visualization of molecule
f = open("optmol_begin.xyz",'w')
f.write('\n' + str(len(bz)) + '\n')
for i in bz:
    f.write(str(i.symbol) + ' ' + str(i.x) + ' ' + str(i.y) + ' ' + str(i.z) + '\n')
f.close()

# Temperature
T = 100.0

# We want to run MD with constant energy using the Langevin algorithm
# with a time step of 5 fs, the temperature T and the friction
# coefficient to 0.02 atomic units.
dyn = Langevin(bz, 0.25 * units.fs, T * units.kB, 0.5)
#dyn = NPTBerendsen(bz, 0.2 * units.fs, temperature=500.,taut=0.1*1000*units.fs, pressure = 1.01325, taup=1.0*1000*units.fs, compressibility=4.57e-5 )
#dyn = NVTBerendsen(bz, 0.5 * units.fs, 200., taut=3.0*1000*units.fs)

mdcrd = open("mdcrd.xyz",'w')
temp = open("temp.dat",'w')

#print(bz.get_positions())

#c=bz.get_positions(wrap=True)
#print(c)
#print(c.x)
#for i in c:
#    print(c.x)
#    print(str(c.x) + ' ' + str(c.y) + ' ' + str(c.z) + '\n')


def printenergy(a=bz,b=mdcrd,d=dyn,t=temp):  # store a reference to atoms in the
    """Function to print the potential, kinetic and total energy."""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Step %i - Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (d.get_number_of_steps(),epot, ekin, ekin / (1.5 * units.kB), epot + ekin))
    t.write(str(d.get_number_of_steps()) + ' ' + str(ekin / (1.5 * units.kB)) + ' ' + str(epot) + ' ' +  str(ekin) + ' ' + str(epot + ekin) + '\n')
    b.write('\n' + str(len(a)) + '\n')
    c=a.get_positions(wrap=True)
    for j,i in zip(a,c):
        b.write(str(j.symbol) + ' ' + str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) + '\n')

dyn.attach(printenergy, interval=50)
#dyn.attach(MDLogger(dyn, bz, 'bz_md_NVT_10ps_1fs.log', header=True, stress=False,
#           peratom=False, mode="w"), interval=50)

#traj = Trajectory('killerapp.traj', 'w', atoms=bz)
#dyn.attach(traj.write, interval=100)

printenergy()

for i in range(0,300,2):
    print("Heating to",float(i), "K...")
    dyn.set_temperature(float(i) * units.kB)
    start_time = time.time()
    dyn.run(4000) # Do 5ps of MD
    print('[ANI Total time:', time.time() - start_time, 'seconds]')

dyn.run(40000)  # Do 5ps of MD

#for i in range(0,5000,1):
#    print("Heating to",float(i), "K...")
#    dyn.set_temperature(float(i) * units.kB)
#    dyn.run(500) # Do 5ps of MD

dyn.run(10000000) # Do 5ps of MD

dyn = LBFGS(bz)
dyn.run(fmax=0.0001)
#dyn = BFGS(bz)
#dyn.run(fmax=0.1)
#print('[ANI Total time:', time.time() - start_time, 'seconds]')

# Write visualization of molecule
f = open("optmol_final.xyz",'w')
f.write('\n' + str(len(bz)) + '\n')
for i in bz:
    f.write(str(i.symbol) + ' ' + str(i.x) + ' ' + str(i.y) + ' ' + str(i.z) + '\n')
f.close()

mdcrd.close()
temp.close()
