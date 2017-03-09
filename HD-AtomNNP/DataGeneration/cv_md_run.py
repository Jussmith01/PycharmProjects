import sys
import time

# Numpy
import numpy as np

# Neuro Chem
from ase_interface import ANI
import pyNeuroChem as pync
import hdnntools as hdt

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


#--------------Parameters------------------
wkdir = '/home/jujuman/Research/CrossValidation/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile = wkdir + 'sae_6-31gd.dat'

At = ['C', 'O', 'N'] # Hydrogens added after check

T = 300.0
dt = 0.25

stdir = '/home/jujuman/Research/CrossValidation/MD_CV/'

#-------------------------------------------

# Construct pyNeuroChem classes
print('Constructing CV network list...')
ncl =  [pync.molecule(cnstfile, saefile, wkdir + 'cv_c08e_ntw_' + str(l) + '/networks/', 0) for l in range(5)]
print('Complete.')

# Set required files for pyNeuroChem
anipath  = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ccdissotest1-ntwk'
cnstfile = anipath + '/rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = anipath + '/sae_6-31gd.dat'
nnfdir   = anipath + '/networks/'

# Construct pyNeuroChem class
print('Constructing MD network...')
nc = ncl[1]
#nc = pync.molecule(cnstfile, saefile, nnfdir, 0)
print('FINISHED')

mol = read('/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_testdata/specialtest/test.xyz')
#mol = read('/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_begdb/begdb-h2oclusters/xyz/4179_water2Cs.xyz')
#mol = read('/home/jujuman/Research/CrossValidation/MD_CV/benzene.xyz')

#L = 16.0
#bz.set_cell(([[L,0,0],[0,L,0],[0,0,L]]))
#bz.set_pbc((True, True, True))

mol.set_calculator(ANI(False))
mol.calc.setnc(nc)

# We want to run MD with constant energy using the Langevin algorithm
# with a time step of 5 fs, the temperature T and the friction
# coefficient to 0.02 atomic units.
dyn = Langevin(mol, dt * units.fs, T * units.kB, 0.01)

mdcrd = open(stdir + "mdcrd.xyz",'w')
temp = open(stdir + "temp.dat",'w')

dyn.get_time()
def printenergy(a=mol,b=mdcrd,d=dyn,t=temp):  # store a reference to atoms in the
    """Function to print the potential, kinetic and total energy."""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Step %i - Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (d.get_number_of_steps(),epot, ekin, ekin / (1.5 * units.kB), epot + ekin))
    t.write(str(d.get_number_of_steps()) + ' ' + str(d.get_time()) + ' ' + str(ekin / (1.5 * units.kB)) + ' ' + str(epot) + ' ' +  str(ekin) + ' ' + str(epot + ekin) + '\n')
    b.write('\n' + str(len(a)) + '\n')
    c=a.get_positions(wrap=True)
    for j,i in zip(a,c):
        b.write(str(j.symbol) + ' ' + str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) + '\n')

dyn.attach(printenergy, interval=50)
dyn.set_temperature(300.0 * units.kB)
start_time2 = time.time()

# get the chemical symbols
spc = mol.get_chemical_symbols()

f = open(stdir + 'md-peptide-cv.dat','w')
for i in range(1000000):
    dyn.run(20)  # Do 5ps of MD

    xyz = np.array(mol.get_positions(), dtype=np.float32).reshape(len(spc), 3)
    energies = np.zeros((5), dtype=np.float64)
    N = 0
    for comp in ncl:
        comp.setMolecule(coords=xyz, types=list(spc))
        energies[N] = comp.energy()[0]
        N = N + 1

    energies = hdt.hatokcal * energies

    f.write("{:.7f}".format((dyn.get_number_of_steps() * dt)/1000.0) +
            ' ' + "{:.7f}".format(energies[0]) +
            ' ' + "{:.7f}".format(energies[1]) +
            ' ' + "{:.7f}".format(energies[2]) +
            ' ' + "{:.7f}".format(energies[3]) +
            ' ' + "{:.7f}".format(energies[4]) +
            ' ' + "{:.7f}".format(np.std(energies)) + '\n')

    output = '  ' + str(i) + ' (' + str(len(spc)) + ') : stps=' + str(dyn.get_number_of_steps()) + ' : ' + str(energies) + ' : std(kcal/mol)=' + str(np.std(energies))
    print(output)

f.close()
end_time2 = time.time()
print('CV MD Total Time:', end_time2 - start_time2)
mdcrd.close()
temp.close()
