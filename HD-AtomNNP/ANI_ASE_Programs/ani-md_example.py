import sys
import time

# Neuro Chem
from ase_interface import ANI

import  ase
from ase.md.langevin import Langevin
from ase import units

#from ase.neb import NEBtools
from ase.io import read, write
from ase.optimize import LBFGS

# Read molecule
bz = read('/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/Fentanyl/Fentanyl_conformers_DFT.xyz')

bz.set_calculator(ANI(gpuid=1))

# Optimize structure
start_time = time.time()
dyn = LBFGS(bz)
dyn.run(fmax=0.0001)
print('[ANI Total time:', time.time() - start_time, 'seconds]')

# MD temperature
T = 300.0

# We want to run MD with constant energy using the Langevin algorithm
# with a time step of 0.25 fs, the temperature T and the friction
# coefficient to 0.01 atomic units.
dyn = Langevin(bz, 0.25 * units.fs, T * units.kB, 0.01)

# Open output files
mdcrd = open("mdcrd.xyz",'w')
temp =  open("temp.dat",'w')

# Define printer function
def printenergy(a=bz,b=mdcrd,d=dyn,t=temp):  # store a reference to atoms in the
    """Function to print the potential, kinetic and total energy."""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    T = ekin / (1.5 * units.kB)
    print('Step %i - Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (d.get_number_of_steps(),epot, ekin, T, epot + ekin))
    t.write(str(d.get_number_of_steps()) + ' ' + str(d.get_time()) + ' ' + str(ekin / (1.5 * units.kB)) + ' ' + str(epot) + ' ' +  str(ekin) + ' ' + str(epot + ekin) + '\n')

    b.write(str(len(a)) + '\n       Temp: ' + str(T) + ' Step: ' + str(d.get_number_of_steps()) + '\n')
    c=a.get_positions(wrap=True)
    for j,i in zip(a,c):
        b.write(str(j.symbol) + ' ' + str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) + '\n')

# Attach printer
dyn.attach(printenergy, interval=10)

# Run dynamics for 20000 steps
start_time2 = time.time()
dyn.run(20000)  # Do 5 ps of MD
end_time2 = time.time()
print('Total Time:', end_time2 - start_time2)

# Close output files
mdcrd.close()
temp.close()
