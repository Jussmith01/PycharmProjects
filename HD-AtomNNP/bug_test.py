#without ase
import os
import pyNeuroChem as pync
import numpy as np
from glob import glob
import hdnntools as hdt # Library for loading conformer data
import time as tm

# ani_energy_OK : Generally, you should only need to load the network once
#                 per run time. Doing it multiple times is unecessarily costly.
#                 However, the memory leak should not be happening, and I suspect
#                 an issue with the python interface not properly calling the
#                 destructor of the C++ classes. I will need to investigate.
#                 further. In any case, DO NOT redeclare your networks as you did
#                 in this function.

wkdir = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile = wkdir + 'sae_6-31gd.dat'
nnfdir = wkdir + 'networks/'

# Delcare NeuroChems molecule class
ncm = pync.molecule(cnstfile, saefile, nnfdir, 0)

# ani_energy_bad : This is not causing any issues for me. My results from this
#                  correlate to the reference DFT calculations with low error
#                  as we expect. Unfortunately I don't have a Quadro to test
#                  the code on. Please double that your inputs into setMolecule
#                  are correctly formatted and ordered.
def ani_energy_bad(xyz,typ):
    # full setup NN for every molecule
    ncm.setMolecule(xyz,list(typ))
    e = ncm.energy()[0]
    return e

def get_eint_m(filename):
    xyz, typ, Eact = hdt.readncdat(filename, np.float32)
    energies = [ani_energy_bad(m,typ) for m in xyz]
    return abs((energies-Eact).sum()/float(len(energies))) * 627.5095

energies1 = []
i = 0
_tb = tm.time()
for f in glob('/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_01/testdata/*.dat'):
    energies1.append(get_eint_m(f))
    i += 1
print('Using molecule class - Time: ' + "{:.4f}".format((tm.time() - _tb) * 1000.0)  + 'ms')

# Instead of using molecule for computing the energies of a bunch of identical conformers,
# we already have a class optimized for doing this. It is the "conformers" class.
# Below I give an example of it in use for a similar problem as yours. I added timers to
# each function call to show you exactly how much time you can save.

# Delcare NeuroChems conformers class
# As with molecule, this loads the NN, so only call it once.
ncc = pync.conformers(cnstfile, saefile, nnfdir, 0)

def ani_energy_conf(xyz,typ):

    # setConformers computes the energies of N conformers for the SAME molecule
    # Arg 1: 3D numpy array of coords (conformers x atoms x coordinates)
    # Arg 2: A single list of the atom types in the correct order.
    ncc.setConformers(xyz,list(typ))

    # returns a numpy array of energies
    e = ncc.energy()
    return e

def get_eint_c(filename):

    # Reads conformers from a file
    xyz, typ, Eact = hdt.readncdat(filename, np.float32)

    energies = ani_energy_conf(xyz,typ)
    return abs((energies-Eact).sum()/float(len(energies))) * 627.5095

energies2 = []
i = 0
_tb = tm.time()
for f in glob('/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_01/testdata/*.dat'):
    energies2.append(get_eint_c(f))
    i += 1
print('Using conformers class - Time: ' + "{:.4f}".format((tm.time() - _tb) * 1000.0)  + 'ms')

#-----------------------CHECK RESULTS----------------------
# These results differ slightly only a few digits in from numerical precision.
# This is expected from our use of single point precision.
print("\nError (kcal/mol) from ref. DFT for nc molecule   class: \n", energies1 )
print("\nError (kcal/mol) from ref. DFT for nc conformers class: \n", energies2 )
