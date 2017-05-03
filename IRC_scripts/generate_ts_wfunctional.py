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

import math3d as m3d

#import matplotlib
#import matplotlib as mpl

#import matplotlib.pyplot as plt

#import seaborn as sns
#%matplotlib inline

# Set required files for pyNeuroChem
anipath  = '/home/jujuman/Dropbox/ChemSciencePaper.AER/networks/ANI-c08f09dd-ntwk-cv/cv_c08e_ntw_0'
cnstfile = anipath + '/../rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = anipath + '/../sae_6-31gd.dat'
nnfdir   = anipath + '/networks/'

reslist = [0,1,2,3,4,5,9]

# Construct pyNeuroChem class
nc = pync.molecule(cnstfile, saefile, nnfdir, 1)

#bz = read('C_100.xyz')
bz = read('/home/jujuman/Research/ReactionGeneration/ts.xyz')

print(bz.get_positions())
print(bz.get_chemical_symbols())
spc = bz.get_chemical_symbols()
spc[11] = 'C'
pos = bz.get_positions()
v1 = m3d.Vector(pos[2])
v2 = m3d.Vector(pos[11])

ve = m3d.Vector(pos[0])

v3 = 1.45 * (v2-v1).get_normalized()
print(v3.get_length())
print(v1 + v3)

print('Cross: ',v3.cross(ve))

vc = (v1 + v3)
pos[11] = vc.get_array_ref()

o = m3d.Orientation()
cv = v3.cross(ve)
o.rotate_b(axis=cv, angle=(180-109.5) * 0.0174533)

vH1 = 1.0*(o * v3) + vc
spc.append('H')
pos = np.vstack([pos,vH1.get_array_ref()])

oc = m3d.Orientation()
oc.rotate_b(axis=v3,angle=2*np.pi - 120.0 * 0.0174533)
cv = oc * cv

o = m3d.Orientation()
o.rotate_b(axis=cv,angle=(180-109.5) * 0.0174533)

vH2 = 1.0*(o * v3) + vc
spc.append('H')
pos = np.vstack([pos,vH2.get_array_ref()])

oc = m3d.Orientation()
oc.rotate_b(axis=v3,angle=2*np.pi -120.0 * 0.0174533)
cv = oc * cv

o = m3d.Orientation()
o.rotate_b(axis=cv,angle=(180-109.5) * 0.0174533)

vH2 = 1.0*(o * v3) + vc
spc.append('H')
pos = np.vstack([pos,vH2.get_array_ref()])

bz = ase.Atoms(positions=pos,symbols=spc)
print(bz.get_chemical_symbols())

f = open("/home/jujuman/Research/ReactionGeneration/optmol_begin.xyz",'w')
f.write(str(len(bz)) + '\n   comment\n')
for i in bz:
    f.write(str(i.symbol) + ' ' + str(i.x) + ' ' + str(i.y) + ' ' + str(i.z) + '\n')
f.close()

bz.set_calculator(ANI(False,1,reslist))
bz.calc.setnc(nc)
#
start_time = time.time()
dyn = LBFGS(bz)
dyn.run(fmax=0.0001)
print('[ANI Total time:', time.time() - start_time, 'seconds]')

# Write visualization of molecule
f = open("/home/jujuman/Research/ReactionGeneration/optmol_end.xyz",'w')
f.write(str(len(bz)) + '\n   comment\n')
for i in bz:
    f.write(str(i.symbol) + ' ' + str(i.x) + ' ' + str(i.y) + ' ' + str(i.z) + '\n')
f.close()
