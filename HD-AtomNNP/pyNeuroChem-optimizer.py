__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import time as tm
import numpy as np
import graphtools as gt

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Dropbox/Research/NeuroChemForceTesting/train_02/'
cnstfile = wkdir + 'rHO-3.0A_4-2.5A_a2-2.params'
saefile  = wkdir + '../sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile,saefile,nnfdir,0)

#xyz,typ,Na = gt.readxyz('/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_testdata/butane/opt_test.xyz')

xyz = [[0.000000,    0.000000,    0.198604, 0.000000,    0.95819,   -0.424415, 0.000000,-0.759819,-0.554415]]#, [0.0, 0.0, 0.118604, 0.0, 0.67086124, -0.40498585, 0.0, -0.75981897, -0.474415]]
xyz1 = [[0.000000,    0.000000,    0.118604, 0.000000,    0.759819,   -0.474415, 0.000000,-0.759819,-0.474415]]#, [0.0, 0.0, 0.118604, 0.0, 0.67086124, -0.40498585, 0.0, -0.75981897, -0.474415]]
xyz2 = [0.000000,    0.000000,    0.119814,0.000000,    0.761299,   -0.479258,0.000000,   -0.761299,   -0.479258]

Na = 3
typ = ['O','H','H']


# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ)

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(nc.getNumConfs()) )

# O of Conformations
print('Optimizing...')
_t1b = tm.time()
F = nc.optimizeGeom(0,conv=0.0000005,step=0.1,dr=0.0066514)
print('Optimization complete. Time: ' + "{:.4f}".format((tm.time() - _t1b) * 1000.0) + 'ms')

print ('-----------------ORIGINAL COORDS---------------')
for i in xyz1:
    for j in range(0,Na):
        print (typ[j] + ' ' + "{:.7f}".format(i[j*3+0]) + ' ' + "{:.7f}".format(i[j*3+1]) + ' ' + "{:.7f}".format(i[j*3+2]))

print ('-----------------OPTIMIZED COORDS---------------')
for i in F:
    for j in range(0,Na):
        print (typ[j] + ' ' + "{:.7f}".format(i[j*3+0]) + ' ' + "{:.7f}".format(i[j*3+1]) + ' ' + "{:.7f}".format(i[j*3+2]))

d0 = np.array(gt.generatedmat(xyz[0],Na))
d1 = np.array(gt.generatedmat(xyz1[0],Na))
print(d1)
d2 = np.array(gt.generatedmat(F[0],Na))
print(d2)
d3 = np.array(gt.generatedmat(xyz2,Na))

print("\nRMSE of distance matrix --")
print("   pre-opt   vs. wb97x opt: " + "{:.7f}".format(gt.calculaterootmeansqrerror(d1,d0)))
print("   ANI opt   vs. wb97x opt: " + "{:.7f}".format(gt.calculaterootmeansqrerror(d1,d2)))
print("   b3lyp opt vs. wb97x opt: " + "{:.7f}".format(gt.calculaterootmeansqrerror(d1,d3)))
