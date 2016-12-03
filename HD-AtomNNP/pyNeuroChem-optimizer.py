__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import time as tm
import numpy as np
import graphtools as gt

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-1-ntwk/'
cnstfile = wkdir + 'rHCNO-4.6A_32-3.1A_a8-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile,saefile,nnfdir,0)

xyz,typ,Na = gt.readxyz('/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_begdb/begdb-h2oclusters/xyz/4216_water10PP2.xyz')
xyz1,typ1,Na1 = gt.readxyz('/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_begdb/begdb-h2oclusters/xyz/4216_water10PP2.xyz')

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ[0])

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(nc.getNumConfs()) )

# O of Conformations
print('Optimizing...')
_t1b = tm.time()
F = nc.optimizeGeomBFGS(conv=0.000001)
print('Optimization complete. Time: ' + "{:.4f}".format((tm.time() - _t1b) * 1000.0) + 'ms')

print ('-----------------ORIGINAL COORDS---------------')
for i in xyz:
    for j in range(0,Na[0]):
        print (typ[0][j] + ' ' + "{:.7f}".format(i[j*3+0]) + ' ' + "{:.7f}".format(i[j*3+1]) + ' ' + "{:.7f}".format(i[j*3+2]))

print ('-----------------OPTIMIZED COORDS---------------')
for i in F:
    for j in range(0,Na[0]):
        print (typ[0][j] + ' ' + "{:.7f}".format(i[j*3+0]) + ' ' + "{:.7f}".format(i[j*3+1]) + ' ' + "{:.7f}".format(i[j*3+2]))


d0 = np.array(gt.generatedmat(xyz[0],Na[0]))
print("Initial (distance matrix):")
print(d0)

d1 = np.array(gt.generatedmat(xyz1[0],Na[0]))
print("wb97x opt (distance matrix):")
print(d1)

d2 = np.array(gt.generatedmat(F[0],Na[0]))
print("ANI opt (distance matrix):")
print(d2)

print("\nRMSE of distance matrix --")
print("   ANI opt   vs. wb97x opt: " + "{:.7f}".format(gt.calculaterootmeansqrerror(d1,d2)))
print("   wb97x opt vs. pre-opt: " + "{:.7f}".format(gt.calculaterootmeansqrerror(d1,d0)))
print("   ANI opt   vs. pre-opt  : " + "{:.7f}".format(gt.calculaterootmeansqrerror(d2,d0)))
