__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
import numpy as np

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Scratch/Research/trainingcases/wB97X-631gd-train-highgarden/train_08-a3.1A_r4.6_dn1/'
cnstfile = wkdir + 'rHCNO-4.6A_32-3.1A_a8-8.params'
saefile  = wkdir + '../sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile,saefile,nnfdir,0)

xyz,typ,Na = gt.readxyz('/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_testdata/butane/opt_test.xyz')

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ[0])

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(nc.getNumConfs()) )

# O of Conformations
dr = 0.0001

F1 = np.array(nc.computeForces(dr=dr))

print ( "{:.7f}".format( dr ) + ': ' + ' VALUES: ' +  "{:.7f}".format(F1[0][0]) )
RMS = np.sqrt(F1**2).sum()/float(F1[0].shape[0])
print (str(RMS))
print ("Forces: ")

N = int(int(F1[0].shape[0])/3)
for i in range(N):
    print (str(i+1) + ' ' + str(typ[0][i]) + ' (' + "{:.7f}".format(F1[0][i*3]) + ',' + "{:.7f}".format(F1[0][i*3+1]) + ',' + "{:.7f}".format(F1[0][i*3+2]) + ')')