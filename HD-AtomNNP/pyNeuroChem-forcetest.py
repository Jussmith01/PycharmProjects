__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
import numpy as np

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Research/NeuroChemForceTesting/train_01/'
cnstfile = wkdir + 'rH-3.0A_4-2.5A_a2-2.params'
saefile  = wkdir + '../sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

xyz = [[0.00000,0.00000,0.37124,0.00000, 0.00000,-0.37124]]
	  #,[0.00000,0.36000,0.00000,0.00000,-0.36000, 0.00000]]
	  #,[0.00000,0.00000,0.41000,0.00000, 0.00000,-0.41000]]

typ = ['H','H']

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile,saefile,nnfdir,0)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ)

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
print F1

N = int(int(F1[0].shape[0])/3)
for i in range(N-1):
    print (str(i+1) + ' ' + str(typ[0][i]) + ' (' + "{:.7f}".format(F1[0][i*3]) + ',' + "{:.7f}".format(F1[0][i*3+1]) + ',' + "{:.7f}".format(F1[0][i*3+2]) + ')')