__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/train_08_3/'
cnstfile = wkdir + 'rHCNO-4.6A_32-3.1A_a8-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile,saefile,nnfdir,0)

xyz,typ,Na = gt.readxyz('/home/jujuman/Dropbox/Research/AceticAcidDimerProtonTransfer/TS.xyz')

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ[0])

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(nc.getNumConfs()) )

# O of Conformations
F = nc.optimizeGeom(0,conv=0.00001,step=0.001,dr=0.0009)

print ('-----------------ORIGINAL COORDS---------------')
for i in xyz:
    for j in range(0,Na[0]):
        print (typ[0][j] + ' ' + "{:.7f}".format(i[j*3+0]) + ' ' + "{:.7f}".format(i[j*3+1]) + ' ' + "{:.7f}".format(i[j*3+2]))

print ('-----------------OPTIMIZED COORDS---------------')
for i in F:
    for j in range(0,Na[0]):
        print (typ[0][j] + ' ' + "{:.7f}".format(i[j*3+0]) + ' ' + "{:.7f}".format(i[j*3+1]) + ' ' + "{:.7f}".format(i[j*3+2]))
