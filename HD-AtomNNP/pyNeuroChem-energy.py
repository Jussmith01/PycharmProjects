__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Research/trainingcases/wB97X-631gd-train-highgarden/'
cnstfile = wkdir + 'train_08-a3.5A_r5.0_dn4/rHCNO-4.5A_32-3.1A_a8-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'train_08-a3.5A_r5.0_dn4/networks/'

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile,saefile,nnfdir,1)

xyz,typ,Na = gt.readxyz('/home/jujuman/python/PycharmProjects/HD-AtomNNP/dipeptide.xyz')

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ[0])

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(nc.getNumConfs()) )

# Compute Forces of Conformations
E1 = nc.computeEnergies()
E2 = nc.computeEnergies()

print ('-----------------ORIGINAL COORDS---------------')
for i in E1:
	print ("{:.30f}".format(i))

for i in E2:
	print ("{:.30f}".format(i))
