__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-1-ntwk/'
cnstfile = wkdir + 'rHCNO-4.6A_32-3.1A_a8-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'


wkdir    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/train_08_5/'
cnstfile = wkdir + 'rHCNO-4.7A_32-3.2A_a8-8.params'
saefile  = wkdir + '../sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

xyz = [[],[],[]]
typ = [[],[],[]]
Na  = [[],[],[]]

xyz[0],typ[0],Na[0] = gt.readxyz('/home/jujuman/Scratch/Research/TestingCases/DimerSeparationTests/dimerstruct.xyz')
xyz[1],typ[1],Na[1] = gt.readxyz('/home/jujuman/Scratch/Research/TestingCases/DimerSeparationTests/monomer1.xyz')
xyz[2],typ[2],Na[2] = gt.readxyz('/home/jujuman/Scratch/Research/TestingCases/DimerSeparationTests/monomer2.xyz')

Ea1 = -496.9336409
Ea2 = -248.4579414
Ea3 = -248.4577671

E = []

for i in range(len(xyz)):
	# Construct pyNeuroChem class
	nc = pync.pyNeuroChem(cnstfile, saefile, nnfdir, 0)

	# Set the conformers in NeuroChem
	nc.setConformers(confs=xyz[i],types=typ[i][0])

	# Print some data from the NeuroChem
	print( 'Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
	print( 'Number of Confs Loaded: ' + str(nc.getNumConfs()) )

	# Compute Energies of Conformations
	E.append(nc.computeEnergies()[0])

print(E)

print ('-----------------DATA---------------')
d1 = gt.hatokcal*(E[0] - (E[1] + E[2]))
d2 = gt.hatokcal*(Ea1 - (Ea2 + Ea3))
print (d1)
print (d2)
print (d1 - d2)