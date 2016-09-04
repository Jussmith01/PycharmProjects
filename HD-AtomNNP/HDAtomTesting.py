__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Research/trainingcases/wB97X-631gd-train-highgarden/'
cnstfile = wkdir + 'nw-64-64-64-32/train_07-a3.1A_r4.5/rHCNO-4.5A_32-3.1A_a8-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'nw-64-64-64-32/train_07-a3.1A_r4.5/networks/'

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile,saefile,nnfdir,1)

# List of atom types of Molecule (H2O)
typ = ['C','H','H','H','H']

# Save list 4 conformers of H2O
xyz = [[0.000000, 0.000000, 0.00000, 0.6418, 0.6418, 0.6418,-0.6418,-0.6418,0.6418,-0.6418,0.6418,-0.6418,0.6418,-0.6418,-0.6418]]

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ)

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(nc.getNumConfs()) )

# Compute Forces of Conformations
F = nc.optimizeGeom(0)

for i in F:
    for j in range(0,5):
        print (typ[j] + ' ' + "{:.7f}".format(i[j*3+0]) + ' ' + "{:.7f}".format(i[j*3+1]) + ' ' + "{:.7f}".format(i[j*3+2]))