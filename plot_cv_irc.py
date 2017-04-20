# Import pyNeuroChem
import pyNeuroChem as pync
import numpy as np
import hdnntools as hdt

# Set required files for pyNeuroChem
anipath  = '/home/jujuman/Dropbox/ChemSciencePaper.AER/networks/ANI-c08f09dd-ntwk-cv/cv_c08e_ntw_0'
cnstfile = anipath + '/../rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = anipath + '/../sae_6-31gd.dat'
nnfdir   = anipath + '/networks/'

file = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_rxns/benz_dbm_rxns/scans_dbm_benz_8/irc.dat'
xyz,typ,Eact = hdt.readncdat(file,np.float32)

# Construct pyNeuroChem class
nc = pync.conformers(cnstfile, saefile, nnfdir, 0)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=list(typ))

# Compute Energies of Conformations
E1 = nc.energy()

print('Energy: ', E1 - Eact)

import matplotlib.pyplot as plt

plt.plot (hdt.hatokcal*(E1), color='blue',  label='ANI',  linewidth=2)
plt.plot (hdt.hatokcal*(Eact), color='black',  label='DFT',  linewidth=2)

#plt.plot (x, morse(x,1.0,1.9,1.53), color='grey',  label='SRC+ANI',  linewidth=2)

plt.title("C-C Dissociation")

plt.ylabel('E (kcal/mol)')
plt.xlabel('Distance $\AA$')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.,fontsize=16)

plt.show()