import pygau09tools as g09
import hdnntools as hdn
import numpy as np
import matplotlib.pyplot as plt
import pyNeuroChem as pync

dir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_rxns/scans_double_bond_migration/'

filef = dir + 'IRC_fwd.log'
fileb = dir + 'IRC_bck.log'

dataf = g09.read_irc(filef)
datab = g09.read_irc(fileb)

xyz = np.concatenate([np.flipud(datab[1]),dataf[1]])
eng = np.concatenate([np.flipud(datab[0]),dataf[0]])

hdn.writexyzfile(dir + 'scan.xyz',xyz,dataf[2])

# Set required files for pyNeuroChem
anipath  = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk'
cnstfile = anipath + '/rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = anipath + '/sae_6-31gd.dat'
nnfdir   = anipath + '/networks/'

# Construct pyNeuroChem class
nc = pync.conformers(cnstfile, saefile, nnfdir, 0)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=dataf[2])

ani = nc.energy()

xv = xyz[:,2,:]-xyz[0,2,:]
x = np.linalg.norm(xv, axis=1)

plt.plot(x,hdn.hatokcal*(eng - eng.max()),'r--', color='blue', label='DFT',linewidth=2)
plt.plot(x,hdn.hatokcal*(ani - ani.max()),'r--', color='red', label='ANN',linewidth=2)

plt.title("Energy Differences Between 8 Random Conformations of a Peptide")
#plt.xlabel('Conformation Pair (Count 49)')
plt.xlabel('Theoretical dE (Hartree)')
plt.ylabel('Calculated dE (Hartree)')
plt.legend(bbox_to_anchor=(0.8, 0.975), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()
