# Import pyNeuroChem
import pyNeuroChem as pync
import numpy as np
import hdnntools as hdt
import matplotlib.pyplot as plt

file = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_rxns/benz_dbm_rxns/scans_dbm_benz_1/irc.dat'
xyz,typ,Eact = hdt.readncdat(file,np.float32)

# Set required files for pyNeuroChem
wkdir  = '/home/jujuman/Dropbox/ChemSciencePaper.AER/networks/ANI-c08f-ntwk-cv/'
cnstfile = 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = 'sae_6-31gd.dat'

nc =  [pync.conformers(wkdir + cnstfile, wkdir + saefile, wkdir + 'cv_c08e_ntw_' + str(l) + '/networks/', 0) for l in range(5)]

rcdir  = '/home/jujuman/Research/ANI-DATASET/RXN1_TNET/ani_benz_rxn1_ntwk/'
ncr1 = pync.conformers(rcdir + '../' + cnstfile, rcdir + '../'  + saefile, rcdir + '/networks/', 0)

# Set the conformers in NeuroChem
ncr1.setConformers(confs=xyz, types=list(typ))

# Compute Energies of Conformations
E1 = ncr1.energy()

# Plot
plt.plot(hdt.hatokcal * (E1), color='red', label='ANI', linewidth=2)

plt.plot (hdt.hatokcal*(Eact), color='black',  label='DFT',  linewidth=2)

for net in nc:
    # Set the conformers in NeuroChem
    net.setConformers(confs=xyz,types=list(typ))

    # Compute Energies of Conformations
    E1 = net.energy()

    # Plot
    plt.plot (hdt.hatokcal*(E1), color='blue',  label='ANI',  linewidth=2)

#plt.plot (x, morse(x,1.0,1.9,1.53), color='grey',  label='SRC+ANI',  linewidth=2)

plt.title("C-C Dissociation")

plt.ylabel('E (kcal/mol)')
plt.xlabel('Distance $\AA$')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.,fontsize=16)

plt.show()