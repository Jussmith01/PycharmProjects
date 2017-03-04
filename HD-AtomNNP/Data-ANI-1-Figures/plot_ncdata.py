import numpy as np
import pyanitools as pyt
import hdnntools as hdn
import pyNeuroChem as pync
import matplotlib.pyplot as plt

file = "/home/jujuman/Research/ANI-DATASET/rxn_db_mig.h5"

al = pyt.anidataloader(file)

al.totalload()
data = al.get_all_data()
al.cleanup()

df_E = hdn.hatokcal * data[1].flatten()
xyz = data[0].reshape(df_E.shape[0],len(data[2][0]),3)

hdn.writexyzfile('/home/jujuman/crds.xyz',xyz,data[2][0])

# Set required files for pyNeuroChem
#wkdir    = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk/'
wkdir    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/ANI-c08e-ntwk_newtrain/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

# Construct pyNeuroChem class
mol = pync.conformers(cnstfile, saefile, nnfdir, 0)

mol.setConformers(confs=xyz,types=list(data[2][0]))

E = hdn.hatokcal * mol.energy()

rmse = hdn.calculaterootmeansqrerror(df_E,E)

x = list(range(0,df_E.shape[0]))
#x = np.linalg.norm(xyz[:,3,:]-xyz[0,3,:],axis=1)

#print(x)

plt.scatter(x, df_E,label='DFT')
plt.scatter(x, E   ,label='ANI err: '+str(rmse)+' kcal/mol')
plt.xlabel('Distance ($\AA$)')
plt.ylabel('Energy (kcal/mol)')
plt.legend(bbox_to_anchor=(0.4, 0.99), loc=2, borderaxespad=0.,fontsize=14)

plt.show()
