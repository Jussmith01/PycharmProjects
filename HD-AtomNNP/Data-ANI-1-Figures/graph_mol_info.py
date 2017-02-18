__author__ = 'jujuman'

import numpy as np
import seaborn as sns
import pandas as pd
import hdnntools as hdn

def convertatomicnumber(X):
    if X == 'H':
        return 1
    elif X == 'C':
        return 6
    elif X == 'N':
        return 7
    elif X == 'O':
        return 8

def getNumberElectrons(types):
    Ne = 0

    for t in types:
        Ne = Ne + convertatomicnumber(t)

    return Ne

def compute_sae(spec):
    sae = 0.0

    for s in spec:
        if s == 'H':
            sae += -0.500607632585
        elif s == 'C':
            sae += -37.8302333826
        elif s == 'N':
            sae += -54.5680045287
        elif s == 'O':
            sae += -75.0362229210
        else:
            print('Error, unknown type: ', s)
            exit(1)
    return sae

path = '/home/jujuman/Python/PycharmProjects/HD-AtomNNP/Data-ANI-1-Figures/ethane_CC_disso_data.h5'

import pyNeuroChem as pync

# Set required files for pyNeuroChem
#wkdir    = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk/'
wkdir    = '/home/jujuman/Research/wB97X-631gd-train-highgarden/train_08e_disso_test/train_08-a3.1A_r4.6_AEV384_newtrain/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = wkdir + '../../sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

con = pync.conformers(cnstfile, saefile, nnfdir, 0)

# opening file
store = pd.HDFStore(path, complib='blosc',complevel=0)

df_E = []
df_Ea = []
df_D = []

maxi = 0.0
mini = 100000000.0

spec_list = []
coord_list = []

for x in store.get_node(""):
    print("Name:", x._v_name)

    index = np.unique(np.asarray(store.select(x._v_name + "/" + x.energy._v_name).index.get_level_values('molecule')))

    for i in index:
        print('Index: ', i)
        energy = store.select(x._v_name + "/" + x.energy._v_name, 'molecule == ' + str(i))
        species = store.select(x._v_name + "/" + x.species._v_name, 'molecule == ' + str(i))
        coords = store.select(x._v_name + "/" + x.xyz._v_name, 'molecule == ' + str(i))

        xyz = np.array(coords,dtype=np.float32).flatten().reshape(len(energy),len(species),3)

        spec_list.append(list(np.array(species).flatten()))
        coord_list.append(xyz)

        con.setConformers(confs=xyz, types=list(np.array(species).flatten()))

        EANI = con.energy()

        #print(species)
        for c in xyz:
            df_D.append(np.linalg.norm(c[0]-c[1]))

        #Ne = getNumberElectrons(np.asarray(species).flatten())
        sae = compute_sae(np.asarray(species).flatten())

        # Subtract interaction energies
        E = np.asarray(energy).flatten()
        EANI = np.asarray(EANI).flatten()

        # Max
        if abs(E.max()) > maxi:
            maxi = abs(E.max())

        # Min
        if abs(E.min()) < mini:
            mini = abs(E.min())

        #
        df_E.append(E)
        df_Ea.append(EANI)
        #print (df_E)

hdn.writexyzfile("coordlist.xyz",np.concatenate(coord_list),spec_list[0])

import matplotlib.pyplot as plt

df_E = hdn.hatokcal * np.array(df_E).flatten()
df_Ea = hdn.hatokcal * np.array(df_Ea).flatten()

rmse = hdn.calculaterootmeansqrerror(df_E,df_Ea)

x = list(range(0,df_E.shape[0]))

plt.scatter(df_D, df_E,label='DFT')
plt.scatter(df_D, df_Ea,color='red',label='ANI RMSE: '+ "{:.2f}".format(rmse) + 'kcal/mol',alpha=0.5)
plt.xlabel('Distance ($\AA$)')
plt.ylabel('Energy (kcal/mol)')
plt.legend(bbox_to_anchor=(0.4, 0.99), loc=2, borderaxespad=0.,fontsize=14)

plt.show()

plt.scatter(x, df_D)
plt.xlabel('Step')
plt.ylabel('Distance ($\AA$)')

plt.show()


store.close()