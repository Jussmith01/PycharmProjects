__author__ = 'jujuman'

import numpy as np
import pandas as pd
import seaborn as sns
import pyanitools as pyt
import hdnntools as hdt

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

path = '/home/jujuman/Research/ANI-DATASET/ANI-1_release/ani-1_data_c08.h5'

# opening file
adl = pyt.anidataloader(path)

gl = adl.get_group_list()[3:]

print(gl)

df_E = []

maxi = 0.0
mini = 100000000.0

for g in gl:
    print(g.name)

    for i,data in enumerate(adl.iter_group(g)):
        #print(g.name.split("_s")[1],' : ',data["name"])

        E = data['energies']
        Ne = getNumberElectrons(np.asarray(data['species']).flatten())
        E = E/float(Ne)
        df_E.append(pd.DataFrame(E, columns=[str(g.name.split("_s")[1])]))

    '''
    for i in index:
        print('Index: ', i)
        energy = store.select(x._v_name + "/" + x.energy._v_name, 'molecule == ' + str(i))
        species = store.select(x._v_name + "/" + x.species._v_name, 'molecule == ' + str(i))

        #sae = compute_sae(np.asarray(species).flatten())

        # Subtract interaction energies
        E = np.asarray(energy).flatten()/float(Ne)

        # Max
        if abs(E.max()) > maxi:
            maxi = abs(E.max())

        # Min
        if abs(E.min()) < mini:
            mini = abs(E.min())

        #

        #print (df_E)
    '''

adl.cleanup()
#print('Max: ', maxi,' Min: ',mini)

df_Ec =  pd.concat(df_E)
sns.set_style("whitegrid")
ax = sns.violinplot(x=df_Ec, scale="area",bw=0.2)

ax.set(xlabel='GDB subset', ylabel='Totale energy / number of electrons (Ha)')

sns.plt.show()