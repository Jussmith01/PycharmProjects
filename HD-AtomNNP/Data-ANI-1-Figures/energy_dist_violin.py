__author__ = 'jujuman'

import numpy as np
import seaborn as sns
import pandas as pd

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

#path = "/home/jujuman/Toshiba460GB/Research/ANI-1-DATA-PAPER-FILES/data/data-ani-1.h5"
path = '/home/jujuman/Python/PycharmProjects/HD-AtomNNP/Data-ANI-1-Figures/ethane_CC_disso_data.h5'
# opening file
store = pd.HDFStore(path, complib='blosc',complevel=0)
print(store)

df_E = []

maxi = 0.0
mini = 100000000.0

for x in store.get_node(""):
    print("Name:", x._v_name)

    index = np.unique(np.asarray(store.select(x._v_name + "/" + x.energy._v_name).index.get_level_values('molecule')))

    for i in index:
        print('Index: ', i)
        energy = store.select(x._v_name + "/" + x.energy._v_name, 'molecule == ' + str(i))
        species = store.select(x._v_name + "/" + x.species._v_name, 'molecule == ' + str(i))
        #Ne = getNumberElectrons(np.asarray(species).flatten())
        sae = compute_sae(np.asarray(species).flatten())

        # Subtract interaction energies
        E = np.asarray(energy).flatten() - sae

        # Max
        if abs(E.max()) > maxi:
            maxi = abs(E.max())

        # Min
        if abs(E.min()) < mini:
            mini = abs(E.min())

        #
        df_E.append(pd.DataFrame(E, columns=[x._v_name]))

        #print (df_E)

store.close()

print('Max: ', maxi,' Min: ',mini)

df_Ec =  pd.concat(df_E)
ax = sns.violinplot(x=df_Ec,scale="count")

sns.plt.show()