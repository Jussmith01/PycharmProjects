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

path = '/home/jujuman/Research/ANI-DATASET/ANI-1_release/ani-1_data_c08_test.h5'

# opening file
adl = pyt.anidataloader(path)

l_E = []

maxi = 0.0
mini = 100000000.0

for i,data in enumerate(adl):
    #print(g.name.split("_s")[1],' : ',data["name"])

    E = data['energies']
    Ea = compute_sae(np.asarray(data['species']).flatten())
    E = E - Ea
    l_E.append(E)

adl.cleanup()
#print('Max: ', maxi,' Min: ',mini)

Ec =  np.concatenate(l_E)
sns.set_style("whitegrid")
ax = sns.distplot(Ec)

ax.set(xlabel='Ea (Ha)', ylabel='Normalized population')

sns.plt.show()