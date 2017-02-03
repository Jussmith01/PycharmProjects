__author__ = 'jujuman'

import numpy as np
import seaborn as sns
import pandas as pd

print('Loading data...')
enrg1 = np.load('data/energy/GDB-01_enrg.npz')['arr_0']
enrg2 = np.load('data/energy/GDB-02_enrg.npz')['arr_0']

df_E = pd.DataFrame(enrg1, np.fill , columns=["e1","e2"])

print (df_E)

#ax = sns.violinplot(x=enrg4,)

#sns.plt.show()