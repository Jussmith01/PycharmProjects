__author__ = 'jujuman'

import numpy as np
import hdnntools as hdt
import matplotlib.pyplot as plt

data = np.loadtxt('/home/jujuman/Research/ANI-DATASET/ANI-SF-TRAIN/cv_c08e_ntw_0/cost.dat',delimiter=' ',dtype=np.float32)
print(data)
plt.plot (data[1:,0],hdt.hatokcal*np.sqrt(data[1:,1]),marker='o', color='blue',  label='Train',  linewidth=2)
plt.plot (data[1:,0],hdt.hatokcal*np.sqrt(data[1:,2]),marker='o', color='red',  label='Valid',  linewidth=2)
plt.plot (data[1:,0],hdt.hatokcal*np.sqrt(data[1:,3]),marker='o', color='green',  label='Best',  linewidth=2)

plt.yscale('log')

plt.title("C10H20 - ANI vs DFT")

plt.ylabel('Error')
plt.xlabel('Epoch')
plt.legend(bbox_to_anchor=(0.8, 0.95), loc=2, borderaxespad=0.,fontsize=16)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

plt.show()