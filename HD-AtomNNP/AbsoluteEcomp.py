__author__ = 'jujuman'

import numpy as np
import statsmodels.api as sm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
import random
import graphtools as gt


# -----------------------
cmap = mpl.cm.brg
# ------------5412.mordor
# AM1 vs Act
# ------------
user = os.environ['USER']

#dir = '/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/train_06AA/'
#dir1 = '/Research/GDB-11-wB97X-6-31gd/train_07/'
#dir1 = '/Research/trainingcases/wB97X-631gd-train-comet/train_07_a2.9A_r5.2A/'
dir1 = '/Research/trainingcases/wB97X-631gd-train-highgarden/train_08-a3.5A_r5.0_dn2/'
#dir2 = '/Research/trainingcases/wB97X-631gd-train-highgarden/nw-64-64-64-32/train_07-a3.1A_r4.5/'

#file = 'pp_01_test.dat_graph'
#file = 'polypep_test.dat_graph'
#file2 = 'polypepPM6_test.dat_graph'
#file = 'pentadecane_test.dat_graph'
#file2 = 'pentadecanePM6_test.dat_graph'
#file = 'C8H16-isomers_test.dat_graph'
#file2 = 'C8H16-isomersPM6_test.dat_graph'
#file = 'benzamide_conformers-0_test.dat_graph'
#file2 = 'benzamide_conformersPM6-0_test.dat_graph'
#file = 'retinolconformer_test.dat_graph'
#file2 = 'retinolconformerPM6_test.dat_graph'
#file='bdpd_test.dat_graph'
#file2='bdpdPM6_test.dat_graph'
#file='ethenedimer_test.dat_graph'
#file2='ethenedimerPM6_test.dat_graph'
#file='formicaciddimer_test.dat_graph'
#file2='formicaciddimerPM6_test.dat_graph'
#file='waterdimer_test.dat_graph'
#file2='waterdimerPM6_test.dat_graph'
#file='ammoniadimer_test.dat_graph'
#file2='ammoniadimerPM6_test.dat_graph'
#file = 'atazanavir_AM1_CLN_test.dat_graph'
file ='h2o2dhl_test.dat_graph'
#file1='c2h2disdata_train.dat_graph'
#file2='h2bondscanR_test.dat_graph'

data1 = gt.getfltsfromfile('/home/' + user + dir1 + file,' ', [0])
data2 = gt.convert * gt.getfltsfromfile('/home/' + user + dir1 + file,' ', [1])
data3 = gt.convert * gt.getfltsfromfile('/home/' + user + dir1 + file,' ', [2])
#data4 = gt.convert * gt.getfltsfromfile('/home/' + user + dir1 + file2,' ', [1])

#data4 = gt.getfltsfromfile('/home/' + user + dir2 + file,' ', [1])
#data5 = gt.getfltsfromfile('/home/' + user + dir2 + file,' ', [2])
#data5 = gt.getfltsfromfile('/home/' + user + dir3 + file, [2])

#mean = np.mean(data3)

#data2 = data2 - np.mean(data2)
#data3 = data3 - np.mean(data3)
#data4 = data4 - np.mean(data4)

rmse1 = gt.calculaterootmeansqrerror(data2,data3)
#rmse2 = gt.calculaterootmeansqrerror(data2,data4)
#rmse3 = gt.calculaterootmeansqrerror(data2,data5)

print('Datasize: ' + str(data1.shape[0]))

font = {'family' : 'Bitstream Vera Sans',
            'weight' : 'normal',
            'size'   : 14}

plt.rc('font', **font)

data1 = data1

#plt.plot(data1, data2, color='black', label='wB97X/6-31G*',linewidth=4)
plt.scatter(data1, data2, color='black',linewidth=3)
plt.scatter(data1, data3, color='blue',linewidth=3)
#plt.plot(data1, data3, '--', color='blue', label='ANN - c08b RMSE: ' + "{:.6f}".format(rmse1) + "kcal/mol",linewidth=2)
#plt.scatter(data1, data4, color='red',linewidth=3)
#plt.plot(data1, data4, '--', color='red', label='ANN - c07b RMSE: ' + "{:.6f}".format(rmse2) + "kcal/mol",linewidth=2)

#plt.plot(data4, data4, color='black',linewidth=2)
#plt.scatter(data4, data5, color='green', label='ANN - GDB-62 RMSE: ' + "{:.6f}".format(rmse2) + "Ha",linewidth=4)
#plt.scatter(data1, data5, color='blue', label='ANN - GDB-61 RMSE: ' + "{:.6f}".format(rmse3) + "Ha",linewidth=4)

#plt.title('300K norm. mode generated conformers of\nH-Gly-Pro-Hyp-Gly-Ala-Gly-OH')
#plt.title("Absolute energy (shifted by mean) of Retinol conformers")
#plt.title("Absolute energy (shifted by mean) of Benzamide\nStructures generated at 300K")
plt.title("Water H-O bond scan")

#plt.xlabel('Conformation Pair (Count 49)')
plt.ylabel('E (kcal/mol)')
#plt.ylabel('Calculated Absolute E (kcal/mol)')
plt.xlabel('Coordinate ($\AA$)')
plt.legend(bbox_to_anchor=(0.3, 0.975), loc=2, borderaxespad=0.,fontsize=12)


# -----
# PLOT
# -----
plt.show()
