import pyanitools as pya
import hdnntools as hdt
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Set the HDF5 file containing the data
hdf5file = '/home/jujuman/Research/ANI-DATASET/ani_data_mp_param.h5'

# Construct the data loader class
adl = pya.anidataloader(hdf5file)

xt_data = []
yt_data = []

xv_data = []
yv_data = []

# Print the species of the data set one by one
for data in adl.getnextdata():

    # Extract the data
    xyz = np.array_split(data['coordinates'], 250)
    erg = np.array_split(data['energies'], 250)
    spc = data['species']

    sae = hdt.compute_sae('/home/jujuman/Research/GDB-11-wB97X-6-31gd/sae_6-31gd.dat', spc)
    print('TShape: ', xyz[0].shape[0],' VShape: ', xyz[1].shape[0])


    xt_data.append(hdt.ncdata(hdt.generatedmatsd3(xyz[0]), spc, xyz[0].shape[1]))
    yt_data.append(erg[0] - sae)

    xv_data.append(hdt.ncdata(hdt.generatedmatsd3(xyz[1]), spc, xyz[1].shape[1]))
    yv_data.append(erg[1] - sae)

yt_data = np.concatenate(yt_data)
yv_data = np.concatenate(yv_data)

print('Training Data Shape: ', yt_data.shape)
print('Testing  Data Shape: ', yv_data.shape)

p0 = (1.1, -0.1,# HH
      1.1, -0.1,# HC
      1.1, -0.1,# HN
      1.1, -0.1,# HO
      1.1, -0.1,# CC
      1.1, -0.1,# CN
      1.1, -0.1,# CO
      1.1, -0.1,# OO
      1.1, -0.1,# ON
      1.1, -0.1,)# NN

bounds = ([0.75, -2.0,
           0.75, -2.0,
           0.75, -2.0,
           0.75, -2.0,
           0.75, -2.0,
           0.75, -2.0,
           0.75, -2.0,
           0.75, -2.0,
           0.75, -2.0,
           0.75, -2.0,],
          [2.0, 2.0,
           2.0, 2.0,
           2.0, 2.0,
           2.0, 2.0,
           2.0, 2.0,
           2.0, 2.0,
           2.0, 2.0,
           2.0, 2.0,
           2.0, 2.0,
           2.0, 2.0,])


popt, pcov = curve_fit(hdt.buckingham_pot, xt_data, yt_data, p0=p0, bounds=bounds)# NN

print(popt)
iEc = hdt.buckingham_pot(xv_data, *p0)
fEc = hdt.buckingham_pot(xv_data, *popt)

irmse = hdt.calculaterootmeansqrerror(iEc, yv_data)
frmse = hdt.calculaterootmeansqrerror(fEc, yv_data)

np.savez('mp_ani_params_gdb06.npz', param=popt)

print('Final RMSE:', hdt.hatokcal * frmse, ' Initial RMSE:', hdt.hatokcal * irmse)

plt.plot(yv_data, yv_data, color='black', label='Act')
plt.scatter(yv_data, iEc, color='red', label='Init')
plt.scatter(yv_data, fEc, color='blue', label='Fit')

plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()

# Closes the H5 data file
adl.cleanup()
