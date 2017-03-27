import pyanitools as pya
import hdnntools as hdt
import numpy as np
import matplotlib.pyplot as plt

scan = hdt.readncdat('/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_testdata/waterdimerscan/trainingData.dat', type=np.float32)

sdat = [hdt.ncdata(hdt.generatedmatsd3(scan[0]), scan[1], scan[0].shape[1])]

sae = hdt.compute_sae('/home/jujuman/Research/GDB-11-wB97X-6-31gd/sae_6-31gd.dat', scan[1])
serg = scan[2] - sae

x = 0.1 * np.array(range(serg.shape[0]), dtype=np.float64) + 0.6
print(len(x))

popt = np.load('mp_ani_params_gdb06.npz')['param']
fsEc = hdt.buckingham_pot(sdat, *popt)

frmse = hdt.calculaterootmeansqrerror(serg, fsEc)

plt.plot(x, serg, color='black', label='Act')
plt.scatter(x, fsEc, color='red', label='Fit Err: ' + str(frmse))

plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()