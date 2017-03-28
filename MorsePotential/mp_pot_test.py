import pyNeuroChem as pync
import pyanitools as pya
import hdnntools as hdt
import numpy as np
import matplotlib.pyplot as plt

#Network 1 Files
wkdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/train_09fsrc/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

#Network 1 Files
wkdir2 = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ccdissotest1-ntwk/'
cnstfile2 = wkdir2 + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile2  = wkdir2 + 'sae_6-31gd.dat'
nnfdir2   = wkdir2 + 'networks/'


# Construct pyNeuroChem classes
nc = pync.conformers(cnstfile, saefile, nnfdir, 1)
nc2 = pync.conformers(cnstfile2, saefile2, nnfdir2, 1)

scan = hdt.readncdat('/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_testdata/waterdimerscan/ethane_diss.dat', type=np.float32)

sdat = [hdt.ncdata(hdt.generatedmatsd3(scan[0]), scan[1], scan[0].shape[1])]

sae = hdt.compute_sae(saefile, scan[1])
serg = scan[2] - sae

# Set the conformers in NeuroChem
nc.setConformers(confs=scan[0],types=list(scan[1]))
nc2.setConformers(confs=scan[0],types=list(scan[1]))

x = 0.1 * np.array(range(serg.shape[0]), dtype=np.float64) + 0.6
print(len(x))

popt = np.load('mp_ani_params_test.npz')['param']
fsEc = hdt.buckingham_pot(sdat, *popt)

aerg = nc.energy() + fsEc - sae
a2erg = nc2.energy() - sae

frmse = hdt.calculaterootmeansqrerror(serg, fsEc)

plt.plot(x, serg, color='black', label='Act')
plt.plot(x, fsEc, color='red', label='Fit Err: ' + str(frmse))
plt.plot(x, aerg, color='green', label='Fit Err: ' + str(frmse))
plt.plot(x, a2erg, color='blue', label='Fit Err: ' + str(frmse))

plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()