import pyNeuroChem as pync
import hdnntools as hdt
import numpy as np
import matplotlib.pyplot as plt

#Network 1 Files
wkdir = '/home/jujuman/Dropbox/ChemSciencePaper.AER/networks/ANI-c08f-ntwk/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

wkdir = '/home/jujuman/Research/CrossValidation/GDB-09-Retrain/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile = wkdir + 'sae_6-31gd.dat'

traj = hdt.readncdat('/home/jujuman/Dropbox/ChemSciencePaper.AER/JustinsDocuments/ACS_april_2017/ANI-MD-vDFT/benzeneMD/benzene_dft.dat',type=np.float32)

nc = pync.conformers(cnstfile, saefile, nnfdir, 0)
ncl =  [pync.conformers(cnstfile, saefile, wkdir + 'cv_c08e_ntw_' + str(l) + '/networks/', 0) for l in range(5)]

nc.setConformers(confs=traj[0],types=list(traj[1]))

Ecv = []
for comp in ncl:
    comp.setConformers(confs=traj[0], types=list(traj[1]))
    E = comp.energy()
    Ecv.append(E-E.min())

Ecv = np.vstack(Ecv)
std = np.std(Ecv, axis=0)
print(std)

Eqm = hdt.hatokcal * traj[2]
Ean = hdt.hatokcal * nc.energy()

Eqm = Eqm - Eqm.min()
Ean = Ean - Ean.min()

print(hdt.calculaterootmeansqrerror(Eqm,Ean))

# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(Eqm, color='black', label='QM')
axarr[0].plot(Ean, color='red', label='ANI')

for p in Ecv:
    axarr[0].plot(hdt.hatokcal * p, color='grey', label='ANI')

axarr[0].set_xlabel('E')

axarr[0].set_title('Sharing X axis')
axarr[1].plot(hdt.hatokcal * std, color='black', label='CV Std. Dev.')
axarr[1].plot(np.abs(Eqm - Ean), color='red', label='Delta E')

axarr[1].set_xlabel('t')
axarr[1].set_ylabel('Std. Dev.')
plt.legend()
plt.show()
