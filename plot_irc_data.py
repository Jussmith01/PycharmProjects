# Import pyNeuroChem
import pyNeuroChem as pync
import numpy as np
import hdnntools as hdt
import matplotlib.pyplot as plt
import os

ddir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_rxns/benz_dbm_rxns/scans_dbm_benz_4/data/'
ircf = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_rxns/benz_dbm_rxns/scans_dbm_benz_4/irc.dat'

# Set required files for pyNeuroChem
rcdir = '/home/jujuman/Research/ANI-DATASET/RXN1_TNET/training/rxn1to6/ani_benz_rxn_ntwk/'
cnstfile = 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile = 'sae_6-31gd.dat'

files = os.listdir(ddir)
nc = pync.conformers(rcdir + '../../' + cnstfile, rcdir + '../../' + saefile, rcdir + '/networks/', 0)

files = list({("_".join(f.split("_")[:-1])) for f in files})
files = sorted(files, key=lambda x: int(x.split('-')[1].split('.')[0]))
files = [i + '_train.dat' for i in files]

datairc = hdt.readncdat(ircf,type=np.float32)

atm = 9
xr = datairc[0][0][atm]
print(xr)

E1 = []
EA = []
dx = []
for f in files:
    print(ddir+f)
    data = hdt.readncdat(ddir+f,type=np.float32)

    # Set the conformers in NeuroChem
    nc.setConformers(confs=data[0], types=list(data[1]))

    dx.append(np.array([np.linalg.norm(m[atm] - xr) for m in data[0]]))

    # Compute Energies of Conformations
    E1t = nc.energy()
    EAt = data[2]

    E1.append(hdt.hatokcal * E1t)
    EA.append(hdt.hatokcal * EAt)

E1 = np.concatenate(E1)
EA = np.concatenate(EA)
dxl = dx
dx = np.concatenate(dx)

for i,x in enumerate(dxl):
    if i%3 == 0:
        plt.hist(x, bins=50, histtype=u'step')
        #plt.plot(np.array(range(0,x.shape[0])),x)

#plt.ylabel('Rc $(\AA)$')
#plt.xlabel('Step')
plt.ylabel('count')
plt.xlabel('$(\AA)$')
plt.show()
# Plot
errn = hdt.calculaterootmeansqrerror(E1, EA)
plt.scatter(dx,E1, color='red', label="{:.2f}".format(errn), linewidth=1)
plt.scatter(dx,EA, color='black', linewidth=1)
plt.plot(np.array([np.linalg.norm(m[atm] - xr) for m in datairc[0]]),hdt.hatokcal * datairc[2],marker='o', color='blue', linewidth=3)

plt.suptitle("Double bond migration IRCs")

#plt.ylabel('E (kcal/mol)')
#plt.xlabel('Distance $\AA$')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.,fontsize=16)

plt.show()