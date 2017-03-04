import numpy as np
import pyNeuroChem as pync
import hdnntools as hdn
import pyanitools as pyt

def MAE (act,pre):
    N = act.shape[0]
    e = (np.abs(pre-act)).sum()
    return e / float(N)

# Set required files for pyNeuroChem
anipath  = '/home/jujuman/Research/QM-7TEST/tester/ANI-QM7-ntwk'
cnstfile = anipath + '/rHCNOS-5.0A_16-3.1A_a4-8.params'
saefile  = anipath + '/../sae_6-31gd.dat'
nnfdir   = anipath + '/networks/'

path = "/home/jujuman/Scratch/Research/QM-7TEST/QM7-test-ho.h5"
datas = pyt.anidataloader(path)
datas.totalload()

# Construct pyNeuroChem class
nc = pync.conformers(cnstfile, saefile, nnfdir, 0)

Ea = np.zeros(datas.size())
Ec = np.zeros(datas.size())

for i in range(datas.size()):

    print(i, ' of ', datas.size())
    data = datas.getdata(i)

    x = data[0]
    e = data[1]
    s = data[2]

    Na = s.shape[0]

    xyz_t = np.array(x, dtype = np.float32, order='C').reshape(1,Na,3)
    spc_t = s

    nc.setConformers(confs=xyz_t,types=list(spc_t))

    # Compute Energies of Conformations
    ec = nc.energy()
    Ec[i] = ec[0]
    Ea[i] = e
    #print(ec,' : ',e)

err = MAE(Ea,Ec)
print('MSE(Ha):       ',err)
print('MSE(Kcal/mol): ',hdn.hatokcal*err)
print('MSE(eV):       ',hdn.evtokcal*err)