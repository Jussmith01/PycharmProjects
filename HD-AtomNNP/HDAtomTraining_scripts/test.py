import pyNeuroChem as pync
import hdnntools as hdt
from os import listdir
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import *

dir = '/home/jujuman/nn_coor/'

wkdir = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk/'

cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

nc = pync.molecule(cnstfile, saefile, nnfdir, 0)

files = listdir(dir)
files = sorted(files, key=lambda x: int(x.split('_')[2].split('.')[0]))
files = sorted(files, key=lambda x: int(x.split('_')[1].split('.')[0]))

x = np.array(range(-180, 190, 10))

print(x.shape)
#exit(0)
phi = []
psi = []
Eact = []

xyz = []
for n,i in enumerate(files):

    phi.append(i.split('_')[1].split('.')[0])
    psi.append(i.split('_')[2].split('.')[0])

    data = hdt.readxyz2(dir + i)
    print(data[0])
    print(data[1])

    xyz.append(data[0])

    nc.setMolecule(coords=data[0][0],types=list(data[1]))

    Eact.append(nc.energy()[0])

xyz = np.vstack(xyz)
#print(xyz)
hdt.writexyzfile('/home/jujuman/Dropbox/ChemSciencePaper.AER/JustinsDocuments/ACS_april_2017/DipeptidePhiPsi/data.xyz', xyz, list(data[1]))

#data2 = hdt.readncdat('/home/jujuman/Dropbox/ChemSciencePaper.AER/JustinsDocuments/ACS_april_2017/DipeptidePhiPsi/data.dat',type=np.float32)

Eact = np.array(Eact)
#Edft = data2[1]


z = np.array(hdt.hatokcal * (Eact-Eact.min()),dtype=np.float32).reshape(x.shape[0],x.shape[0])
#z2 = np.array(hdt.hatokcal * (Edft-Edft.min()),dtype=np.float32).reshape(x.shape[0],x.shape[0])

Spline=scipy.interpolate.RectBivariateSpline(x,x,z)
XY_New=np.linspace(-180,180,200)

fig=plt.figure()
fig.set_size_inches(6,6)
ax=fig.add_subplot(111)
#plt.pcolormesh(XY_New,XY_New,Spline(XY_New,XY_New).T,cmap=plt.get_cmap('Spectral_r'))
plt.pcolormesh(XY_New,XY_New,Spline(XY_New,XY_New).T,cmap=plt.get_cmap('jet'))
#plt.colorbar()
CS=plt.contour(XY_New,XY_New,Spline(XY_New,XY_New).T,10,colors='k',linewidth=0.5)
plt.clabel(CS,inline=1, fontsize=7)
plt.axis([-180, 180, -180, 180])

plt.show()