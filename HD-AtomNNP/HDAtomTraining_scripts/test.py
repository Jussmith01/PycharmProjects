import pyNeuroChem as pync
import hdnntools as hdt
from os import listdir
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import *

dir = '/home/jujuman/nn_coor/'

wkdir = '/home/jujuman/Dropbox/ChemSciencePaper.AER/networks/ANI-c08f-R1-ntwk-cv/cv_c08e_ntw_3/'

cnstfile = wkdir + '../rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = wkdir + '../sae_6-31gd.dat'
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

data2 = hdt.readncdat('/home/jujuman/Dropbox/ChemSciencePaper.AER/JustinsDocuments/ACS_april_2017/DipeptidePhiPsi/data.dat',type=np.float32)

Eact = np.array(Eact)
Edft = data2[2]

z = np.array(hdt.hatokcal * (Eact-Eact.min()),dtype=np.float32).reshape(x.shape[0],x.shape[0])
z2 = np.array(hdt.hatokcal * (Edft-Edft.min()),dtype=np.float32).reshape(x.shape[0],x.shape[0])


rmse = hdt.calculaterootmeansqrerror(z,z2)
print('RMSE: ',rmse)

Spline1=scipy.interpolate.RectBivariateSpline(x,x,z)
Spline2=scipy.interpolate.RectBivariateSpline(x,x,z2)
Spline3=scipy.interpolate.RectBivariateSpline(x,x,abs(z-z2))
XY_New=np.linspace(-180,180,200)

f, ax = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True)
f.set_size_inches(12,6)

maxi = np.concatenate([z,z2]).max()
mini = np.concatenate([z,z2]).min()

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 18}


ax.flat[0].pcolormesh(XY_New,XY_New,Spline1(XY_New,XY_New).T,cmap=plt.get_cmap('jet'), vmin = mini, vmax = maxi)
CS1=ax.flat[0].contour(XY_New,XY_New,Spline1(XY_New,XY_New).T,10,colors='k',linewidth=0.5)
ax.flat[0].clabel(CS1,inline=1, fontsize=7)
ax.flat[0].axis([-180, 180, -180, 180])
ax.flat[0].set_title("ANI")
ax.flat[0].set_xlabel('$\phi$', fontdict=font)
ax.flat[0].set_ylabel('$\psi$', fontdict=font)

test = ax.flat[1].pcolormesh(XY_New,XY_New,Spline2(XY_New,XY_New).T,cmap=plt.get_cmap('jet'), vmin = mini, vmax = maxi)
CS2=ax.flat[1].contour(XY_New,XY_New,Spline2(XY_New,XY_New).T,10,colors='k',linewidth=0.5)
ax.flat[1].clabel(CS2,inline=1, fontsize=7)
ax.flat[1].axis([-180, 180, -180, 180])
ax.flat[1].set_title("DFT")
ax.flat[1].set_xlabel('$\phi$', fontdict=font)
ax.flat[1].set_ylabel('$\psi$', fontdict=font)

ax.flat[2].pcolormesh(XY_New,XY_New,Spline3(XY_New,XY_New).T,cmap=plt.get_cmap('jet'), vmin = mini, vmax = maxi)
CS3=ax.flat[2].contour(XY_New,XY_New,Spline3(XY_New,XY_New).T,10,colors='k',linewidth=0.5,)
ax.flat[2].clabel(CS3,inline=1, fontsize=7)
ax.flat[2].axis([-180, 180, -180, 180])
ax.flat[2].set_title("$|\Delta E|$ ANI vs DFT - " + 'RMSE: ' + "{:.3f}".format(rmse) +' RMSE/atom: ' + "{:.3f}".format(rmse/float(len(list(data[1])))))
ax.flat[2].set_xlabel('$\phi$', fontdict=font)
ax.flat[2].set_ylabel('$\psi$', fontdict=font)

#plt.suptitle("Phi Psi rotation of serine-alanine dipeptide (kcal/mol)", fontsize=18, fontweight='bold')

cbaxes = f.add_axes([0.92, 0.105, 0.015, 0.8])
ch = f.colorbar(test, ax=ax.ravel().tolist(),cax=cbaxes, orientation='vertical')
ch.set_label('$\Delta$E kcal/mol',fontsize=14)

plt.show()