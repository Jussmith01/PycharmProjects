__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate
import time as tm

def getmaxmin(data):
    x, y, z, d = gt.calculateelementdiff2D(data)
    z = gt.convert * z
    return z.max(), z.min()


def graphEdiff2D(ax, data1, data2, title):
    x, y, z, d = gt.calculatecompareelementdiff2D(data1, data2)

    #z = gt.convert * z

    # Set up a regular grid of interpolation points
    xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
    xi, yi = np.meshgrid(xi, yi)

    # Interpolate
    rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
    zi = rbf(xi, yi)

    im = ax.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
                   extent=[x.min(), x.max(), y.min(), y.max()])

    ax.set_title(title,fontsize=14)
    # plt.xlabel('Structure')
    # plt.ylabel('Structure')

    ax.scatter(x, y, c=z, vmin=z.min(), vmax=z.max())
    ax.plot([x.min(),x.max()],[y.min(),x.max()],'--',color='red',linewidth=4,alpha=0.8)

    ax.set_xlim([-0.02*x.max(),x.max()+0.02*x.max()])
    ax.set_ylim([-0.02*y.max(),y.max()+0.02*y.max()])

    ax.grid(True)

    return im


def graphEdiffDelta2D(ax, data1, data2, Na):
    #data1 = gt.convert * data1
    #data2 = gt.convert * data2

    x, y, z, d = gt.calculateelementdiff2D(data1)
    x2, y2, z2, d2 = gt.calculateelementdiff2D(data2)

    RMSE = gt.calculaterootmeansqrerror(d,d2) / float(Na)

    print ('dataz1:',z)
    print ('dataz2:',z2)

    z = np.abs(z - z2)

    print ('zdiff:',z)

    # Set up a regular grid of interpolation points
    xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
    xi, yi = np.meshgrid(xi, yi)

    # Interpolate
    rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
    zi = rbf(xi, yi)

    im = ax.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
                   extent=[x.min(), x.max(), y.min(), y.max()])

    ax.set_title('Difference Plot (symmetric)\nRMSE: ' + "{:.5f}".format(RMSE) + 'kcal/mol/atom Atoms: ' + str(Na),fontsize=14)
    # plt.xlabel('Structure')
    # plt.ylabel('Structure')

    ax.scatter(x, y, c=z, vmin=z.min(), vmax=z.max())

    ax.set_xlim([-0.02*x.max(),x.max()+0.02*x.max()])
    ax.set_ylim([-0.02*y.max(),y.max()+0.02*y.max()])

    #ax.plot([x.min(),x.max()],[y.min(),x.max()],color='red',linewidth=4)
    ax.grid(True)

    return im


# Set required files for pyNeuroChem
wkdir1    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/train_08_2/'
wkdir2    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/train_08_3/'

#Network 1 Files
cnstfile1 = wkdir1 + 'rHCNO-4.6A_32-3.1A_a8-8.params'
saefile1  = wkdir1 + 'sae_6-31gd.dat'
nnfdir1   = wkdir1 + 'networks/'

# Network 2 Files
cnstfile2 = wkdir2 + 'rHCNO-4.6A_32-3.1A_a8-8.params'
saefile2  = wkdir2 + 'sae_6-31gd.dat'
nnfdir2   = wkdir2 + 'networks/'

# Construct pyNeuroChem classes
nc1 = pync.pyNeuroChem(cnstfile1,saefile1,nnfdir1,0)
nc2 = pync.pyNeuroChem(cnstfile2,saefile2,nnfdir2,0)

xyz,typ,Eact = gt.readncdat('/home/jujuman/Dropbox/Research/ChemSciencePaper/TestCases/C10H20Isomers/isomer_structures_DFT.dat')

Eact = np.array(Eact)

# Set the conformers in NeuroChem
nc1.setConformers(confs=xyz,types=typ)
nc2.setConformers(confs=xyz,types=typ)

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(nc1.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(nc1.getNumConfs()) )

# Compute Forces of Conformations
print('Computing energies 1...')
_t1b = tm.time()
Ecmp1 = np.array( nc1.computeEnergies() )
print('Computation complete 1. Time: ' + "{:.4f}".format((tm.time() - _t1b) * 1000.0)  + 'ms')

print('Computing energies 2...')
_t2b = tm.time()
Ecmp2 = np.array( nc2.computeEnergies() )
print('Computation complete 2. Time: ' + "{:.4f}".format((tm.time() - _t2b) * 1000.0) + 'ms')

Ecmp1 = gt.hatokcal * Ecmp1
Ecmp2 = gt.hatokcal * Ecmp2
Eact = gt.hatokcal * Eact

#Ecmp1 = Ecmp1 - gt.calculatemean(Ecmp1)
#Ecmp2 = Ecmp2 - gt.calculatemean(Ecmp2)
#Eact  = Eact  - gt.calculatemean(Eact)

IDX = np.arange(0,Eact.shape[0],1,dtype=float)

fig, axes = plt.subplots(nrows=2, ncols=2)

im1 = graphEdiff2D(axes.flat[0], Eact, Ecmp1, 'Top left: ANN - c07b\nBottom right: wB97x/6-31g*')
im2 = graphEdiffDelta2D(axes.flat[1], Eact, Ecmp1, nc1.getNumAtoms())
im3 = graphEdiff2D(axes.flat[2], Eact, Ecmp2, 'Top left: ANN - c08b\nBottom right: wB97x/6-31g*')
im4 = graphEdiffDelta2D(axes.flat[3], Eact, Ecmp2, nc2.getNumAtoms())

#th = fig.suptitle('300K norm. mode generated conformers of Benzamide\n(Rcr=4.5$\AA$)')
#th = fig.suptitle('300K norm. mode generated conformers of\nH-Gly-Pro-Hyp-Gly-Ala-Gly-OH')
th = fig.suptitle('Energy difference between conformers\nof Retinol')
#th = fig.suptitle('Linearly interpolated energy difference\nbetween Isomers of $\\rm C_8 \\rm H_{16}$')
th.set_fontsize(26)
th.set_position([0.5,0.98])

#fig.text(0.3, 0.515, 'Network tag: ' + networktag1,fontsize=10,verticalalignment='center',horizontalalignment='center')
#fig.text(0.3, 0.075, 'Network tag: ' + networktag2,fontsize=10,verticalalignment='center',horizontalalignment='center')

cbaxes = fig.add_axes([0.4, 0.54, 0.01, 0.36])
ch1 = fig.colorbar(im1, ax=axes.ravel().tolist(),cax=cbaxes, orientation='vertical')
ch1.set_label('$\Delta$E kcal/mol',fontsize=12)

cbaxes = fig.add_axes([0.6, 0.54, 0.01, 0.36])
ch2 = fig.colorbar(im2, ax=axes.ravel().tolist(),cax=cbaxes, orientation='vertical')
ch2.set_label('|$\Delta\Delta$E| kcal/mol',fontsize=12)
ch2.ax.yaxis.set_ticks_position('left')

cbaxes = fig.add_axes([0.4, 0.1, 0.01, 0.36])
ch3 = fig.colorbar(im3, ax=axes.ravel().tolist(),cax=cbaxes, orientation='vertical')
ch3.set_label('$\Delta$E kcal/mol',fontsize=12,)

cbaxes = fig.add_axes([0.6, 0.1, 0.01, 0.36])
ch4 = fig.colorbar(im4, ax=axes.ravel().tolist(),cax=cbaxes, orientation='vertical')
ch4.set_label('|$\Delta\Delta$E| kcal/mol',fontsize=12)
ch4.ax.yaxis.set_ticks_position('left')

plt.show()