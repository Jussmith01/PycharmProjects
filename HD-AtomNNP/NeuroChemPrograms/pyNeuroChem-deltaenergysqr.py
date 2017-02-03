__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate
import time as tm

def graphEdiff2D(ax, data1, data2, title):
    x, y, z, d = gt.calculatecompareelementdiff2D(data1, data2)

    C = data1.shape[0]

    mat = np.ndarray(shape=(C,C),dtype=float)

    for i in x:
        for j in y:
            I = int(i)
            J = int(j)
            mat[I,J] = z[I+J*C]

    #get discrete colormap
    cmap = plt.get_cmap('RdBu', np.max(mat)-np.min(mat)+1)

    # Show mat
    im = ax.matshow(mat,vmin = np.min(mat), vmax = np.max(mat))

    th = ax.set_title(title,fontsize=16)
    th.set_position([0.5,1.005])

    cmap = plt.cm.jet
    norm = plt.Normalize(mat.min(), mat.max())
    rgba = cmap(norm(mat))
    rgba[range(C), range(C), :3] = 1, 1, 1
    ax.imshow(rgba, interpolation='nearest')

    # Plot center line
    ax.plot([x.min()-0.5,x.max()+0.5],[y.min()-0.5,x.max()+0.5],'--',color='red',linewidth=4,alpha=0.8)

    # Set Limits
    ax.set_xlim([-0.02*x.max(),x.max()+0.02*x.max()])
    ax.set_ylim([-0.02*y.max(),y.max()+0.02*y.max()])

    ax.xaxis.tick_bottom()

    return im


def graphEdiffDelta2D(ax, data1, data2, Na):
    #data1 = gt.convert * data1
    #data2 = gt.convert * data2

    x, y, z, d = gt.calculateelementdiff2D(data1)
    x2, y2, z2, d2 = gt.calculateelementdiff2D(data2)

    RMSE = gt.calculaterootmeansqrerror(d,d2) / float(Na)

    print ('Number of atoms: ' + str(Na))
    print ('RMSE: ' + str(RMSE) + ' kcal/mol/atom')
    print ('RMSE: ' + str(float(Na)*RMSE) + ' kcal/mol')

    z = np.abs(z - z2)

    C = data1.shape[0]

    mat = np.ndarray(shape=(C,C),dtype=float)

    for i in x:
        for j in y:
            I = int(i)
            J = int(j)
            mat[J,I] = z[J+I*C]

    #get discrete colormap
    cmap = plt.get_cmap('RdBu', np.max(mat)-np.min(mat)+1)

    # Show mat
    im = ax.matshow(mat,vmin = np.min(mat), vmax = np.max(mat))

    cmap = plt.cm.jet
    norm = plt.Normalize(mat.min(), mat.max())
    rgba = cmap(norm(mat))
    rgba[range(C), range(C), :3] = 1, 1, 1
    ax.imshow(rgba, interpolation='nearest')

    # Plot center line
    #ax.plot([x.min()-0.5,x.max()+0.5],[y.min()-0.5,x.max()+0.5],color='black',linewidth=4,alpha=0.9)

    # Set Limits
    ax.set_xlim([-0.02*x.max(),x.max()+0.02*x.max()])
    ax.set_ylim([-0.02*y.max(),y.max()+0.02*y.max()])

    ax.xaxis.tick_bottom()

    return im

font = {'family': 'Bitstream Vera Sans',
        'weight': 'normal',
        'size': 14}

plt.rc('font', **font)

# Set required files for pyNeuroChem
wkdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/train_08_5/'

# Network  Files
cnstfile = wkdir + 'rHCNO-4.7A_32-3.2A_a8-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

# Construct pyNeuroChem classes
nc = pync.pyNeuroChem(cnstfile,saefile,nnfdir,0)

# Read nc DATA
xyz,typ,Eact = gt.readncdat('/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/Retinol/data/retinolconformer_DFT.dat')
xyz1,typ1,Eact1 = gt.readncdat('/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/Retinol/data/retinolconformer_DFTB.dat')

Eact = np.array(Eact)
Eact1 = np.array(Eact1)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=typ)

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(nc.getNumConfs()) )

# Compute Forces of Conformations
print('Computing energies...')
_t1b = tm.time()
Ecmp = np.array( nc.computeEnergies() )
print('Computation complete. Time: ' + "{:.4f}".format((tm.time() - _t1b) * 1000.0)  + 'ms')

Ecmp = gt.hatokcal * Ecmp
Eact = gt.hatokcal * Eact

Eact1 = gt.hatokcal * Eact1

IDX = np.arange(0,Eact.shape[0],1,dtype=float)

fig, axes = plt.subplots(nrows=1, ncols=2)

im1 = graphEdiff2D(axes.flat[0], Eact, Ecmp, 'Top left: wB97x/6-31g* - Bottom right: AM1')
im2 = graphEdiffDelta2D(axes.flat[1], Eact, Ecmp, nc.getNumAtoms())
#im3 = graphEdiff2D(axes.flat[2], Eact, Ecmp2, 'Top left: ANN - c08b\nBottom right: wB97x/6-31g*')
#im4 = graphEdiffDelta2D(axes.flat[3], Eact, Ecmp2, nc2.getNumAtoms())

th = fig.suptitle('$\Delta$E between conformers of Retinol')
th.set_fontsize(32)
th.set_position([0.5,0.95])

#tell the colorbar to tick at integers
cbaxes = fig.add_axes([0.48, 0.15, 0.015, 0.7])
ch1 = fig.colorbar(im1, ax=axes.ravel().tolist(),cax=cbaxes, orientation='vertical')
ch1.set_label('$\Delta$E kcal/mol',fontsize=14)

cbaxes = fig.add_axes([0.905, 0.15, 0.015, 0.7])
ch1 = fig.colorbar(im2, ax=axes.ravel().tolist(),cax=cbaxes, orientation='vertical')
ch1.set_label('|$\Delta \Delta$E| kcal/mol',fontsize=14)

plt.show()