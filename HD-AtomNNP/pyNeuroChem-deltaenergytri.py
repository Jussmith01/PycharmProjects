__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Grid
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate
from scipy import stats as st
import time as tm

def flipmat(mat,C):

    mat2 = np.ndarray(shape=(C, C), dtype=float)

    for i in range(C):
        for j in range(C):
            mat2[j, i] = mat[i,C-j-1]

    return mat2

def graphEdiff2D(ax, data, title, min, max):
    x, y, z, d = gt.calculateelementdiff2D(data)

    C = data.shape[0]

    mat = np.ndarray(shape=(C,C),dtype=float)

    for i in range(C):
        for j in range(C):
            I = int(i)
            J = int(j)
            mat[I,J] = z[I+J*C]

    mat = np.transpose(mat)
    mat = flipmat(mat,C)

    #get discrete colormap
    cmap = plt.get_cmap('RdBu', np.max(mat)-np.min(mat)+1)

    # Show mat
    #im = ax.matshow(mat,vmin = np.min(mat), vmax = np.max(mat))
    im = ax.matshow(mat, vmin=min, vmax=max)

    th = ax.set_title(title,fontsize=16)
    th.set_position([0.5,1.005])

    cmap = plt.cm.jet
    #norm = plt.Normalize(mat.min(), mat.max())
    norm = plt.Normalize(min, max)
    rgba = cmap(norm(mat))
    for i in range(C):
        rgba[range(i+1,C), C-i-1, :3] = 1, 1, 1
    ax.imshow(rgba, interpolation='nearest')

    # Plot center line
    ax.plot([x.max()+1-0.5,x.min()-1+0.5],[y.min()-0.5,x.max()+0.5],'--',color='red',linewidth=4,alpha=0.8)

    # Set Limits
    ax.set_xlim([-0.06*x.max(),x.max()+0.06*x.max()])
    ax.set_ylim([-0.06*y.max(),y.max()+0.06*y.max()])

    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)

    #ax.xaxis.tick_bottom()

    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    return im, z


def graphEdiffDelta2D(ax, title, data1, data2, Na, min, max):
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

    mat = np.transpose(mat)
    mat = flipmat(mat, C)

    #get discrete colormap
    cmap = plt.get_cmap('RdBu', np.max(mat)-np.min(mat)+1)

    # Show mat
    im = ax.matshow(mat,vmin = min, vmax = max)

    th = ax.set_title(title,fontsize=16)
    th.set_position([0.5,1.005])

    cmap = plt.cm.jet
    norm = plt.Normalize(min, max)
    rgba = cmap(norm(mat))
    for i in range(C):
        rgba[range(i+1,C), C-i-1, :3] = 1, 1, 1
    ax.imshow(rgba, interpolation='nearest')

    # Plot center line
    ax.plot([x.max() + 1 - 0.5, x.min() - 1 + 0.5], [y.min() - 0.5, x.max() + 0.5], '--', color='red', linewidth=4,alpha=0.8)

    # Set Limits
    ax.set_xlim([-0.06*x.max(),x.max()+0.06*x.max()])
    ax.set_ylim([-0.06*y.max(),y.max()+0.06*y.max()])

    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)

    #ax.xaxis.tick_bottom()
    #ax.yaxis.tick_right()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    return im

font = {'family': 'Bitstream Vera Sans',
        'weight': 'normal',
        'size': 14}

plt.rc('font', **font)

# Set required files for pyNeuroChem
wkdir = '/home/jujuman/Scratch/Research/trainingcases/wB97X-631gd-train-highgarden/train_08-a3.1A_r4.6_dn1/'

# Network  Files
cnstfile = wkdir + 'rHCNO-4.6A_32-3.1A_a8-8.params'
saefile  = wkdir + '../sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

# Construct pyNeuroChem classes
nc = pync.pyNeuroChem(cnstfile,saefile,nnfdir,0)

# Read nc DATA
xyz,typ,Eact,tmp = gt.readncdat('/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/Retinol/data/retinolconformer_DFT.dat')
xyz1,typ1,Eact1,tmp = gt.readncdat('/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/Retinol/data/retinolconformer_DFTB.dat')

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

fig, axes = plt.subplots(nrows=2, ncols=3)
fig.delaxes(axes.flat[3])

x, y, z, d = gt.calculateelementdiff2D(Eact1)

#th = axes.flat[3].set_title("TGM vs. DFT\n$E_T$ correlation", fontsize=16)
#th.set_position([0.5, 1.005])

#axes.flat[3].plot(Eact, Eact, color="black", label= "DFT", linewidth=3)
#axes.flat[3].scatter(Eact, Ecmp, marker=r'o', color="red", label= "TGM", linewidth=3)
#axes.flat[3].legend(bbox_to_anchor=(0.1, 0.8), loc=3, borderaxespad=0.,fontsize=14)

#axes.flat[3].set_ylabel('$E_{cmp}$ (kcal/mol)')
#axes.flat[3].set_xlabel('$E_{ref}$ (kcal/mol)')

#axes.flat[3].set_xlim([Eact.min()-0.5, Ecmp.max()+0.5])
#axes.flat[3].set_ylim([Eact.min()-0.5, Ecmp.max()+0.5])

print (st.linregress(Eact,Ecmp))
print (Eact)
print (Ecmp)

print (np.abs(Eact) - np.abs(Ecmp))

min1 = z.min()
max1 = z.max()

im1,del1 = graphEdiff2D(axes.flat[0], Eact, 'DFT $\Delta E$', min1, max1)
im2,del2 = graphEdiff2D(axes.flat[1], Ecmp, 'ANI-1 $\Delta E$', min1, max1)
im3,del3 = graphEdiff2D(axes.flat[2], Eact1, 'DFTB $\Delta E$', min1, max1)

x, y, z, d = gt.calculateelementdiff2D(Eact)
x2, y2, z2, d2 = gt.calculateelementdiff2D(Eact1)
z = np.abs(z - z2)

min2 = z.min()
max2 = z.max()

im4 = graphEdiffDelta2D(axes.flat[4], 'DFT vs ANI-1\n$|\Delta \Delta E|$', Eact, Ecmp, nc.getNumAtoms(), min2, max2)
im5 = graphEdiffDelta2D(axes.flat[5], 'DFT vs DFTB\n$|\Delta \Delta E|$', Eact, Eact1, nc.getNumAtoms(), min2, max2)

#th = fig.suptitle('$\Delta$E between conformers of Retinol')
#th.set_fontsize(32)
#th.set_position([0.5,0.97])

#tell the colorbar to tick at integers
#cbaxes = fig.add_axes([0.335, 0.55, 0.01, 0.33])
#ch1 = fig.colorbar(im1, ax=axes.ravel().tolist(),cax=cbaxes, orientation='vertical')
#ch1.set_label('$\Delta$E kcal/mol',fontsize=14)

#tell the colorbar to tick at integers
#cbaxes = fig.add_axes([0.605, 0.55, 0.01, 0.33])
#ch1 = fig.colorbar(im2, ax=axes.ravel().tolist(),cax=cbaxes, orientation='vertical')
#ch1.set_label('$\Delta$E kcal/mol',fontsize=14)

#tell the colorbar to tick at integers
cbaxes = fig.add_axes([0.92, 0.54, 0.015, 0.36])
ch1 = fig.colorbar(im3, ax=axes.ravel().tolist(),cax=cbaxes, orientation='vertical')
ch1.set_label('$\Delta$E kcal/mol',fontsize=14)

#cbaxes = fig.add_axes([0.605, 0.115, 0.01, 0.33])
#ch1 = fig.colorbar(im4, ax=axes.ravel().tolist(),cax=cbaxes, orientation='vertical')
#ch1.set_label('|$\Delta \Delta$E| kcal/mol',fontze=14)

cbaxes = fig.add_axes([0.92, 0.105, 0.015, 0.36])
ch1 = fig.colorbar(im5, ax=axes.ravel().tolist(),cax=cbaxes, orientation='vertical')
ch1.set_label('|$\Delta \Delta$E| kcal/mol',fontsize=14)

#fig.tight_layout()

plt.show()