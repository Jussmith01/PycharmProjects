__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import hdnntools as gt
import numpy as np
import matplotlib.pyplot as plt
import time as tm
from scipy import stats as st

def sortbyother(Y, X):
    xy = zip(X, Y)
    xy = sorted(xy, key=lambda x: x[0])
    X, Y = zip(*xy)
    return np.array(Y)

def produce_scan(ax,title,xlabel,cnstfile,saefile,nnfdir,dtdir,dt1,dt2,dt3,smin,smax,iscale,ishift):
    xyz, typ, Eact, chk = gt.readncdat(dtdir + dt1)
    xyz2, typ2, Eact2, chk = gt.readncdat(dtdir + dt2)
    xyz3, typ3, Eact3, chk = gt.readncdat(dtdir + dt3)

    Eact = np.array(Eact)
    Eact2 = np.array(Eact2)
    Eact3 = np.array(Eact3)

    # Construct pyNeuroChem classes
    nc1 = pync.pyNeuroChem(cnstfile, saefile, nnfdir, 1)

    # Set the conformers in NeuroChem
    nc1.setConformers(confs=xyz, types=typ)

    # Print some data from the NeuroChem
    print('1) Number of Atoms Loaded: ' + str(nc1.getNumAtoms()))
    print('1) Number of Confs Loaded: ' + str(nc1.getNumConfs()))

    # Compute Forces of Conformations
    print('Computing energies 1...')
    _t1b = tm.time()
    Ecmp1 = np.array(nc1.computeEnergies())
    print('Computation complete 1. Time: ' + "{:.4f}".format((tm.time() - _t1b) * 1000.0) + 'ms')

    n = smin
    m = smax
    Ecmp1 = gt.hatokcal * Ecmp1
    Eact  = gt.hatokcal * Eact
    Eact2 = gt.hatokcal * Eact2
    Eact3 = gt.hatokcal * Eact3

    IDX = np.arange(0, Eact.shape[0], 1, dtype=float) * iscale + ishift

    IDX = IDX[n:m]
    Eact = Eact[n:m]
    Eact2 = Eact2[n:m]
    Eact3 = Eact3[n:m]
    Ecmp1 = Ecmp1[n:m]

    Ecmp1 = Ecmp1 - Ecmp1.min()
    Eact  = Eact  - Eact.min()
    Eact2 = Eact2 - Eact2.min()
    Eact3 = Eact3 - Eact3.min()

    rmse1 = gt.calculaterootmeansqrerror(Eact, Ecmp1)
    rmse3 = gt.calculaterootmeansqrerror(Eact, Eact2)
    rmse4 = gt.calculaterootmeansqrerror(Eact, Eact3)

    print("Spearman corr. 1: " + "{:.3f}".format(st.spearmanr(Ecmp1, Eact)[0]))
    print("Spearman corr. 2: " + "{:.3f}".format(st.spearmanr(Eact2, Eact)[0]))
    print("Spearman corr. 3: " + "{:.3f}".format(st.spearmanr(Eact3, Eact)[0]))

    ax.plot(IDX, Eact, '-', marker=r'o', color='black', label='DFT',
             linewidth=2, markersize=7)
    ax.plot(IDX, Ecmp1, ':', marker=r'D', color='red', label='ANI-1 RMSE: ' + '%s' % float('%.3g' % rmse1) + ' kcal/mol',
             linewidth=2, markersize=5)
    ax.plot(IDX, Eact2, ':', marker=r'v', color='blue', label='DFTB  RMSE: ' + '%s' % float('%.3g' % rmse3) + ' kcal/mol',
             linewidth=2, markersize=5)
    ax.plot(IDX, Eact3, ':', marker=r'*', color='orange', label='PM6   RMSE: ' + '%s' % float('%.3g' % rmse4) + ' kcal/mol',
             linewidth=2, markersize=7)

    #ax.plot(IDX, Eact, color='black', label='DFT', linewidth=3)
    #ax.scatter(IDX, Eact, marker='o', color='black', linewidth=4)

    th = ax.set_title(title,fontsize=16)
    th.set_position([0.5,1.005])

    # Set Limits
    ax.set_xlim([ IDX.min(),IDX.max()])
    ax.set_ylim([Eact.min()-1.0,Eact.max()+1.0])

    ax.set_ylabel('$\Delta$E calculated (kcal/mol)')
    ax.set_xlabel(xlabel)
    ax.legend(bbox_to_anchor=(0.2, 0.98), loc=2, borderaxespad=0., fontsize=14)

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Research/NeuroChemForceTesting/train_01/'
cnstfile = wkdir + 'rH-3.0A_4-2.5A_a2-2.params'
saefile  = wkdir + '../sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

dtdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_testdata/'

fig, axes = plt.subplots(nrows=1, ncols=1)

produce_scan(axes,'Fentanyl NC bond stretch','Bond distance ($\AA$)'                ,cnstfile,saefile,nnfdir,dtdir,'h2bondscan_test.dat' ,'h2bondscan_test.dat','h2bondscan_test.dat' ,0,998,0.00005,0.6)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

plt.show()
