__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
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
    xyz, typ, Eact, chk = gt.readncdat(dtdir + dt1,np.float32)
    xyz2, typ2, Eact2, chk = gt.readncdat(dtdir + dt2)
    xyz3, typ3, Eact3, chk = gt.readncdat(dtdir + dt3)

    #gt.writexyzfile("/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/Dihedrals/4-Cyclohexyl-1-butanol/optimization/dihedral_"+dt1+".xyz",xyz,typ)

    #Eact = np.array(Eact)
    #Eact2 = np.array(Eact2)
    #Eact3 = np.array(Eact3)

    # Construct pyNeuroChem classes
    nc1 = pync.pyNeuroChem(cnstfile, saefile, nnfdir, 0)

    # Set the conformers in NeuroChem
    nc1.setConformers(confs=xyz, types=typ)

    # Print some data from the NeuroChem
    print('1) Number of Atoms Loaded: ' + str(nc1.getNumAtoms()))
    print('1) Number of Confs Loaded: ' + str(nc1.getNumConfs()))

    # Compute Forces of Conformations
    print('Computing energies 1...')
    _t1b = tm.time()
    Ecmp1 = nc1.energy()
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

#Network 1 Files
wkdir    = '/home/jujuman/Research/wB97X-631gd-train-highgarden/train_08-a3.1A_r4.6_AEV384_1/'
cnstfile1 = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile1  = wkdir + '../sae_6-31gd.dat'
nnfdir1   = wkdir + 'networks/'

dtdir1 = '/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/Bonds/Fentanyl/'
dtdir2 = '/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/Bonds/Fentanyl/'
dtdir3 = '/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/Dihedrals/4-Cyclohexyl-1-butanol/'
dtdir4 = '/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/Dihedrals/Lisdexamfetamine/'

#dtdir4 = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_testdata/butane/'

fig, axes = plt.subplots(nrows=2, ncols=2)

produce_scan(axes.flat[0],'Fentanyl NC bond stretch'                     ,'Bond distance ($\AA$)',cnstfile1,saefile1,nnfdir1,dtdir1,'data_bond_scan_DFT.dat'  ,'data_bond_scan_DFTB.dat' ,'data_bond_scan_PM6.dat' ,10,50,0.0025,1.3)
produce_scan(axes.flat[1],'Fentanyl CCC angle bend'                      ,'Angle ($^\circ$)'     ,cnstfile1,saefile1,nnfdir1,dtdir2,'data_angle_scan_DFT.dat' ,'data_angle_scan_DFTB.dat','data_angle_scan_PM6.dat',10,180,0.5,80.0)
produce_scan(axes.flat[2],'4-Cyclohexyl-1-butanol CCCC dihedral rotation','Angle ($^\circ$)'     ,cnstfile1,saefile1,nnfdir1,dtdir3,'dhl_scan_DFT.dat'        ,'dhl_scan_DFTB.dat'       ,'dhl_scan_PM6.dat'       ,0,144,2.5,0.0)
produce_scan(axes.flat[3],'Lisdexamfetamine NCCC dihedral rotation'      ,'Angle ($^\circ$)'     ,cnstfile1,saefile1,nnfdir1,dtdir4,'data_dhl_scan_DFT.dat'   ,'data_dhl_scan_DFTB.dat'  ,'data_dhl_scan_PM6.dat'  ,0,144,2.5,0.0)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

plt.show()
