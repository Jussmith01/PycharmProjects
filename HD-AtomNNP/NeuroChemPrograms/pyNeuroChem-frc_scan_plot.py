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

def produce_scan(title,xlabel,cnstfile,saefile,nnfdir,dtdir,dt1,smin,smax,iscale,ishift,atm):
    xyz, frc, typ, Eact, chk = gt.readncdatwforce(dtdir + dt1)

    xyz = np.asarray(xyz,dtype=np.float32)
    xyz = xyz.reshape((xyz.shape[0],len(typ),3))

    print(xyz)

    frc = np.asarray(frc,dtype=np.float32)
    frc = frc.reshape((frc.shape[0],len(typ),3))

    Eact = np.array(Eact)

    # Construct pyNeuroChem classes
    nc1 = pync.pyNeuroChem(cnstfile, saefile, nnfdir, 0)

    # Set the conformers in NeuroChem
    nc1.setConformers(confs=xyz, types=typ)

    # Print some data from the NeuroChem
    print('1) Number of Atoms Loaded: ' + str(nc1.getNumAtoms()))
    print('1) Number of Confs Loaded: ' + str(nc1.getNumConfs()))

    # Compute Energies of Conformations
    print('Computing energies...')
    _t1b = tm.time()
    Ecmp1 = nc1.energy()
    print('Energy computation complete. Time: ' + "{:.4f}".format((tm.time() - _t1b) * 1000.0) + 'ms')

    # Compute Forces of Conformations
    print('Compute forces...')
    _t1b = tm.time()
    F = nc1.force()
    print('Force computation complete. Time: ' + "{:.4f}".format((tm.time() - _t1b) * 1000.0) + 'ms')

    #Fn = np.array(nc1.computeNumericalForces(dr=0.0001))

    n = smin
    m = smax
    Ecmp1 = gt.hatokcal * Ecmp1
    Eact  = gt.hatokcal * Eact

    IDX = np.arange(0, Eact.shape[0], 1, dtype=float) * iscale + ishift

    IDX = IDX
    Eact = Eact
    Ecmp1 = Ecmp1

    Ecmp1 = Ecmp1 - Ecmp1.min()
    Eact  = Eact  - Eact.min()

    rmse1 = gt.calculaterootmeansqrerror(Eact, Ecmp1)

    print("Spearman corr. 1: " + "{:.3f}".format(st.spearmanr(Ecmp1, Eact)[0]))

    fig, axes = plt.subplots(nrows=2, ncols=2)

    axes.flat[0].plot(IDX, Eact, '-', color='black', label='DFT',
             linewidth=6)
    axes.flat[0].plot(IDX, Ecmp1, '--', color='red', label='ANI-1',
             linewidth=6)

    #ax.plot(IDX, Eact, color='black', label='DFT', linewidth=3)
    #ax.scatter(IDX, Eact, marker='o', color='black', linewidth=4)

    th = axes.flat[0].set_title(title,fontsize=13)
    th.set_position([0.5,1.00])

    # Set Limits
    axes.flat[0].set_xlim([ IDX.min(),IDX.max()])
    #axes.flat[0].set_ylim([Eact.min()-1.0,Eact.max()+1.0])

    axes.flat[0].set_ylabel('$\Delta$E calculated (kcal/mol)')
    axes.flat[0].set_xlabel(xlabel)
    axes.flat[0].legend(bbox_to_anchor=(0.2, 0.98), loc=2, borderaxespad=0., fontsize=12)

    #print (Fn)

    for i in range(0,3):
        Fq = F[:,:,i][:,atm]
        #Fnq = Fn[:,i]
        Faq = (1.8897259885789*np.array(frc)[:,:,i][:,atm])
        print (Fq)
        print (Faq)

        th = axes.flat[i+1].set_title("Force for atom " + str(atm) + ": coordinate " + str(i),fontsize=13)
        th.set_position([0.5,1.00])

        axes.flat[i+1].set_xlim([IDX.min(), IDX.max()])

        axes.flat[i+1].plot(IDX, Faq, '-', color='black', label='DFT',
            linewidth=6)
        #axes.flat[i+1].plot(IDX, Fnq, '-', color='blue', label='ANI Numerical',
        #    linewidth=6)
        axes.flat[i+1].plot(IDX, Fq, '--', color='red', label='ANI Analytical',
            linewidth=6)

        axes.flat[i+1].set_ylabel('Force (Ha/A)')
        axes.flat[i+1].set_xlabel(xlabel)
        axes.flat[i+1].legend(bbox_to_anchor=(0.2, 0.98), loc=2, borderaxespad=0., fontsize=12)

    font = {'family' : 'Bitstream Vera Sans',
            'weight' : 'normal',
            'size'   : 12}

    plt.rc('font', **font)

    plt.show()

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Research/NeuroChemForceTesting/train_02/'
cnstfile = wkdir + 'rHO-3.0A_4-2.5A_a2-2.params'
saefile  = wkdir + '../sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

test = np.array([[[1.0,2.0,3.0],[3.0,4.0,5.0]]
                ,[[5.0,6.0,7.0],[7.0,8.0,9.0]]],dtype=np.float32)

print(test[:,:,0][:,1])

dtdir = '/home/jujuman/Research/NeuroChemForceTesting/'

produce_scan('H2O bond stretch (O-H)','Bond distance ($\AA$)'                ,cnstfile,saefile,nnfdir,dtdir,'trainingData.dat' ,0,249,0.001,0.85,0)
#produce_scan('H2O bond stretch (O-H)','Bond distance ($\AA$)'                ,cnstfile,saefile,nnfdir,dtdir,'OH_bond_h2o.dat' ,0,249,0.001,0.85,1)
#produce_scan('H2O bond stretch (O-H)','Bond distance ($\AA$)'                ,cnstfile,saefile,nnfdir,dtdir,'OH_bond_h2o.dat' ,0,249,0.001,0.85,2)
