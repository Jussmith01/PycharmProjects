__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time as tm
from scipy import stats as st
from os import listdir
import os.path

def sortbyother(Y, X):
    xy = zip(X, Y)
    xy = sorted(xy, key=lambda x: x[0])
    X, Y = zip(*xy)
    return np.array(Y)

def setmaxE(X,Y,E):
    m = min(X)
    newlist = []
    for i in range(0,len(X)):
        if gt.hatokcal * abs(X[i]-m) <= E:
            newlist.append(Y[i])

    return newlist

def combinenparray (X):
    tlist = []
    for i in X:
        for j in i:
            tlist.append(j)

    return np.array(tlist,dtype=float)

def combinenparrayrange (X,a,b,Eidx):
    tlist = []
    for i in range(a,b):
        for j in X[Eidx[i]]:
            tlist.append(j)

    return np.array(tlist,dtype=float)

def plt_by_index(Edat,idx):
    IDX = np.arange(0, Edat.shape[0], 1, dtype=float)

    plt.plot(IDX, Edat, color='black', label='DFT', linewidth=3)
    plt.scatter(IDX, Edat, marker='o', color='black', s=10)

    # plt.title("300K NMS structures of\nNME-Gly-Pro-Hyp-Gly-Ala-Gly-ACE")
    plt.title("Errors for molecule: " + str(i))

    plt.ylabel('E calculated (Ha)')
    plt.xlabel('Structure Index')
    plt.legend(bbox_to_anchor=(0.5, 0.25), loc=2, borderaxespad=0., fontsize=14)

    font = {'family': 'Bitstream Vera Sans',
            'weight': 'normal',
            'size': 14}

    plt.rc('font', **font)

    plt.show()

def corrEtotplot(ax,Eact,Ecmp,a,b,label):

    d1 = combinenparrayrange(Eact, a, b, Eidx)
    d2 = combinenparrayrange(Ecmp, a, b, Eidx)

    rmse = gt.calculaterootmeansqrerror(d1, d2)
    slope, intercept, r_value, p_value, std_err = st.linregress(d1, d2)

    lin = 'Slope: ' + '%s' % float('%.3g' % slope) + ' Int.: ' + '%s' % float('%.3g' % intercept) + ' $r^2$: ' + '%s' % float('%.3g' % r_value**2)

    ax.plot(d1,d1,color='black',label='DFT',linewidth=2,)
    ax.scatter(d1, d2, color='red',label=label + ' RMSE ' + '%s' % float('%.3g' % rmse) + ' kcal/mol ' + lin, marker=r"o",s=50)

    # Set Limits
    ax.set_xlim([d1.min(),d1.max()])
    ax.set_ylim([d1.min(),d1.max()])

    font = {'family': 'Bitstream Vera Sans',
            'weight': 'heavy',
            'size': 24}

    ax.set_ylabel('$E_{cmp}$',fontdict=font)
    ax.set_xlabel('$E_{ref}$',fontdict=font)
    ax.legend(bbox_to_anchor=(0.01, 0.99), loc=2, borderaxespad=0., fontsize=10)

    t_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(t_formatter)
    ax.yaxis.set_major_formatter(t_formatter)


# Set data fields
dtdir =  '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_10/testdata/'
fpref = 'gdb11_10-'
fpost = '_test.dat'
rng = [0,140]

files = listdir(dtdir)
#print (files)
# Set required files for pyNeuroChem

#Network 1 Files
wkdir1    = '/home/jujuman/Scratch/Research/trainingcases/wB97X-631gd-train-highgarden/train_08-a3.1A_r4.6_dn1/'
cnstfile1 = wkdir1 + 'rHCNO-4.6A_32-3.1A_a8-8.params'
saefile1  = wkdir1 + '../sae_6-31gd.dat'
nnfdir1   = wkdir1 + 'networks/'

# Construct pyNeuroChem classes
nc = pync.pyNeuroChem(cnstfile1, saefile1, nnfdir1, 0)

Ecmp = []
Eact = []

Eidx = []
Emin = []
Efle = []

Eerr = []
time = 0.0

ld = [0,0.0]
sd = [0,100000.0]

err = []
sze = []

N = 0

for i in files:
#for i in range(rng[0],rng[1]):
    #xyz,typ,Eact_t,readf    = gt.readncdat(dtdir + fpref + str(i) + fpost)
    xyz,typ,Eact_t,readf    = gt.readncdat(dtdir + i)

    if readf:
        # Set the conformers in NeuroChem
        nc.setConformers(confs=xyz,types=typ)

        #print('FILE: ' + dtdir + fpref + str(i) + fpost)
        print('FILE: ' + dtdir + i)

        # Print some data from the NeuroChem
        print( '1) Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
        print( '1) Number of Confs Loaded: ' + str(nc.getNumConfs()) )

        # Compute Energies of Conformations
        print('Computing energies...')
        _t1b = tm.time()
        Ecmp_t = nc.computeEnergies()
        _t2b = (tm.time() - _t1b) * 1000.0
        print('Computation complete. Time: ' + "{:.4f}".format(_t2b)  + 'ms')

        Ecmp_t = setmaxE(Eact_t, Ecmp_t, 300.0)
        Eact_t = setmaxE(Eact_t, Eact_t, 300.0)

        tNa = nc.getNumAtoms()
        err.append(gt.hatokcal * gt.calculaterootmeansqrerror(np.array(Eact_t, dtype=float),np.array(Ecmp_t, dtype=float)) / float(tNa))
        sze.append(float(len(Eact_t)))

        time += _t2b

        Eidx.append(N)
        Emin.append(gt.hatokcal * np.array( Ecmp_t ).min())
        Efle.append(i)

        N += 1
        Ecmp.append( gt.hatokcal * np.array( Ecmp_t ) )
        Eact.append( gt.hatokcal * np.array( Eact_t ) )

err = np.array(err, dtype=float)
sze = np.array(sze, dtype=float)

Eidx = sortbyother(Eidx,Emin)
Efle = sortbyother(Efle,Emin)
Emin = sortbyother(Emin,Emin)

#slope1, intercept1, r_value1, p_value1, std_err1 = st.linregress(Eact,Ecmp)

fig, axes = plt.subplots(nrows=2, ncols=2)

corrEtotplot(axes.flat[0],Eact,Ecmp,0,132,'ANI-1')
corrEtotplot(axes.flat[1],Eact,Ecmp,46,55,'ANI-1')
print(Emin[46:55])
print(Efle[46:55])

corrEtotplot(axes.flat[2],Eact,Ecmp,50,51,'ANI-1')

#a=50
#b=51
#for i in range(len(Emin[a:b])):
#    axes.flat[0].annotate(Efle[a:b][i], xy=(Emin[a:b][i], Emin[a:b][i]), xytext=(Emin[a:b][i]+0.10, Emin[a:b][i]-0.3),fontsize=20)
#    axes.flat[1].annotate(Efle[a:b][i], xy=(Emin[a:b][i], Emin[a:b][i]), xytext=(Emin[a:b][i]+0.10, Emin[a:b][i]-0.3),fontsize=20)

#plt.plot (Eact, Eact, color='black',  label='DFT', linewidth=1)
#plt.scatter (Eact, Ecmp, marker=r'x', color='red',  label='TGM RMSE: ' + "{:.3f}".format(rmse1),  linewidth=1)

#plt.ylabel('$E_{cmp}$ (kcal/mol)')
#plt.xlabel('$E_{ref}$ (kcal/mol)')

#plt.title("Asolute E - GDB-10 test data correlation")
#plt.legend(bbox_to_anchor=(0.01, 0.99), loc=2, borderaxespad=0.,fontsize=14)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 12}

plt.rc('font', **font)

plt.show()