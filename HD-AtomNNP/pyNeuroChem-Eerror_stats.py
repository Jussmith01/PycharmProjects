__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
import numpy as np
import matplotlib.pyplot as plt
import time as tm
from scipy import stats as st
from os import listdir
import os.path

def setmaxE(X,Y,E):
    m = min(X)
    newlist = []
    for i in range(0,len(X)):
        if gt.hatokcal * abs(X[i]-m) <= E:
            newlist.append(Y[i])

    return newlist

def setminE(X,Y):
    m = min(X)
    newlist = []
    newlist.append(Y.min())

    return newlist

def MAE (act,pre):
    N = act.shape[0]
    e = (np.abs(pre-act)).sum()
    return e / float(N)

def MAPE (act,pre):
    N = act.shape[0]
    e = 100.0 * (np.abs((act-pre)/act)).sum()
    return e / float(N)

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

def corrEplot(ax,d1,d2,shr1,shr2):
    ax.plot(d1,d1,color='black',linewidth=2,)
    ax.scatter(d1, d2, color='red', marker=r"o",s=50)

    # Set Limits
    ax.set_xlim([shr1,shr2])
    ax.set_ylim([shr1,shr2])

    font = {'family': 'Bitstream Vera Sans',
            'weight': 'heavy',
            'size': 24}

    ax.set_ylabel('$E_{cmp}$',fontdict=font)
    ax.set_xlabel('$E_{ref}$',fontdict=font)


# Set data fields
#dtdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/testdata/'
#dtdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/peptidetestdata/'
dtdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_10/testdata/'
#fpref = 'gdb11_10-'
#fpost = '_test.dat'
#rng = [0,140]

files = listdir(dtdir)
#print (files)
# Set required files for pyNeuroChem

#Network 1 Files
#wkdir1    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dataset_size_testing/train08b_05p_4/'
wkdir1    = '/home/jujuman/Research/wB97X-631gd-train-highgarden/train_08-a3.1A_r4.6_AEV384_1/'

cnstfile1 = wkdir1 + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile1  = wkdir1 + '../sae_6-31gd.dat'
nnfdir1   = wkdir1 + 'networks/'

E_max = 200.0 # an energy cuttoff for error considerations in kcal/mol

# Construct pyNeuroChem classes
nc = pync.pyNeuroChem(cnstfile1, saefile1, nnfdir1, 1)

Ecmp = []
Eact = []
Eerr = []
time = 0.0

ld = [0,0.0]
sd = [0,100000.0]

err = []
sze = []

Herror = 0.0
Wfile = ''

Lerror = 100.0
Bfile = ''

Nf = len(files)
cnt = 0

mNa = 100

for i in files:
    cnt += 1
#for i in range(rng[0],rng[1]):
    #xyz,typ,Eact_t,readf    = gt.readncdat(dtdir + fpref + str(i) + fpost)
    xyz,typ,Eact_W,readf    = gt.readncdat(dtdir + i)

    xyz = np.asarray(xyz,dtype=np.float32)
    xyz = xyz.reshape((xyz.shape[0],len(typ),3))

    if readf and xyz.shape[0] > 0:

        #print('FILE: ' + dtdir + fpref + str(i) + fpost)
        print('FILE: ' + str(cnt) + ' of ' + str(Nf) + ' ' + i)

        Nm = xyz.shape[0]
        Na = xyz.shape[1]

        if Na < mNa:
            mNa = Na

        Nat = Na * Nm

        Nit = int(np.ceil(Nat/65000.0))
        Nmo = int(65000/Na)
        Nmx = Nm

        for j in range(0,Nit):
            #print("Index: " + str(j*Nmo) + " " + str(min(j*Nmo+Nmo,Nm)) )

            # Setup idicies
            i1 = j*Nmo
            i2 = min(j*Nmo+Nmo,Nm)

            #copy array subset
            Eact_t = Eact_W[i1:i2]

            # Set the conformers in NeuroChem
            nc.setConformers(confs=xyz[i1:i2],types=typ)

            # Print some data from the NeuroChem
            #print( '1) Number of Atoms Loaded: ' + str(nc.getNumAtoms()) )
            #print( '1) Number of Confs Loaded: ' + str(nc.getNumConfs()) )

            # Compute Energies of Conformations
            #print('Computing energies...')
            _t1b = tm.time()
            Ecmp_t = nc.energy()
            _t2b = (tm.time() - _t1b) * 1000.0
            #print('Computation complete. Time: ' + "{:.4f}".format(_t2b)  + 'ms')

            Ecmp_t = Ecmp_t - Ecmp_t.min()
            Eact_t = Eact_t - Eact_t.min()

            Ecmp_t = setmaxE(Eact_t, Ecmp_t, E_max)
            Eact_t = setmaxE(Eact_t, Eact_t, E_max)

            #Ecmp_t = setminE(Eact_t, Ecmp_t)
            #Eact_t = setminE(Eact_t, Eact_t)

            deltas = gt.hatokcal * np.abs(Ecmp_t - np.array(Eact_t, dtype=float))
            Me = max (deltas)
            if Me > Herror:
                Herror = Me
                Wfile = i

            Le = min (deltas)
            if Le < Lerror:
                Lerror = Le
                Bfile = i

            print (gt.hatokcal * gt.calculaterootmeansqrerror(np.array(Eact_t, dtype=float),Ecmp_t))

            tNa = nc.getNumAtoms()
            err.append(gt.hatokcal * gt.calculaterootmeansqrerror(np.array(Eact_t, dtype=float),Ecmp_t) / float(tNa))
            sze.append(float(len(Eact_t)))

            time += _t2b

            Ecmp += Ecmp_t
            Eact += Eact_t

#plt_by_index(np.array(Eerr),-1)

Ndps = len(Ecmp)

print ('\nMax Delta (kcal/mol): ' + str(Herror) + ' FILE: ' + Wfile)
print ('Min Delta (kcal/mol): ' + str(Lerror) + ' FILE: ' + Bfile)
print('\nMAXE')
print(ld)
print('MINE')
print(sd)

print('Min Na: ' + str(mNa))

Ecmp = gt.hatokcal * np.array(Ecmp, dtype=float)
Eact = gt.hatokcal * np.array(Eact, dtype=float)

err = np.array(err, dtype=float)
sze = np.array(sze, dtype=float)

rmse1 = gt.calculaterootmeansqrerror(Eact,Ecmp)
mae1 = MAE(Eact,Ecmp)
slope1, intercept1, r_value1, p_value1, std_err1 = st.linregress(Eact,Ecmp)

mx = Eact.max()
mn = Eact.min()

print('Max: ' + str(mx) + ' Min: ' + str(mn))

print ("{:.7f}".format(mae1))
print ("{:.7f}".format(100.0*(mae1/(mx - mn))))
print ("{:.7f}".format(rmse1))
print ("{:.7f}".format(100.0*(rmse1/(mx - mn))))
print ("{:.7f}".format(MAPE(Eact,Ecmp)))
print ('Per atom RMSE (kcal/mol) = ' + "{:.7f}".format( ((sze / float(Ndps)) * err).sum() ))
print ("{:.7f}".format(slope1))
print ("{:.7f}".format(intercept1))
print ("{:.7f}".format(r_value1))
print ("{:.7f}".format(r_value1**2))
print ("{:.7f}".format(time))
print (str(Ndps))

#slope1, intercept1, r_value1, p_value1, std_err1 = st.linregress(Eact,Ecmp)

#fig, axes = plt.subplots(nrows=2, ncols=2)

#corrEplot(axes.flat[0],Eact,Ecmp,Eact.min(),Eact.max())
#corrEplot(axes.flat[1],Eact,Ecmp,Eact.min(),Eact.max())
#corrEplot(axes.flat[2],Eact,Ecmp,Eact.min(),Eact.max())
#corrEplot(axes.flat[3],Eact,Ecmp,Eact.min(),Eact.max())

#plt.plot (Eact, Eact, color='black',  label='DFT', linewidth=1)
#plt.scatter (Eact, Ecmp, marker=r'x', color='red',  label='TGM RMSE: ' + "{:.3f}".format(rmse1),  linewidth=1)

#plt.ylabel('$E_{cmp}$ (kcal/mol)')
#plt.xlabel('$E_{ref}$ (kcal/mol)')

#plt.title("Asolute E - GDB-10 test data correlation")
#plt.legend(bbox_to_anchor=(0.01, 0.99), loc=2, borderaxespad=0.,fontsize=14)

#font = {'family' : 'Bitstream Vera Sans',
#        'weight' : 'normal',
#        'size'   : 18}

#plt.rc('font', **font)

plt.show()
