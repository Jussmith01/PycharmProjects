__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import hdnntools as gt
import numpy as np
import matplotlib.pyplot as plt
import time as tm
from scipy import stats as st
import os.path

def shiftlsttomin(X):
    m = min(X)
    for i in range(0,len(X)):
        X[i] = X[i] - m

    return X

def setmaxE(X,Y,E):
    m = min(X)
    newlist = []
    for i in range(0,len(X)):
        if X[i] <= E:
            newlist.append(Y[i])

    return np.array(newlist,dtype=float)

def sortbyother(Y, X):
    xy = zip(X, Y)
    xy.sort()
    X, Y = zip(*xy)
    return np.array(Y)

def pyNCcomputeTestSet(cnstfile1,saefile1,nnfdir1,dtdir,P=1.0):
    # Construct pyNeuroChem classes
    nc = pync.conformers(cnstfile1,saefile1,nnfdir1,0)

    files = os.listdir(dtdir)

    Eact = []
    Ecmp = []

    Ndat = 0
    Nmol = 0
    t = 0.0

    for i in files:
        print ('|----- File ' + str(i) + ' -----|')
        print ('Name: ' + dtdir + i)

        rv = bool(np.random.binomial(1,P))
        if rv:

            xyz,species,energy = gt.readncdat(dtdir + i, np.float32)

            shiftlsttomin(energy)
            Eact.extend(energy)

            Nmol += 1
            Ndat += energy.shape[0]

            # Set the conformers in NeuroChem
            nc.setConformers(confs=xyz,types=list(species))

            # Compute Forces of Conformations
            print('Computing energies...')
            _t1b = tm.time()
            energy_comp = nc.energy()
            _t2b = (tm.time() - _t1b) * 1000.0
            t += _t2b
            print('Computation complete. Time: ' + "{:.4f}".format(_t2b)  + 'ms')

            shiftlsttomin(energy_comp)
            Ecmp.extend(energy_comp)

#            for j,k in zip(energy,energy_comp):
#                print('  ',j, ':', k)

    Eact = np.concatenate([Eact])
    Ecmp = np.concatenate([Ecmp])

    return Eact,Ecmp,Ndat,Nmol,t

# Set data fields
dtdir =  '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_10/testdata/'
dtdftpref = 'gdb11_s10-'

# Set required files for pyNeuroChem
#Network 1 Files
wkdir1    = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-1-ntwk/'
cnstfile1 = wkdir1 + 'rHCNO-4.6A_32-3.1A_a8-8.params'
saefile1  = wkdir1 + 'sae_6-31gd.dat'
nnfdir1   = wkdir1 + 'networks/'

Eact,Ecmp1,Ndat,Nmol,t1 = pyNCcomputeTestSet(cnstfile1,saefile1,nnfdir1,dtdir)

print(Eact)

for j, k in zip(Eact, Ecmp1):
    print('  ', j, ':', k)

print ("Time 1: " + "{:.4f}".format(t1/1000.0)  + 's')

Ecmp1 = gt.hatokcal * Ecmp1
Eact  = gt.hatokcal * Eact

Emax = 100.0
Ecmp1 = setmaxE(Eact,Ecmp1,Emax)
Eact =  setmaxE(Eact,Eact, Emax)

print('NMOL: ',Ecmp1.shape[0])

rmse2 = gt.calculaterootmeansqrerror(Eact,Ecmp1)

mx = Eact.max()
mn = Eact.min()

#<<<<<<< HEAD
#plt.scatter(IDX, Eact, marker='o' , color='black',  linewidth=3)

#=======
#print ( "Spearman corr. DFTB: " + "{:.3f}".format(st.spearmanr(Eotr,Eact)[0]) )
#>>>>>>> 1ec96245cdd1c6d1647ffeed1a6bec3a4b7e4bb4
print ( "Spearman corr. TGM 08: " + "{:.3f}".format(st.spearmanr(Ecmp1,Eact)[0]) )

slope2, intercept2, r_value2, p_value2, std_err2 = st.linregress(Eact,Ecmp1)

fig = plt.figure()

ax1 = fig.add_subplot(111)
ax1.plot   ((mn,mx), (mn,mx), color='blue',  label='DFT',linewidth=3)
#ax1.scatter (Eact, Ecmp1, marker=r'o', color='#7BAFD4',  label='ANI-1    RMSE: ' + "{:.1f}".format(rmse2) + ' kcal/mol',  linewidth=1)
ax1.scatter (Eact, Ecmp1, marker=r'o', color='#7BAFD4',  label='ANI-1',  linewidth=1)

#ax1.set_title("GDB-10 test data correlation (" + str(Ndat) + " data points from " + str(Nmol) + " molecules)")
ax1.set_title("Extensibility test set")
ax1.set_ylabel('$E_{cmp}$ (kcal/mol)')
ax1.set_xlabel('$E_{ref}$ (kcal/mol)')
ax1.legend(bbox_to_anchor=(0.01, 0.99), loc=2, borderaxespad=0.,fontsize=20)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 24}

plt.rc('font', **font)

plt.show()
