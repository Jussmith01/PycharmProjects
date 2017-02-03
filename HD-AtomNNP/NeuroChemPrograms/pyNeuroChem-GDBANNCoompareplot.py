__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
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

def pyNCcomputeTestSet(cnstfile1,saefile1,nnfdir1,dtdir,dtdftpref,dtpm6dir,dtpm6pref,N,P=1.0):
    # Construct pyNeuroChem classes
    nc = pync.pyNeuroChem(cnstfile1,saefile1,nnfdir1,0)

    Eact = []
    Ecmp = []
    Eotr = []
    Ndat = 0
    Nmol = 0
    t = 0.0
    for i in range(0,N):
        print ('|----- File ' + str(i) + ' -----|')
        print ('Name: ' + dtdir + dtdftpref + str(i) + '_test.dat')

        rv = bool(np.random.binomial(1,P))
        if (os.path.isfile(dtdir + dtdftpref + str(i) + '_test.dat') and rv):

            xyz,typ,Eact_t    = gt.readncdat(dtdir + dtdftpref + str(i) + '_test.dat')
            xyz1,typ1,Eotr_t  = gt.readncdat(dtpm6dir + dtpm6pref + str(i) + '_test.dat')

            if len(Eact_t) == len(Eotr_t):

                Eact += shiftlsttomin( Eact_t )
                Eotr += shiftlsttomin( Eotr_t )

                #Eact +=  Eact_t
                #Eotr +=  Eotr_t

                Nmol += 1
                Ndat += len( Eact_t )

                # Set the conformers in NeuroChem
                nc.setConformers(confs=xyz,types=typ)

                # Compute Forces of Conformations
                print('Computing energies 1...')
                _t1b = tm.time()
                Ecmp_t = nc.computeEnergies()
                _t2b = (tm.time() - _t1b) * 1000.0
                t += _t2b
                print('Computation complete 1. Time: ' + "{:.4f}".format(_t2b)  + 'ms')
                Ecmp += shiftlsttomin(  Ecmp_t )
                #Ecmp +=  Ecmp_t

    Eact = np.array(Eact,dtype=float)
    Eotr = np.array(Eotr,dtype=float)
    Ecmp = np.array(Ecmp,dtype=float)

    return Eact,Ecmp,Eotr,Ndat,Nmol,t

# Set data fields
dtdir =  '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_10/testdata/'
dtdftpref = 'gdb11_s10-'

dtpm6dir =  '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_10/DFTBdata/'
dtpm6pref = 'gdb11_s10-'

# Set required files for pyNeuroChem

#Network 1 Files
wkdir1    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/train_08_3/'
cnstfile1 = wkdir1 + 'rHCNO-4.6A_32-3.1A_a8-8.params'
saefile1  = wkdir1 + '../sae_6-31gd.dat'
nnfdir1   = wkdir1 + 'networks/'

# Network 2 Files
wkdir2    = '/home/jujuman/Research/trainingcases/wB97X-631gd-train-highgarden/nw-64-64-64-32/train_07-a3.1A_r4.5/'
cnstfile2 = wkdir2 + 'rHCNO-4.5A_32-3.1A_a8-8.params'
saefile2  = wkdir2 + '../../sae_6-31gd.dat'
nnfdir2   = wkdir2 + 'networks/'

# Network 3 Files
wkdir3    = '/home/jujuman/Research/trainingcases/wB97X-631gd-train-highgarden/train_06/'
cnstfile3 = wkdir3 + 'rHCNO-4.3A_32-3.0A_a8-8.params'
saefile3  = wkdir3 + '../sae_6-31gd.dat'
nnfdir3   = wkdir3 + 'networks/'

# Network 3 Files
wkdir4    = '/home/jujuman/Research/trainingcases/wB97X-631gd-train-highgarden/train_05/'
cnstfile4 = wkdir4 + 'rHCNO-4.2A_32-2.9A_a8-8.params'
saefile4  = wkdir4 + '../sae_6-31gd.dat'
nnfdir4   = wkdir4 + 'networks/'

# Network 3 Files
wkdir5    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/train_04/'
cnstfile5 = wkdir5 + 'rHCNO-3.5A_32-2.8A_a8-8.params'
saefile5  = wkdir5 + 'sae_6-31gd.dat'
nnfdir5   = wkdir5 + 'networks/'

N = 120

Eact,Ecmp1,Eotr,Ndat,Nmol,t1 = pyNCcomputeTestSet(cnstfile1,saefile1,nnfdir1,dtdir,dtdftpref,dtpm6dir,dtpm6pref,N)
Eact,Ecmp2,Eotr,Ndat,Nmol,t2 = pyNCcomputeTestSet(cnstfile2,saefile2,nnfdir2,dtdir,dtdftpref,dtpm6dir,dtpm6pref,N)
Eact,Ecmp3,Eotr,Ndat,Nmol,t3 = pyNCcomputeTestSet(cnstfile3,saefile3,nnfdir3,dtdir,dtdftpref,dtpm6dir,dtpm6pref,N)
Eact,Ecmp4,Eotr,Ndat,Nmol,t4 = pyNCcomputeTestSet(cnstfile4,saefile4,nnfdir4,dtdir,dtdftpref,dtpm6dir,dtpm6pref,N)
Eact,Ecmp5,Eotr,Ndat,Nmol,t5 = pyNCcomputeTestSet(cnstfile5,saefile5,nnfdir5,dtdir,dtdftpref,dtpm6dir,dtpm6pref,N)

print ("Time 1: " + "{:.4f}".format(t1/1000.0)  + 's')
print ("Time 2: " + "{:.4f}".format(t2/1000.0)  + 's')
print ("Time 3: " + "{:.4f}".format(t3/1000.0)  + 's')
print ("Time 4: " + "{:.4f}".format(t4/1000.0)  + 's')
print ("Time 5: " + "{:.4f}".format(t5/1000.0)  + 's')

Ecmp1 = gt.hatokcal * Ecmp1
Ecmp2 = gt.hatokcal * Ecmp2
Ecmp3 = gt.hatokcal * Ecmp3
Ecmp4 = gt.hatokcal * Ecmp4
Ecmp5 = gt.hatokcal * Ecmp5
Eact  = gt.hatokcal * Eact
Eotr  = gt.hatokcal * Eotr

Emax = 400.0
Ecmp1 = setmaxE(Eact,Ecmp1,Emax)
Ecmp2 = setmaxE(Eact,Ecmp2,Emax)
Ecmp3 = setmaxE(Eact,Ecmp3,Emax)
Ecmp4 = setmaxE(Eact,Ecmp4,Emax)
Ecmp5 = setmaxE(Eact,Ecmp5,Emax)
Eotr =  setmaxE(Eact,Eotr, Emax)
Eact =  setmaxE(Eact,Eact, Emax)

rmse1 = gt.calculaterootmeansqrerror(Eact,Eotr)
rmse2 = gt.calculaterootmeansqrerror(Eact,Ecmp1)
rmse3 = gt.calculaterootmeansqrerror(Eact,Ecmp2)
rmse4 = gt.calculaterootmeansqrerror(Eact,Ecmp3)
rmse5 = gt.calculaterootmeansqrerror(Eact,Ecmp4)
rmse6 = gt.calculaterootmeansqrerror(Eact,Ecmp5)

mx = Eact.max()
mn = Eact.min()

#plt.scatter(IDX, Eact, marker='o' , color='black',  linewidth=3)

print ( "Spearman corr. DFTB: " + "{:.3f}".format(st.spearmanr(Eotr,Eact)[0]) )
print ( "Spearman corr. TGM 08: " + "{:.3f}".format(st.spearmanr(Ecmp1,Eact)[0]) )
#print ( "Spearman corr. TGM 07: " + "{:.3f}".format(st.spearmanr(Ecmp2,Eact)[0]) )
#print ( "Spearman corr. TGM 06: " + "{:.3f}".format(st.spearmanr(Ecmp3,Eact)[0]) )
#print ( "Spearman corr. TGM 05: " + "{:.3f}".format(st.spearmanr(Ecmp4,Eact)[0]) )

slope1, intercept1, r_value1, p_value1, std_err1 = st.linregress(Eact,Eotr)
slope2, intercept2, r_value2, p_value2, std_err2 = st.linregress(Eact,Ecmp1)
#slope3, intercept3, r_value3, p_value3, std_err3 = st.linregress(Eact,Ecmp2)
#slope4, intercept4, r_value4, p_value4, std_err4 = st.linregress(Eact,Ecmp3)
#slope5, intercept5, r_value5, p_value5, std_err5 = st.linregress(Eact,Ecmp4)

fig = plt.figure()

ax1 = fig.add_subplot(221)
ax1.plot   ((mn,mx), (mn,mx), color='blue',  label='DFT',linewidth=3)
ax1.scatter (Eact, Eotr, marker=r'o', color='#7BAFD4',  label='DFTB    RMSE: ' + "{:.3f}".format(rmse1) + ' kcal/mol slope: ' + "{:.3f}".format(slope1) + " intercept: " + "{:.3f}".format(intercept1) + " $r^2$: " + "{:.3f}".format(r_value1**2),  linewidth=1)

ax1.set_title("GDB-10 test data correlation (" + str(Ndat) + " data points from " + str(Nmol) + " molecules)")
ax1.set_ylabel('$E_{cmp}$ (kcal/mol)')
ax1.set_xlabel('$E_{ref}$ (kcal/mol)')
ax1.legend(bbox_to_anchor=(0.01, 0.99), loc=2, borderaxespad=0.,fontsize=14)

ax1 = fig.add_subplot(222)
ax1.plot   ((mn,mx), (mn,mx), color='blue',  label='DFT',linewidth=3)
ax1.scatter (Eact, Ecmp1, marker=r'x', color='orange',  label='TGM08 RMSE: ' + "{:.3f}".format(rmse2) + ' kcal/mol slope: ' + "{:.3f}".format(slope2) + " intercept: " + "{:.3f}".format(intercept2) + " $r^2$: " + "{:.3f}".format(r_value2**2),  linewidth=1)

ax1.set_title("GDB-10 test data correlation (" + str(Ndat) + " data points from " + str(Nmol) + " molecules)")
ax1.set_ylabel('$E_{cmp}$ (kcal/mol)')
ax1.set_xlabel('$E_{ref}$ (kcal/mol)')
ax1.legend(bbox_to_anchor=(0.01, 0.99), loc=2, borderaxespad=0.,fontsize=14)

ax1 = fig.add_subplot(223)
ax1.plot   ((mn,mx), (mn,mx), color='blue',  label='DFT',linewidth=3)
ax1.scatter (Eact, Ecmp2, marker=r'x', color='orange',  label='TGM08 RMSE: ' + "{:.3f}".format(rmse2) + ' kcal/mol slope: ' + "{:.3f}".format(slope2) + " intercept: " + "{:.3f}".format(intercept2) + " $r^2$: " + "{:.3f}".format(r_value2**2),  linewidth=1)

ax1.set_title("GDB-10 test data correlation (" + str(Ndat) + " data points from " + str(Nmol) + " molecules)")
ax1.set_ylabel('$E_{cmp}$ (kcal/mol)')
ax1.set_xlabel('$E_{ref}$ (kcal/mol)')
ax1.legend(bbox_to_anchor=(0.01, 0.99), loc=2, borderaxespad=0.,fontsize=14)

ax1 = fig.add_subplot(224)
ax1.plot   ((mn,mx), (mn,mx), color='blue',  label='DFT',linewidth=3)
ax1.scatter (Eact, Ecmp3, marker=r'x', color='orange',  label='TGM08 RMSE: ' + "{:.3f}".format(rmse2) + ' kcal/mol slope: ' + "{:.3f}".format(slope2) + " intercept: " + "{:.3f}".format(intercept2) + " $r^2$: " + "{:.3f}".format(r_value2**2),  linewidth=1)

ax1.set_title("GDB-10 test data correlation (" + str(Ndat) + " data points from " + str(Nmol) + " molecules)")
ax1.set_ylabel('$E_{cmp}$ (kcal/mol)')
ax1.set_xlabel('$E_{ref}$ (kcal/mol)')
ax1.legend(bbox_to_anchor=(0.01, 0.99), loc=2, borderaxespad=0.,fontsize=14)

#ax2 = fig.add_subplot(122)
#ar1 = [int(4),int(5),int(6),int(7),int(8)]
#ar2 = [rmse6,rmse5,rmse4,rmse3,rmse2]
#ax2.plot (ar1, ar2,marker=r'o', color='blue',linewidth=4,markersize=10)

#ax2.set_title("GDB-10 test data error computed on different TGM networks")
#ax2.set_ylabel('RMSE (kcal/mol)')
#ax2.set_xlabel('Max GDB set size')

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

plt.show()
