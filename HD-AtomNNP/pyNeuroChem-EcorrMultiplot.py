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
        if abs(X[i]) <= E:
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

        rv = bool(np.random.binomial(1,P))
        if (os.path.isfile(dtdir + dtdftpref + str(i) + '_test.dat') and rv):

            xyz,typ,Eact_t,chk    = gt.readncdat(dtdir + dtdftpref + str(i) + '_test.dat')
            xyz1,typ1,Eotr_t,chk  = gt.readncdat(dtpm6dir + dtpm6pref + str(i) + '_test.dat')

            if len(Eact_t) == len(Eotr_t):

                #print ('|----- File ' + str(i) + ' -----|')
                #print ('Name: ' + dtdir + dtdftpref + str(i) + '_test.dat')

                Eact += shiftlsttomin( Eact_t )
                Eotr += shiftlsttomin( Eotr_t )

                #Eact +=  Eact_t
                #Eotr +=  Eotr_t

                Nmol += 1
                Ndat += len( Eact_t )

                # Set the conformers in NeuroChem
                nc.setConformers(confs=xyz,types=typ)

                # Compute Forces of Conformations
                print(' ' + str(Nmol) + ') Computing ' + str(len( Eact_t )) + ' energies...')
                _t1b = tm.time()
                Ecmp_t = nc.computeEnergies()
                _t2b = (tm.time() - _t1b) * 1000.0
                t += _t2b
                #print('Computation complete. Time: ' + "{:.4f}".format(_t2b)  + 'ms')
                Ecmp += shiftlsttomin(  Ecmp_t )
                #Ecmp +=  Ecmp_t99
            else:
                print (str(len(Eact_t)) + '!=' + str(len(Eotr_t)) + ' File: ' + dtdir + dtdftpref + str(i) + '_test.dat')
        else:
            print('File not found: ' + dtdir + dtdftpref + str(i) + '_test.dat')
    Eact = np.array(Eact,dtype=float)
    Eotr = np.array(Eotr,dtype=float)
    Ecmp = np.array(Ecmp,dtype=float)

    return Eact,Ecmp,Eotr,Ndat,Nmol,t

def Ecorrplot (ax1,Eact,Ecmp,mlbl,color,lab=False):
    mx = Eact.max()
    mn = Eact.min()

    if lab:
        ax1.plot((mn, mx), (mn, mx), color='black',label='DFT', linewidth=5)
    else:
        ax1.plot((mn, mx), (mn, mx), color='black', linewidth=5)

    rmse = gt.calculaterootmeansqrerror(Eact, Ecmp)
    ax1.scatter(Eact, Ecmp, marker=r'o', color=color, label= mlbl + ' RMSE: ' + "{:.3f}".format(rmse) + ' kcal/mol', linewidth=1)

    ax1.set_xlim([mn,mx])
    ax1.set_ylim([mn,mx])

    #ax1.set_title("title)
    ax1.set_ylabel('$\Delta E_{cmp}$ (kcal/mol)')
    ax1.set_xlabel('$\Delta E_{ref}$ (kcal/mol)')
    ax1.legend(bbox_to_anchor=(0.01, 0.99), loc=2, borderaxespad=0., fontsize=16)


# Set data fields
dtdir =  '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_10/testdata/'
dtdftpref = 'gdb11_s10-'

dtdftbdir =  '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_10/DFTBdata/'
dtpm6dir  =  '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_10/PM6data/'
dtam1dir  =  '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_10/AM1data/'

# Set required files for pyNeuroChem

#Network 1 Files
wkdir1    = '/home/jujuman/Scratch/Research/trainingcases/wB97X-631gd-train-highgarden/train_08-a3.1A_r4.6_dn2/'
cnstfile1 = wkdir1 + 'rHCNO-4.6A_32-3.1A_a8-8.params'
saefile1  = wkdir1 + '../sae_6-31gd.dat'
nnfdir1   = wkdir1 + 'networks/'

# Network 2 Files
wkdir2    = '/home/jujuman/Research/trainingcases/wB97X-631gd-train-highgarden/nw-64-64-64-32/train_07-a3.1A_r4.5/'
cnstfile2 = wkdir2 + 'rHCNO-4.5A_32-3.1A_a8-8.params'
saefile2  = wkdir2 + '../../sae_6-31gd.dat'
nnfdir2   = wkdir2 + 'networks/'

# Network 3 Files
wkdir3    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/train_06/'
cnstfile3 = wkdir3 + 'rHCNO-4.5A_32-3.1A_a8-8.params'
saefile3  = wkdir3 + '../sae_6-31gd.dat'
nnfdir3   = wkdir3 + 'networks/'

# Network 3 Files
wkdir4    = '/home/jujuman/Research/trainingcases/wB97X-631gd-train-highgarden/train_05/'
cnstfile4 = wkdir4 + 'rHCNO-4.2A_32-2.9A_a8-8.params'
saefile4  = wkdir4 + '../sae_6-31gd.dat'
nnfdir4   = wkdir4 + 'networks/'

# Network 3 Files
wkdir5    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/train_04/'
cnstfile5 = wkdir5 + 'rHCNO-4.5A_32-3.1A_a8-8.params'
saefile5  = wkdir5 + 'sae_6-31gd.dat'
nnfdir5   = wkdir5 + 'networks/'

N = 140
Eact,Ecmp1,Eotr1,Ndat1,Nmol1,t1 = pyNCcomputeTestSet(cnstfile1,saefile1,nnfdir1,dtdir,dtdftpref,dtdftbdir,dtdftpref,N)
Eact,Ecmp2,Eotr1,Ndat1,Nmol1,t1 = pyNCcomputeTestSet(cnstfile2,saefile2,nnfdir2,dtdir,dtdftpref,dtdftbdir,dtdftpref,N)
Eact,Ecmp3,Eotr1,Ndat1,Nmol1,t1 = pyNCcomputeTestSet(cnstfile3,saefile3,nnfdir3,dtdir,dtdftpref,dtdftbdir,dtdftpref,N)
Eact,Ecmp4,Eotr1,Ndat1,Nmol1,t1 = pyNCcomputeTestSet(cnstfile4,saefile4,nnfdir4,dtdir,dtdftpref,dtdftbdir,dtdftpref,N)
Eact,Ecmp5,Eotr1,Ndat1,Nmol1,t1 = pyNCcomputeTestSet(cnstfile5,saefile5,nnfdir5,dtdir,dtdftpref,dtdftbdir,dtdftpref,N)

Eact,Ecmp1,Eotr2,Ndat2,Nmol2,t1 = pyNCcomputeTestSet(cnstfile1,saefile1,nnfdir1,dtdir,dtdftpref,dtpm6dir, dtdftpref,N)
Eact,Ecmp1,Eotr3,Ndat3,Nmol3,t1 = pyNCcomputeTestSet(cnstfile1,saefile1,nnfdir1,dtdir,dtdftpref,dtam1dir, dtdftpref,N)

print('Ndat1: ' + str(Nmol1) + ':' + str(Ndat1) + ' Ndat2: ' + str(Nmol2) + ':' + str(Ndat2) + ' Ndat3: ' + str(Nmol3) + ':' + str(Ndat3))

print ("Time 1: " + "{:.4f}".format(t1/1000.0)  + 's')

Ecmp1 = gt.hatokcal * Ecmp1
Ecmp2 = gt.hatokcal * Ecmp2
Ecmp3 = gt.hatokcal * Ecmp3
Ecmp4 = gt.hatokcal * Ecmp4
Ecmp5 = gt.hatokcal * Ecmp5
Eact  = gt.hatokcal * Eact
Eotr1  = gt.hatokcal * Eotr1
Eotr2  = gt.hatokcal * Eotr2
Eotr3  = gt.hatokcal * Eotr3

Emax = 300.0
Ecmp1 = setmaxE(Eact,Ecmp1,Emax)
Ecmp2 = setmaxE(Eact,Ecmp2,Emax)
Ecmp3 = setmaxE(Eact,Ecmp3,Emax)
Ecmp4 = setmaxE(Eact,Ecmp4,Emax)
Ecmp5 = setmaxE(Eact,Ecmp5,Emax)
Eotr1 = setmaxE(Eact,Eotr1,Emax)
Eotr2 = setmaxE(Eact,Eotr2,Emax)
Eotr3 = setmaxE(Eact,Eotr3,Emax)
Eact  = setmaxE(Eact,Eact, Emax)

print('Act count: ' + str(Eact.shape[0]))

rmse1 = gt.calculaterootmeansqrerror(Eact,Eotr1)
rmse2 = gt.calculaterootmeansqrerror(Eact,Eotr2)
rmse3 = gt.calculaterootmeansqrerror(Eact,Eotr3)
rmse4 = gt.calculaterootmeansqrerror(Eact,Ecmp1)
rmse5 = gt.calculaterootmeansqrerror(Eact,Ecmp2)
rmse6 = gt.calculaterootmeansqrerror(Eact,Ecmp3)
rmse7 = gt.calculaterootmeansqrerror(Eact,Ecmp4)
rmse8 = gt.calculaterootmeansqrerror(Eact,Ecmp5)

#plt.scatter(IDX, Eact, marker='o' , color='black',  linewidth=3)

print ( "Spearman corr. DFTB:  " + "{:.7f}".format(st.spearmanr(Eotr1,Eact)[0]) )
print ( "Spearman corr. PM6:   "  + "{:.7f}".format(st.spearmanr(Eotr2,Eact)[0]) )
print ( "Spearman corr. AM1:   "  + "{:.7f}".format(st.spearmanr(Eotr3,Eact)[0]) )
print ( "Spearman corr. ANI-1: "  + "{:.7f}".format(st.spearmanr(Ecmp1,Eact)[0]) )

slope1, intercept1, r_value1, p_value1, std_err1 = st.linregress(Eact,Eotr1)
slope2, intercept2, r_value2, p_value2, std_err2 = st.linregress(Eact,Eotr2)
slope3, intercept3, r_value3, p_value3, std_err3 = st.linregress(Eact,Eotr3)
slope4, intercept4, r_value4, p_value4, std_err4 = st.linregress(Eact,Ecmp1)

fig = plt.figure()

#fig.suptitle("Correlation plots")

ax = fig.add_subplot(231)
Ecorrplot(ax,Eact,Ecmp1,"ANI-1","red",True)

ax = fig.add_subplot(232)
Ecorrplot(ax,Eact,Eotr1,"DFTB","green")

ax = fig.add_subplot(234)
Ecorrplot(ax,Eact,Eotr2,"PM6","blue")

ax = fig.add_subplot(235)
Ecorrplot(ax,Eact,Eotr3,"AM1","orange")

ax = fig.add_subplot(133)
ar1 = [int(4),int(5),int(6),int(7),int(8)]
ar2 = [rmse8,rmse7,rmse6,rmse5,rmse4]
ax.plot (ar1, np.array(ar2),marker=r'o', color='blue',linewidth=4,markersize=15)

ax.annotate("{:.3f}".format(rmse8), xy=(ar1[0], ar2[0]), xytext=(ar1[0]+0.10, ar2[0]-0.3),fontsize=20)
ax.annotate("{:.3f}".format(rmse7), xy=(ar1[1], ar2[1]), xytext=(ar1[1]+0.10, ar2[1]-0.1),fontsize=20)
ax.annotate("{:.3f}".format(rmse6), xy=(ar1[2], ar2[2]), xytext=(ar1[2]+0.10, ar2[2]-0.1),fontsize=20)
ax.annotate("{:.3f}".format(rmse5), xy=(ar1[3], ar2[3]), xytext=(ar1[3]-0.75, ar2[3]-1.0),fontsize=20)
ax.annotate("{:.3f}".format(rmse4), xy=(ar1[4], ar2[4]), xytext=(ar1[4]-0.75, ar2[ 4]-1.0),fontsize=20)

#ax.set_xlim([min(ar2)-0.25, max(ar2)+0.25])

ax.set_xlabel('Maximum GDB training set size')
ax.set_ylabel('RMSE (kcal/mol)')

ar1 = [int(4),'',int(5),'',int(6),'',int(7),'',int(8)]
ax.set_xticklabels(ar1)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 20}

plt.rc('font', **font)

plt.show()
