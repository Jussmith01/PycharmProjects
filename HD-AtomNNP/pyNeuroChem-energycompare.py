__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
import numpy as np
import matplotlib.pyplot as plt
import time as tm
#from scipy import stats as st

def sortbyother(Y, X):
    xy = zip(X, Y)
    xy.sort()
    X, Y = zip(*xy)
    return np.array(Y)

# Set required files for pyNeuroChem
wkdir1    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/train_08_3/'
wkdir2    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/train_08_5/'
#wkdir2    = '/home/jujuman/Research/trainingcases/wB97X-631gd-train-highgarden/nw-64-64-64-32/train_07-a3.1A_r4.5/'

#Network 1 Files
cnstfile1 = wkdir1 + 'rHCNO-4.6A_32-3.1A_a8-8.params'
saefile1  = wkdir1 + 'sae_6-31gd.dat'
nnfdir1   = wkdir1 + 'networks/'

# Network 2 Files
cnstfile2 = wkdir2 + 'rHCNO-4.7A_32-3.2A_a8-8.params'
saefile2  = wkdir2 + 'sae_6-31gd.dat'
nnfdir2   = wkdir2 + 'networks/'

<<<<<<< HEAD
# Construct pyNeuroChem classes
nc1 = pync.pyNeuroChem(cnstfile1,saefile1,nnfdir1,0)
#nc2 = pync.pyNeuroChem(cnstfile2,saefile2,nnfdir2,1)

dtdir = '/home/jujuman/Dropbox/Research/ChemSciencePaper/TestCases/Bonds/Fentanyl/'
=======
dtdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_testdata/Atazanavir/data/'
>>>>>>> 353a9faba101b20434eaf1228a2883f24bcf077a

#xyz,typ,Eact = gt.readncdat('../data_irc.dat')
xyz,typ,Eact    = gt.readncdat(dtdir + 'atazanavir_DFT2_test.dat')
xyz2,typ2,Eact2 = gt.readncdat(dtdir + 'atazanavir_AM1_test.dat')
xyz3,typ3,Eact3 = gt.readncdat(dtdir + 'atazanavir_PM6_test.dat')

Eact = np.array(Eact)
Eact2 = np.array(Eact2)
Eact3 = np.array(Eact3)

# Construct pyNeuroChem classes
nc1 = pync.pyNeuroChem(cnstfile1,saefile1,nnfdir1,1)

# Set the conformers in NeuroChem
nc1.setConformers(confs=xyz,types=typ)

# Print some data from the NeuroChem
print( '1) Number of Atoms Loaded: ' + str(nc1.getNumAtoms()) )
print( '1) Number of Confs Loaded: ' + str(nc1.getNumConfs()) )

# Compute Forces of Conformations
print('Computing energies 1...')
_t1b = tm.time()
Ecmp1 = np.array( nc1.computeEnergies() )
print('Computation complete 1. Time: ' + "{:.4f}".format((tm.time() - _t1b) * 1000.0)  + 'ms')

# Construct pyNeuroChem classes
nc2 = pync.pyNeuroChem(cnstfile2,saefile2,nnfdir2,1)

# Set the conformers in NeuroChem
nc2.setConformers(confs=xyz,types=typ)

# Print some data from the NeuroChem
print( '2) Number of Atoms Loaded: ' + str(nc2.getNumAtoms()) )
print( '2) Number of Confs Loaded: ' + str(nc2.getNumConfs()) )

# Compute Forces of Conformations
print('Computing energies 1...')
_t2b = tm.time()
Ecmp2 = np.array( nc2.computeEnergies() )
print('Computation complete 1. Time: ' + "{:.4f}".format((tm.time() - _t2b) * 1000.0)  + 'ms')

Ecmp1 = gt.hatokcal * Ecmp1
Ecmp2 = gt.hatokcal * Ecmp2
Eact  = gt.hatokcal * Eact
Eact2 = gt.hatokcal * Eact2
Eact3 = gt.hatokcal * Eact3

IDX = np.arange(0,Eact.shape[0],1,dtype=float)
#IDX2 = sortbyother(IDX, Eact)
#print IDX2

Ecmp1 = sortbyother(Ecmp1, Eact)
Ecmp2 = sortbyother(Ecmp2, Eact)
Eact2 = sortbyother(Eact2, Eact)
Eact3 = sortbyother(Eact3, Eact)
Eact = np.sort( Eact )

#print (Eact2.min())

Ecmp1 = Ecmp1 - Ecmp1.min()
Ecmp2 = Ecmp2 - Ecmp2.min()
Eact  = Eact  - Eact.min()
Eact2 = Eact2 - Eact2.min()
Eact3 = Eact3 - Eact3.min()

rmse1 = gt.calculaterootmeansqrerror(Eact,Ecmp1)
rmse2 = gt.calculaterootmeansqrerror(Eact,Ecmp2)
rmse3 = gt.calculaterootmeansqrerror(Eact,Eact2)
rmse4 = gt.calculaterootmeansqrerror(Eact,Eact3)

<<<<<<< HEAD
plt.plot   (IDX, Eact, color='black',  label='DFT',linewidth=3)
plt.scatter(IDX, Eact, marker='o' , color='black',  linewidth=3)

#print ( "Spearman corr. 1: " + "{:.3f}".format(st.spearmanr(Ecmp1,Eact)[0]) )
#print ( "Spearman corr. 1: " + "{:.3f}".format(st.spearmanr(Eact2,Eact)[0]) )
#print ( "Spearman corr. 1: " + "{:.3f}".format(st.spearmanr(Eact3,Eact)[0]) )
=======
print ( "Spearman corr. 1: " + "{:.7f}".format(st.spearmanr(Ecmp1,Eact)[0]) )
print ( "Spearman corr. 2: " + "{:.7f}".format(st.spearmanr(Ecmp2,Eact)[0]) )
print ( "Spearman corr. 3: " + "{:.7f}".format(st.spearmanr(Eact2,Eact)[0]) )
print ( "Spearman corr. 4: " + "{:.7f}".format(st.spearmanr(Eact3,Eact)[0]) )
>>>>>>> 353a9faba101b20434eaf1228a2883f24bcf077a

plt.plot   (IDX, Ecmp1, ':', marker=r'D', color='red',  label='ANN - c08b RMSE: ' + "{:.3f}".format(rmse1) + ' kcal/mol',  linewidth=2, markersize=8)
plt.plot   (IDX, Ecmp2, ':', marker=r'D', color='orange',  label='ANN - c08c RMSE: ' + "{:.3f}".format(rmse2) + ' kcal/mol',  linewidth=2, markersize=8)
plt.plot   (IDX, Eact2, ':', marker=r'v', color='blue', label='AM1   RMSE: ' + "{:.2f}".format(rmse3) + ' kcal/mol',  linewidth=2, markersize=8)
plt.plot   (IDX, Eact3, ':', marker=r'*', color='green',label='PM6   RMSE: ' + "{:.2f}".format(rmse4) + ' kcal/mol',  linewidth=2, markersize=10)

plt.plot   (IDX, Eact, color='black',  label='DFT',linewidth=3)
plt.scatter(IDX, Eact, marker='o' , color='black',  linewidth=4)

plt.title("300K NMS structures of\nNME-Gly-Pro-Hyp-Gly-Ala-Gly-ACE")

plt.ylabel('E calculated (kcal/mol)')
plt.xlabel('Structure')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.,fontsize=14)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

plt.show()
