from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import normalize
from sklearn.model_selection import train_test_split

from sklearn.neural_network import MLPRegressor

from sklearn.externals import joblib

from os import listdir
import time as tm
import random

import numpy as np
import graphtools as gt
import matplotlib.pyplot as plt

import coulomb_matrix as cm

def read_data_files_convert_cm(file,N):

    # Get training molecules
    xyz_tr, typ_tr, Eact_tr, readf = gt.readncdat(file, np.float32)

    # Compute energy of atoms at infinite separation
    ise = cm.computerISE(typ_tr)

    Eact_tr = (Eact_tr - ise)

    Nm = Eact_tr.shape[0]
    cmat = cm.GenCMatData2(xyz_tr,typ_tr,N)

    return cmat,Eact_tr,ise

def read_data_indexed(input_file,N,M,Nd,idx,idx_arr,cmats):

    linesize = 4*(N*N+1)

    c_idx = 0
    for i in idx_arr[idx*M:idx*M+M]:
        input_file.seek( i * linesize )
        cmats[c_idx] = np.fromfile(input_file, dtype=np.float32, count=(N*N+1))
        c_idx += 1

    return cmats[0:,0:N*N],cmats[0:,N*N]

def read_data_random(input_file,N,M,Nd,cmats):

    linesize = 4*(N*N+1)

    rint = np.random.randint(int(Nd), size=M)

    idx = 0
    for i in rint:
        input_file.seek( i * linesize )
        cmats[idx] = np.fromfile(input_file, dtype=np.float32, count=(N*N+1))

        idx += 1

    return cmats[0:,0:N*N],cmats[0:,N*N]

# Get the number of CM datapoints in file
def getNdat(file, N):
    file.seek(0,2)
    Nd = int(file.tell() / (4*(N*N+1)))
    file.seek(0)
    return Nd

def setmaxE(X,Y,E):
    m = min(X)
    newlist = []
    for i in range(0,len(X)):
        if gt.hatokcal * abs(X[i]-m) <= E:
            newlist.append(Y[i])

    return newlist

# Training params
wkdir = '/home/jujuman/Research/CMatrixBaseline_data/data05p/'
tdatafn = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/testdata/'
#tdatafn = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_10/testdata/'
paramfile = wkdir + 'parameters.pkl'
scalerfile = wkdir + 'scaler.pkl'

N = 32
#M = 100
P = 1.0

# Set this based on training sets ymin
y_min = -3.95532

# ------- Begin MLP Regressor Prediction -------
print ('Building scalers...')
# Transform data
X_scaler = joblib.load(scalerfile)

print ('Y-min: ' + str(y_min))

# Set timer
_t1b = tm.time()

# Prepare MLP model
clf = joblib.load(paramfile)

# Get files
files = listdir(tdatafn)

print ('# of test files: ' + str(len(files)))
print ('Begin prediction...')

Ecmp = []
Eact = []

#---------BEGIN VALIDATION LOOP---------
# Set square diff container
sqd = 0.0

random.shuffle(files)

batch = 0
dpts = 0
for i in files:
    batch += 1

    # Get training data
    X, y, ise = read_data_files_convert_cm(tdatafn+i,N)

    Nm = y.shape[0]

    # Transform data
    X_train = X_scaler.transform(X)

    # Compute training delta E
    y_tmp = y_min * clf.predict(X_train)# + ise

    # Undo ISE shift
    y = y# + ise
    dE = y_tmp - y

    #y_tmp = y_tmp / y_min
    #y = y / y_min

    # Restrict
    Ecmp_t = setmaxE(y, y_tmp, 300.0)
    Eact_t = setmaxE(y, y, 300.0)

    Ecmp += Ecmp_t
    Eact += Eact_t

    dpts += len(Eact_t)

    # Compute sum of squared diff
    btcherr = (dE * dE).sum()
    sqd += btcherr

    print ('File ' + str(batch) + ' of ' + str(len(files)) + ' complete. RMSE(kcal/mol): ' + "{:.5f}".format(gt.hatokcal*np.sqrt(btcherr/float(Nm))))

mse_t = sqd / float(dpts)

print('|-----------Prediction Complete-----------|')
print('  Epoch error MSE(Ha) -- ')
print('     Test: ' + "{:.5f}".format(mse_t))

print('  Epoch error RMSE(kcal/mol) -- ')
print('     Test: ' + "{:.5f}".format(gt.hatokcal*np.sqrt(mse_t)))

print('|-----------------------------------------|')

Ecmp = gt.hatokcal * np.array(Ecmp, dtype=float)
Eact = gt.hatokcal * np.array(Eact, dtype=float)

# Compute RMSE in kcal/mol
rmse = gt.calculaterootmeansqrerror(Ecmp,Eact)

# End timer
_t1e = tm.time()
print('Computation complete. Time: ' + "{:.4f}".format((_t1e - _t1b)) + 's')

# Output model information
print ('RMSE: ' + str(rmse))

# Plot
plt.scatter(Eact, Ecmp, label='Baseline RMSE: ' + '%s' % float('%.3g' % rmse) + ' kcal/mol', color='red')
plt.scatter(Eact, Eact, label='DFT', color='blue')

plt.title('GDB-10 testset - CM/MLP relative energy correlation to DFT')
plt.xlabel('E act (kcal/mol)')
plt.ylabel('E cmp (kcal/mol)')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()