from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import normalize
from sklearn.model_selection import train_test_split

from sklearn.neural_network import MLPRegressor

from sklearn.externals import joblib

from os import listdir
import time as tm

import numpy as np
import graphtools as gt
import matplotlib.pyplot as plt

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

# Training params
wkdir = '/home/jujuman/Research/CMatrixBaseline_data/data50p/'
tdatafn = wkdir + 'test_cm_data.dat'
paramfile = wkdir + 'parameters.pkl'
scalerfile = wkdir + 'scaler.pkl'

N = 32
M = 10000
P = 1.0

# Set this based on training sets ymin
y_min = -3.96939

# ------- Begin MLP Regressor Prediction -------
t_file = open(tdatafn, 'r')

Ntd = getNdat(t_file,N)
print ('Test data: ' + str(Ntd))

print ('Building scalers...')
# Transform data
X_scaler = joblib.load(scalerfile)

print ('Y-min: ' + str(y_min))

# Set timer
_t1b = tm.time()

# Initilize working memory
t_cmats = np.empty([M,N*N+1],dtype=np.float32)

# Prepare MLP model
clf = joblib.load(paramfile)

# Create data index
t_Ainx = np.arange(Ntd)

# Rand Init
np.random.shuffle(t_Ainx)

# Get number of batches
Ntb = int( P * np.floor(Ntd/M) )

print ('# of test batches: ' + str(Ntb))
print ('Begin prediction...')

#---------BEGIN VALIDATION LOOP---------
# Set square diff container
sqd = 0.0

batch = 0
for i in range(0,Ntb):
    batch += 1

    # Get training data
    X, y = read_data_indexed(t_file,N,M,Ntd,i,t_Ainx,t_cmats)

    # Transform data
    X_train = X_scaler.transform(X)

    # Compute training delta E
    dE = y_min * clf.predict(X_train) - y

    # Compute sum of squared diff
    btcherr = (dE * dE).sum()
    sqd += btcherr

    print ('Batch ' + str(batch) + ' of ' + str(Ntb) + ' complete. RMSE(kcal/mol): ' + "{:.5f}".format(gt.hatokcal*np.sqrt(btcherr/float(M))))

mse_t = sqd / float(Ntb*M)

print('|-----------Prediction Complete-----------|')
print('  Epoch error MSE(Ha) -- ')
print('     Test: ' + "{:.5f}".format(mse_t))

print('  Epoch error RMSE(kcal/mol) -- ')
print('     Test: ' + "{:.5f}".format(gt.hatokcal*np.sqrt(mse_t)))

print('|-----------------------------------------|')

# Compute and print r^2 score
#print( clf.score(X_test,y_test) )

# Store predicted energies
#Ecmp = clf.predict(X_test)

#Ecmp = gt.hatokcal*(y_min*Ecmp)
#Eact = gt.hatokcal*(y_test)

# Compute RMSE in kcal/mol
#rmse = gt.calculaterootmeansqrerror(Ecmp,Eact)

# End timer
_t1e = tm.time()
print('Computation complete. Time: ' + "{:.4f}".format((_t1e - _t1b)) + 's')

# Output model information
#print ('RMSE: ' + str(rmse))
#print(clf.coef_)
#print(clf.intercept_)

# Plot
#plt.scatter(Eact, Ecmp, label='Baseline RMSE: ' + '%s' % float('%.3g' % rmse) + ' kcal/mol', color='red')
#plt.scatter(Eact, Eact, label='DFT', color='blue')

#plt.title('GDB-2 - CM/MLP energy correlation to DFT')
#plt.xlabel('E act (kcal/mol)')
#plt.ylabel('E cmp (kcal/mol)')
#plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
#plt.show()