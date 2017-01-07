from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import normalize
from sklearn.model_selection import train_test_split

from sklearn.neural_network import MLPRegressor

from sklearn.externals import joblib

from os import listdir
import time as tm
import sys

import numpy as np
import graphtools as gt
import matplotlib.pyplot as plt

def progress(count, total, suffix=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', suffix))
    sys.stdout.flush()  # As suggested by Rom Ruben

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

#-------------------------------------------------
# Training params
wkdir = '/home/jujuman/Research/CMatrixBaseline_data/data100p/'
tdatafn = wkdir + 'train_cm_data.dat'
vdatafn = wkdir + 'valid_cm_data.dat'
paramfile = wkdir + 'parameters.pkl'
scalerfile = wkdir + 'scaler.pkl'

N = 32
M = 20 * 1024
P = 0.5
tol_set = 3
load_model = False
verbose = False
#------------------------------------------------


# ------- Begin MLP Regressor -------
t_file = open(tdatafn, 'r')
v_file = open(vdatafn, 'r')

Ntd = getNdat(t_file,N)
Nvd = getNdat(v_file,N)
print ('Training   data: ' + str(Ntd))
print ('Validation data: ' + str(Nvd))

Mt = int( 0.05 * Ntd )
t_cmats = np.empty([Mt,N*N+1],dtype=np.float32)

print ('Scaling with: ' + str(Mt))

# Read scaling data
X, y = read_data_random(t_file,N,Mt,Ntd,t_cmats)

print ('Building scalers...')
# Transform data
X_scaler = StandardScaler().fit(X)

# Dump scaler for testing
joblib.dump(X_scaler, scalerfile)

y_max = y.max()
y_min = y.min()
print ('Max/min: ' + str(y_max) + '/' + str(y_min))

# Set timer
_t1b = tm.time()

# Initilize working memory
t_cmats = np.empty([M,N*N+1],dtype=np.float32)
X = np.empty([M,N*N],dtype=np.float32)
y = np.empty([M],dtype=np.float32)

# Prepare MLP model
if load_model:
    clf = joblib.load(paramfile)
else:
    clf = MLPRegressor(hidden_layer_sizes=(128,128,64,),activation='tanh',solver='adam',verbose=False,batch_size=1024,max_iter=1,tol=1.0E-7)

# Create data index
t_Ainx = np.arange(Ntd)
v_Ainx = np.arange(Nvd)

# Rand Init
np.random.shuffle(v_Ainx)

# Get number of batches
Ntb = int( P * np.floor(Ntd/M) )
Nvb = int( P * np.floor(Nvd/M) )

#Ntb = 10
#Nvb = 10

print ('# of train batches: ' + str(Ntb))
print ('# of valid batches: ' + str(Nvb))

Epoch = 0
train = True

tol = tol_set

best_mse_v = 100.0

print ('Begin regression...')
while train == True:

    # Set Epoch timer
    _e1b = tm.time()

    Epoch += 1

    #---------BEGIN TRAINING LOOP---------
    # Shuffle data index
    np.random.shuffle(t_Ainx)

    # Set square diff container
    sqd = 0.0

    bcnt = 0
    for i in range(0,Ntb):
        bcnt += 1

        # Get training data
        _tet1b = tm.time()
        X, y = read_data_indexed(t_file,N,M,Ntd,i,t_Ainx,t_cmats)
        _tet1e = tm.time()

        # Transform data
        _tet2b = tm.time()
        X = X_scaler.transform(X)
        y = y/y_min
        _tet2e = tm.time()

        # Fit model
        _tet3b = tm.time()
        clf.partial_fit(X, y)
        _tet3e = tm.time()

        # Compute training delta E
        _tet4b = tm.time()
        dE = y_min * clf.predict(X) - y_min * y
        _tet4e = tm.time()

        # Compute sum of squared diff
        btcherr = (dE * dE).sum()
        sqd += btcherr

        progress(bcnt,Ntb,'Training progress')

        if verbose:
            print('\n  Batch: ' + str(bcnt) + ' of ' + str(Ntb))
            print('  Read Time: ' + "{:.4f}".format((_tet1b - _tet1e)) + 's')
            print('  Scal Time: ' + "{:.4f}".format((_tet2b - _tet2e)) + 's')
            print('  Pfit Time: ' + "{:.4f}".format((_tet3b - _tet3e)) + 's')
            print('  Pred Time: ' + "{:.4f}".format((_tet4b - _tet4e)) + 's')
            print('  loss: ' + str(y_min * clf.loss_) + ' mse: ' + str(btcherr/float(M)))

    mse_t = sqd / (Ntb*M)

    print('\n')

    #---------BEGIN VALIDATION LOOP---------
    # Set square diff container
    sqd = 0.0
    bcnt = 0
    for i in range(0,Nvb):
        bcnt += 1

        # Get training data
        X, y = read_data_indexed(v_file,N,M,Nvd,i,v_Ainx,t_cmats)

        # Transform data
        X = X_scaler.transform(X)

        # Compute training delta E
        dE = y_min * clf.predict(X) - y

        # Compute sum of squared diff
        sqd += (dE * dE).sum()

        progress(bcnt,Nvb,'Validation progress')

    mse_v = sqd / (Nvb*M)

    # End timer
    _e1e = tm.time()

    tol -= 1

    if mse_v < best_mse_v:
        best_mse_v = mse_v

        tol = tol_set

        # Offload model to disk
        joblib.dump(clf, paramfile)

    print('\n|-----------Epoch ' + str(Epoch) + '-----------|')
    print('  Tolerance: ' + str(tol) + ' of ' + str(tol_set))
    print('\n  Epoch error MSE(Ha) -- ')
    print('     Train: ' + "{:.7f}".format(mse_t))
    print('     Valid: ' + "{:.7f}".format(mse_v))

    print('\n  Epoch error RMSE(kcal/mol) -- ')
    print('     Train: ' + "{:.5f}".format(gt.hatokcal*np.sqrt(mse_t)))
    print('     Valid: ' + "{:.5f}".format(gt.hatokcal*np.sqrt(mse_v)))

    print('\n  Epoch Time: ' + "{:.4f}".format((_e1e - _e1b)) + 's')

    print('|------------------------------|\n')

    if tol == 0:
        train = False


print ('y-min: ' + str(y_min))

# End timer
_t1e = tm.time()
print('Computation complete. Time: ' + "{:.4f}".format((_t1e - _t1b)) + 's')
