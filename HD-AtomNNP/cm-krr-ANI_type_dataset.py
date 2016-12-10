from sklearn import linear_model
from sklearn.kernel_ridge import KernelRidge
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import normalize
from sklearn.model_selection import train_test_split

from sklearn.neural_network import MLPRegressor

from sklearn import svm

from os import listdir
import time as tm

import numpy as np
import graphtools as gt
import coulomb_matrix as cm
import matplotlib.pyplot as plt

dtdir='/home/jujuman/Dropbox/Research/LinearModelTesting/data/'
#dtdir='/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_02/testdata/'
files = listdir(dtdir)

N = 17

X = np.empty([1,N*N],dtype=float)
y = np.empty([1],dtype=float)

print('Computing...')
_t1b = tm.time()

cnt=0

for i in files:
    cnt += 1
    print('FILE: ' + str(cnt) + ' of ' + str(len(files)))

    # Set file
    file = dtdir + i

    # Get training molecules
    xyz_tr, typ_tr, Eact_tr, readf = gt.readncdat(file, np.float32)

    # Compute energy of atoms at infinite separation
    ise = cm.computerISE(typ_tr)

    Nm = Eact_tr.shape[0]
    xyz_tr= xyz_tr[0:int(1.0*Nm)]
    Eact_tr = Eact_tr[0:int(1.0*Nm)]

    # Compute training coulomb matrices
    Xtr = np.asfortranarray(cm.GenCMatData(xyz_tr, typ_tr, N))

    Eact_tr = (Eact_tr - ise)

    # Add data to pot
    X = np.concatenate((X, Xtr), axis=0)
    y = np.concatenate((y, Eact_tr), axis=0)

X = X[1:]
y = y[1:]



print(X.shape)

# For speed set as FORTRAN contiguous
#X = np.asfortranarray( X )
#y = np.asfortranarray( y )

# Transform data
X = StandardScaler().fit_transform(X)

#y_scaler = StandardScaler().fit(y)
#y = y_scaler.transform(y)

#y = normalize(y,axis=0)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.1, random_state=42)

# Prepare linear model
#clf = linear_model.SGDRegressor()
#clf = svm.SVR(kernel='linear',verbose=True)ru
clf = MLPRegressor(hidden_layer_sizes=(16,16,16,),activation='tanh',solver='adam',verbose=True,tol=1.0E-10,early_stopping=True)

'''
font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 9}

plt.rc('font', **font)

fig, axes = plt.subplots(nrows=1, ncols=1)
axes.set_title("Data: " + file)
axes.set_ylabel('Normalized distant count')
axes.set_xlabel('Distance ($\AA$)')

axes.hist(y_train, 150, color='blue',normed=True, label='plot',linewidth=2,alpha=1.0)
plt.show()
'''

# Fit model
clf.fit(X_train, y_train)

# Compute and print r^2 score
print( clf.score(X_test,y_test) )

# Store predicted energies
Ecmp = clf.predict(X_test)

#Ecmp = gt.hatokcal*(y_scaler.inverse_transform(Ecmp))
#Eact = gt.hatokcal*(y_scaler.inverse_transform(y_test))
Ecmp = gt.hatokcal*(Ecmp)
Eact = gt.hatokcal*(y_test)

# Compute RMSE in kcal/mol
rmse = gt.calculaterootmeansqrerror(Ecmp,Eact)

# End timer
_t1e = tm.time()
print('Computation complete. Time: ' + "{:.4f}".format((_t1e - _t1b)) + 's')

# Output model information
print ('RMSE: ' + str(rmse))
#print(clf.coef_)
#print(clf.intercept_)

# Plot
plt.scatter(Eact, Ecmp, label='CM/SGDR RMSE: ' + '%s' % float('%.3g' % rmse) + ' kcal/mol', color='red')
plt.scatter(Eact, Eact, label='DFT', color='blue')

plt.title('GDB-2 - CM/MLP energy correlation to DFT')
plt.xlabel('E act (kcal/mol)')
plt.ylabel('E cmp (kcal/mol)')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()