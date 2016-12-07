from sklearn import linear_model
from sklearn.kernel_ridge import KernelRidge
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

from os import listdir

import numpy as np
import graphtools as gt
import coulomb_matrix as cm
import matplotlib.pyplot as plt

dtdir='/home/jujuman/Research/LinearModelTesting/data/'

files = listdir(dtdir)

N = 11

X = np.empty([1,N*N],dtype=float)
y = np.empty([1],dtype=float)

for i in files:

    # Get training molecules
    xyz_tr,typ_tr,Eact_tr,readf = gt.readncdat(dtdir + i)

    # Compute energy of atoms at infinite separation
    ise = cm.computerISE(typ_tr)

    # Compute training coulomb matrices
    Xtr = np.asfortranarray( cm.GenCMatData(xyz_tr,typ_tr,N) )

    X = np.concatenate( (X, Xtr) , axis=0)
    y = np.concatenate( (y, Eact_tr - ise) , axis=0)


X = X[1:]
y = y[1:]

print(X.shape)

# For speed set as FORTRAN contiguous
X = np.asfortranarray( X )
y = np.asfortranarray( y )

# Transform data
X = StandardScaler().fit_transform(X)

y_scaler = StandardScaler().fit(y)
y = y_scaler.transform(y)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2, random_state=42)

# Prepare linear model
clf = linear_model.SGDRegressor()
#clf = KernelRidge(alpha=0.1)

# Fit model
clf.fit(X_train, y_train)

# Compute and print r^2 score
print( clf.score(X_test,y_test) )

# Store predicted energies
Ecmp = clf.predict(X_test)

Ecmp = gt.hatokcal*(y_scaler.inverse_transform(Ecmp))
Eact = gt.hatokcal*(y_scaler.inverse_transform(y_test))

#Ecmp = Ecmp - min(Ecmp)
#Eact = Eact - min(Eact)

# Compute RMSE in kcal/mol
rmse = gt.calculaterootmeansqrerror(Ecmp,Eact)

# Output model information

print ('RMSE: ' + str(rmse))
#print(clf.coef_)
#print(clf.intercept_)

# Plot
plt.scatter(Eact, Ecmp, label='CM/SGDR RMSE: ' + '%s' % float('%.3g' % rmse) + ' kcal/mol', color='red')
plt.scatter(Eact, Eact, label='DFT', color='blue')

plt.title('GDB-2 - CM/SGDR energy correlation to DFT')
plt.xlabel('E act (kcal/mol)')
plt.ylabel('E cmp (kcal/mol)')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()