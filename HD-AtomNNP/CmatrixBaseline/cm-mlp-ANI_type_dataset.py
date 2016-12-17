from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

from sklearn.neural_network import MLPRegressor

from os import listdir
import time as tm

import numpy as np
import graphtools as gt
import matplotlib.pyplot as plt

def read_data(datadir,N,P=1.0):
    files = listdir(datadir)

    cmats = np.empty([1,N*N+1],dtype=np.float32)

    cnt = 0
    for i in files:
        cnt += 1
        print ('File ' + str(cnt) + ' of ' + str(len(files)))
        input_file = open(datadir + i, 'r')

        float_array = np.fromfile(input_file, dtype=np.float32)
        float_array = float_array.reshape(int(float_array.shape[0]/(N*N+1)),N*N+1)

        float_array = float_array[0:int(P*float_array.shape[0])]

        cmats = np.concatenate((cmats,float_array), axis=0)

    return cmats[1:,0:N*N],cmats[1:,N*N]

cmdfiledir = '/home/jujuman/Research/CMatrixBaseline_data/'

N = 30

X, y = read_data(cmdfiledir,N,0.05)

print(X.shape)

_t1b = tm.time()

# Transform data
X = StandardScaler().fit_transform(X)

#y_scaler = StandardScaler().fit(y)
#y = y_scaler.transform(y)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.1, random_state=42)

# Prepare linear model
clf = MLPRegressor(hidden_layer_sizes=(128,128,64,),activation='tanh',solver='adam',verbose=True,tol=1.0E-10,early_stopping=True)

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
plt.scatter(Eact, Ecmp, label='Baseline RMSE: ' + '%s' % float('%.3g' % rmse) + ' kcal/mol', color='red')
plt.scatter(Eact, Eact, label='DFT', color='blue')

plt.title('GDB-2 - CM/MLP energy correlation to DFT')
plt.xlabel('E act (kcal/mol)')
plt.ylabel('E cmp (kcal/mol)')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()