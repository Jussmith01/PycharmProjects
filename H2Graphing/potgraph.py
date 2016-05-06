import numpy as np
import statsmodels.api as sm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
import random

# -----------------------
# readfile into np array
# -----------------------
def getfltsfromfile(file, cols):
    # Open and read the data file
    infile = open(file, 'r')

    infile_s = []

    for line in infile:
        row = line.strip().split(" ")
        infile_s.append(row)

    # Truncate and convert to numpy array
    nparray = np.array(infile_s)
    data = nparray[:, cols]
    data = np.array(data, dtype='f8')
    return data

# -----------------------

# ----------------------------
# Calculate Mean Squared Diff
# ----------------------------

def calculatemeansqrerror(data1, data2):
    data = np.power(data1 - data2, 2)
    return np.mean(data)

def calculaterootmeansqrerror(data1, data2):
    data = np.power(data1 - data2, 2)
    return np.sqrt(np.mean(data))

# ----------------------------
# Calculate Mean Squared Diff
# ----------------------------

def calculatenumderiv(data1, dx):
    C = np.shape(data1)[0]
    data = np.zeros((C, 2))
    for i in range(1, C-1):
        data[i,0] = i
        data[i,1] = (data1[i-1] - data1[i+1]) / (2.0*dx)

    return data

def calculatemean(data1):
    C = np.shape(data1)[0]
    print (C)
    sum = 0.0
    for i in data1:
        sum += i

    return sum/float(C)

def calculateabsdiff(data1):
    C = int(np.shape(data1)[0]/2)
    data = np.zeros((C, 2))
    for i in range(1, C-1):
        data[i,0] = i
        data[i,1] = data1[i*2+1] - data1[i*2]

    return data

# -----------------------

# ----------------------------
# Linear-ize Data
# ----------------------------
def makedatalinear(datain):
    data = np.zeros((np.shape(datain)[0]*np.shape(datain)[1],1))
    #data = datain[:,0]

    C = np.shape(datain)[0]
    R = np.shape(datain)[1]
    print (C, "X", R)

    i = j = 0
    for i in range(0, C):
        for j in range(0, R):
            data[i*3+j] = datain[i, j]

    return data

# -----------------------
cmap = mpl.cm.brg
# ------------
# AM1 vs Act
# ------------
user = os.environ['USER']
#dir = '/Research/ANN-Test-Data/GDB-11/train5/'
#dir = '/Research/ANN-Test-Data/GDB-11/train3/'
'''
N = 8
for i in range(0,N):
    file = 'gdb11_s02-' + str(i) + '_train.dat_graph.dat'

    data1 = getfltsfromfile('/home/' + user + dir + file, [0])
    data2 = getfltsfromfile('/home/' + user + dir + file, [1])rm
    data3 = getfltsfromfile('/home/' + user + dir + file, [2])

    data2 = (data3 - data2)*(data3 - data2)

    print('Datasize: ' + str(data2.shape[0]))

    font = {'family' : 'Bitstream Vera Sans',
            'weight' : 'normal',
            'size'   : 14}

    plt.rc('font', **font)

    s = 288

    plt.plot(data1, data2, color=cmap((i+1)/float(N)), label=str(i),linewidth=1)
    plt.scatter(data1, data2, color=cmap((i+1)/float(N)), label=str(i),linewidth=2)
    #plt.scatter(data2, data3, color=cmap((i+1)/float(N)), label=str(i),linewidth=1)
'''

dir = '/Research/ANN-Test-Data/GDB-11/train3/'
#file = 'L-glutamic-acid.dat_graph.dat'
file = 'peptide-1.dat_graph.dat'

data1 = getfltsfromfile('/home/' + user + dir + file, [0])
data2 = getfltsfromfile('/home/' + user + dir + file, [1])
data3 = getfltsfromfile('/home/' + user + dir + file, [2])

dir = '/Research/ANN-Test-Data/GDB-11/train4/'
data4 = getfltsfromfile('/home/' + user + dir + file, [2])

dir = '/Research/ANN-Test-Data/GDB-11/train5/'
data5 = getfltsfromfile('/home/' + user + dir + file, [2])

dir = '/Research/ANN-Test-Data/GDB-11/train6/'
data6 = getfltsfromfile('/home/' + user + dir + file, [2])

mean2 = np.mean(data2)
mean3 = np.mean(data3)
mean4 = np.mean(data4)
mean5 = np.mean(data5)
mean6 = np.mean(data6)

data3 = data3 - (mean6 - mean2)
data4 = data4 - (mean6 - mean2)
data5 = data5 - (mean6 - mean2)
data6 = data6 - (mean6 - mean2)

print ((mean6 - mean2))
print('RMSE: ',calculaterootmeansqrerror(data3,data2),' MSE: ', calculatemeansqrerror(data3,data2))
print('RMSE: ',calculaterootmeansqrerror(data4,data2),' MSE: ', calculatemeansqrerror(data4,data2))
print('RMSE: ',calculaterootmeansqrerror(data5,data2),' MSE: ', calculatemeansqrerror(data5,data2))
print('RMSE: ',calculaterootmeansqrerror(data6,data2),' MSE: ', calculatemeansqrerror(data6,data2))

#data2 = np.log10(data2)
#data3 = np.log10(data3)

#data2 = (data3 - data2)*(data3 - data2)
#data3 = np.log10(data3)

print('Datasize: ' + str(data1.shape[0]))

font = {'family' : 'Bitstream Vera Sans',
            'weight' : 'normal',
            'size'   : 14}

plt.rc('font', **font)

plt.plot(data2, data2, color='blue', label='B3LYP/6-31G*',linewidth=2)
plt.scatter(data2, data3, color='red', label='ANN - GDB3',linewidth=4)
plt.scatter(data2, data4, color='orange', label='ANN - GDB4',linewidth=4)
plt.scatter(data2, data5, color='purple', label='ANN - GDB5',linewidth=4)
plt.scatter(data2, data6, color='green', label='ANN - GDB6',linewidth=4)

#plt.title("H-Gly-Pro-Hyp-Gly-Ala-Gly-OH")
plt.title("H-Gly-Pro-Hyp-Gly-Ala-Gly-OH")
plt.xlabel('Target E B3LYP/6-31Gd (Hartree)')
plt.ylabel('Actual E (Hartree)')
plt.legend(bbox_to_anchor=(0.025, 0.975), loc=2, borderaxespad=0.)


# -----
# PLOT
# -----
plt.show()
