import numpy as np
import statsmodels.api as sm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os


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
    data = nparray[:-1, cols]
    data = np.array(data, dtype=float)
    return data


# -----------------------

# ----------------------------
# Calculate Mean Squared Diff
# ----------------------------

def calculatemeansqrdiff(data1, data2):
    data = np.power(data1 - data2, 2)
    print(np.mean(data))
    return

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

# ------------
# AM1 vs Act
# ------------
user = os.environ['USER']
dir = '/Research/ANN-Test-Data/FormaldehydeFrag/fragTrain/'

data1 = getfltsfromfile('/home/' + user + dir + 'testgraph.dat', [0])
data1 = data1 * 0.005 + 0.6
data2 = getfltsfromfile('/home/' + user + dir + 'testgraph.dat', [1])

data3 = getfltsfromfile('/home/' + user + dir + 'testgraph.dat', [0])
data3 = data3 * 0.005 + 0.6
data4 = getfltsfromfile('/home/' + user + dir + 'testgraph.dat', [2])


#data7 = getfltsfromfile('/home/'+user+'/Research/ANN-Test-Data/FormaldehydeFrag/bku2_FPNTest/afterh20.dat', [0])
#data8 = getfltsfromfile('/home/'+user+'/Research/ANN-Test-Data/FormaldehydeFrag/bku2_FPNTest/afterh20.dat', [2])


font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

# --------------
# Setup 2D Plot
# --------------
plt.scatter(data2, data4, color='blue', label='UB3LYP/6-31g* vs. ANN')
#plt.scatter(data1, data2, color='blue', label='UB3LYP/6-31g*')
#plt.scatter(data3, data4, color='red', label='NNP',linewidth=3)

#plt.scatter(data7, data8, color='orange', label='NNP After Retrain with H2O')
#plt.scatter(data1, data4, color='green', label='NNP O-H2')

#plt.scatter(data1, data4, color='red', label='M=3 UB3LYP/6-31g*')
#plt.scatter(data1, data3, color='green', label='MLNN[(2:5:2:5:2)-32]')

plt.title('UB3LYP/6-31g* vs. ANN Energy for $H_2 CO$')
#plt.title('SCAN: Formic Acid Energy vs. H-O-H Angle')
plt.xlabel('ANN Energy (Hartrees)')
plt.ylabel('UB3LYP/6-31g* Energy (Hartrees)')
#plt.legend(bbox_to_anchor=(0.1, 0.95), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()
