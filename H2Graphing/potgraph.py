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
    #print (nparray)
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
dir = '/Research/ANN-Test-Data/GDB-11-B3LYP-6-31gd/train4/'

file = 'fixmolecule-2_valid.dat_graph.dat'

data1 = getfltsfromfile('/home/' + user + dir + file, [0])
data2 = getfltsfromfile('/home/' + user + dir + file, [1])
data3 = getfltsfromfile('/home/' + user + dir + file, [2])


font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)


# --------------
# Setup 2D Plot
# --------------
plt.plot(data1, data2, color='blue', label='AM1',linewidth=1)
plt.scatter(data1, data2, color='blue', label='AM1',linewidth=4)
plt.plot(data1, data3, color='orange', label='AM1',linewidth=1)
plt.scatter(data1, data3, color='orange', label='AM1',linewidth=4)

plt.title('Formic Acid Reaction Scan H-O1 -> H-O2')
#plt.title('SCAN: Formic Acid Energy vs. H-O-H Angle')
plt.xlabel('Reaction Coordinate (Angstroms)')
plt.ylabel('Energy (Hartrees)')
plt.legend(bbox_to_anchor=(0.7, 0.3), loc=2, borderaxespad=0.)


# -----
# PLOT
# -----
plt.show()
