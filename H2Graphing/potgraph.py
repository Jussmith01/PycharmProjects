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
dir = '/Research/ANN-Test-Data/GDB-11/train2/'

file = 'graph.dat'

data1 = getfltsfromfile('/home/' + user + dir + file, [0])
data1 = data1 * 2.0 + -180.0
data2 = getfltsfromfile('/home/' + user + dir + file, [1])
data3 = getfltsfromfile('/home/' + user + dir + file, [2])

dir = '/Research/ANN-Test-Data/GDB-11/train3/'
data4 = getfltsfromfile('/home/' + user + dir + file, [2])
#AM1 = getfltsfromfile('/home/' + user + dir + file, [1])

dir = '/Research/ANN-Test-Data/GDB-11/train4/'
data5 = getfltsfromfile('/home/' + user + dir + file, [2])
#PM6 = getfltsfromfile('/home/' + user + dir + file, [1])

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

#print(AM1.shape[0])
#AM1 = AM1 - AM1[AM1.shape[0]-1]
#PM6 = PM6 - PM6[AM1.shape[0]-1]
#data2 = data2 - data2[AM1.shape[0]-1]
#data3 = data3 - data3[AM1.shape[0]-1]
#data4 = data4 - data4[AM1.shape[0]-1]
#data5 = data5 - data5[AM1.shape[0]-1]

# --------------
# Setup 2D Plot
# --------------
#plt.plot(data1, AM1, color='black', label='AM1',linewidth=4)
#plt.plot(data1, PM6, color='grey', label='PM6',linewidth=4)
plt.plot(data2, data2, color='blue', label='B3LYP/6-31g*',linewidth=2)
plt.scatter(data2, data3, color='red', label='ANN - up to GDB-2',linewidth=4)
plt.scatter(data2, data4, color='orange', label='ANN - up to GDB-3',linewidth=4)
plt.scatter(data2, data5, color='green', label='ANN - up to GDB-4',linewidth=4)
#plt.plot(data1, data3, color='red', label='ANN - up to GDB-3',linewidth=4)
#plt.plot(data1, data4, color='green', label='ANN - up to GDB-4',linewidth=4)

#plt.plot(data1, data4, color='green', label='ANN - only $H_3 C H_2 C H_2 N$',linewidth=4)
#plt.plot(data1, data4, color='orange', label='ANN - up to GDB-3',linewidth=4)
#plt.scatter(data1, data6, color='black', label='ANN GDB-3',linewidth=4)

#plt.title('Random points on the glycine potential surface')
#plt.title('SCAN: Formic Acid Energy vs. H-O-H Angle')
#plt.xlabel('E Target (Hartrees)')
#plt.ylabel('E Computed (Hartrees)')
#plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.)


# -----
# PLOT
# -----
plt.show()
