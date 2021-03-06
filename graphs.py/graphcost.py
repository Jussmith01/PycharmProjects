import numpy as np
import statsmodels.api as sm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


# -----------------------
# readfile into np array
# -----------------------
def getfltsfromfile(file, cols):
    # Open and read the data file
    infile = open(file, 'r')

    infile_s = []

    for line in infile:
        row = line.strip().split(",")
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
    for i in range(0, C-1):
        data[i,0] = i
        data[i,1] = (data1[i] - data1[i+1]) / dx

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
data1 = getfltsfromfile('/home/jujuman-home/Gits/ForcePredictionNetwork/g09DNNTSData/H2O2Energy/GPU-1/costgraph.out', [0])
data2 = getfltsfromfile('/home/jujuman-home/Gits/ForcePredictionNetwork/g09DNNTSData/H2O2Energy/GPU-1/costgraph.out', [1])
data3 = getfltsfromfile('/home/jujuman-home/Gits/ForcePredictionNetwork/g09DNNTSData/H2O2Energy/GPU-1/costgraph.out', [2])
#data4 = getfltsfromfile('/home/jujuman-home/Gits/ForcePredictionNetwork/g09DNNTSData/H2OEnergy/GPU-1/costgraphsmax.out', [1])
#data5 = getfltsfromfile('/home/jujuman-home/Gits/ForcePredictionNetwork/g09DNNTSData/H2OEnergy/GPU-1/costgraphsmax.out', [2])

# --------------
# Setup 2D Plot
# --------------
plt.scatter(data1, data2, color='green', label='Training Cost - GFL')
plt.scatter(data1, data3, color='red', label='Validation Cost - GFL')
#plt.scatter(data1, data4, color='blue', label='Training Cost - Smax')
#plt.scatter(data1, data5, color='orange', label='Validation Cost - Smax')

plt.title('Cost Comparison')
plt.xlabel('Epoch')
plt.ylabel('Cost')
plt.legend(bbox_to_anchor=(0.7, 0.95), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()
