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
data1 = getfltsfromfile('/home/jujuman/ServerAccess/homegtx770s/Gits/ForcePredictionNetwork/g09DNNTSData/H2O2Energy/GPU-1-Int-e/scanbondOO.dat', [4])
data2 = getfltsfromfile('/home/jujuman/ServerAccess/homegtx770s/Gits/ForcePredictionNetwork/g09DNNTSData/H2O2Energy/GPU-1-Int-e/scanbondOO.dat', [9])
data3 = getfltsfromfile('/home/jujuman/ServerAccess/homegtx770s/Gits/ForcePredictionNetwork/g09DNNTSData/H2O2Energy/GPU-1-Int-e/graph_act.dat', [0])

#data3 = getfltsfromfile('/home/jujuman/ServerAccess/homegtx770s/Gits/ForcePredictionNetwork/g09DNNTSData/H2O2Energy/GPU-1-Int-e/graph_act.dat', [0])
#data4 = getfltsfromfile('/home/jujuman/Gits/ForcePredictionNetwork/g09DNNTSData/H2O2Energy/GPU-1-b/graph_actOH2.dat', [0])

#data4 = getfltsfromfile('/home/jujuman/Gits/g09DNNTSBuilder/bin/Release/c5uhf.dat', [15])

# --------------
# Setup 2D Plot
# --------------
plt.scatter(data1, data2, label='UB3LYP/6-31g*')
plt.scatter(data1, data3, color='red', label='NNP O-H1')
#plt.scatter(data1, data4, color='green', label='NNP O-H2')

#plt.scatter(data1, data4, color='red', label='M=3 UB3LYP/6-31g*')
#plt.scatter(data1, data3, color='green', label='MLNN[(2:5:2:5:2)-32]')

plt.title('SCAN: H202 Energy vs. O-O Bond')
plt.xlabel('O-O Bond (Angstroms)')
plt.ylabel('Energy (Hartrees)')
plt.legend(bbox_to_anchor=(0.5, 0.95), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()

# --------------
# Setup 3D Plot
# --------------
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(data1,data2,data3,depthshade=True)
#ax.scatter(data1,data2,data4,color='red',depthshade=True)

#ax.set_xlabel('Bond O-H1 (Angstroms)')
#ax.set_ylabel('Bond O-H2 (Angstroms)')
#ax.set_zlabel('HF/6-31g* Force on H1 (Hartree/Bohr)')

# ------------
# Numer Deriv
# ------------
derivHF = calculatenumderiv(data2,0.5)
plt.scatter(derivHF[:, 0], derivHF[:, 1], label='UB3LYP/6-31g*')

derivMLNN1 = calculatenumderiv(data3,0.5)
plt.scatter(derivMLNN1[:, 0], derivMLNN1[:, 1],color='red', label='NNP')

#derivMLNN2 = calculatenumderiv(data4,0.5)
#plt.scatter(derivMLNN2[:, 0]/2.0, derivMLNN2[:, 1],color='red', label='Deriv MLNN[(2:5:2:5:2)-32]')

#derivMLNN = calculatenumderiv(data5,0.002)
#plt.scatter(derivMLNN[:, 0], derivMLNN[:, 1],color='red', label='Deriv MLNN')

plt.title('SCAN: H202 Force vs. O-O Bond')
plt.xlabel('O-O Bond Scan Step')
plt.ylabel('Force')
plt.legend(bbox_to_anchor=(0.5, 0.95), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()

# ----------------
# Numer Second Deriv
# ----------------
#derivHF2 = calculatenumderiv(derivHF[:, 1],1.0)
#plt.scatter(derivHF2[:, 0], derivHF2[:, 1], label='B3LYP/6-31g*')

#derivMLNN12 = calculatenumderiv(derivMLNN1[:, 1],1.0)
#plt.scatter(derivMLNN12[:, 0], derivMLNN12[:, 1],color='green', label='MLNN[(2:4:2:4:4:2)-32]')

#derivMLNN22 = calculatenumderiv(derivMLNN2[:, 1],1.0)
#plt.scatter(derivMLNN22[:, 0], derivMLNN22[:, 1],color='red', label='Deriv HF')

#derivMLNN = calculatenumderiv(data5,0.002)
#plt.scatter(derivMLNN[:, 0], derivMLNN[:, 1],color='red', label='Deriv MLNN')

#plt.title('SCAN: H20 O-H1 Energy Second Numer. Derivative per Scan Step')
#plt.xlabel('Scan Step')
#plt.ylabel('Energy Second Derivative')
#plt.legend(bbox_to_anchor=(0.7, 0.95), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
#plt.show()