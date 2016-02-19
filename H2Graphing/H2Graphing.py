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
dir = '/Research/ANN-Test-Data/FormaldehydeFrag/graphs/'

data1 = getfltsfromfile('/home/' + user + dir + 'nnphcoangleretrainfrag.dat', [0])
data1 = data1 + 27
data2 = getfltsfromfile('/home/' + user + dir + 'nnphcoangleretrainfrag.dat', [1])

data3 = getfltsfromfile('/home/' + user + dir + 'nnphcoangleretrainfrag.dat', [0])
data3 = data3 + 27
data4 = getfltsfromfile('/home/' + user + dir + 'nnphcoangleretrainfrag.dat', [2])

data5 = getfltsfromfile('/home/' + user + dir + 'nnphcoanglenofixfrag.dat', [0])
data5 = data5 * 2 + 27
data6 = getfltsfromfile('/home/' + user + dir + 'nnphcoanglenofixfrag.dat', [2])

data7 = getfltsfromfile('/home/' + user + dir + 'nnphcoangleretrainnofrag.dat', [0])
data7 = data7 + 27
data8 = getfltsfromfile('/home/' + user + dir + 'nnphcoangleretrainnofrag.dat', [2])


font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

# --------------
# Setup 2D Plot
# --------------
plt.scatter(data1, data2, color='blue', label='UB3LYP/6-31g*')
plt.plot(data3, data4, color='red', label='ANN H-C-O Angle Retrain - Frag',linewidth=3)
plt.plot(data5, data6, color='green', label='ANN Before Data set Fix - Frag',linewidth=2)
plt.plot(data7, data8, color='orange', label='ANN H-C-O Angle Retrain - No Frag',linewidth=2)

plt.title(r'SCAN: $\mathrm{H_2CO}$ Energy vs. H-C-O Angle')
plt.xlabel('H-C-O Angle (Degrees)')
plt.ylabel('Energy (Hartrees)')
plt.legend(bbox_to_anchor=(0.38, 0.95), loc=2, borderaxespad=0.)

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
#derivHF = calculatenumderiv(data2,0.5)
#plt.scatter(derivHF[:, 0], derivHF[:, 1], label='UB3LYP/6-31g*')

#derivMLNN1 = calculatenumderiv(data3,0.5)
#plt.scatter(derivMLNN1[:, 0], derivMLNN1[:, 1],color='red', label='NNP')

#derivMLNN2 = calculatenumderiv(data4,0.5)
#plt.scatter(derivMLNN2[:, 0]/2.0, derivMLNN2[:, 1],color='red', label='Deriv MLNN[(2:5:2:5:2)-32]')

#derivMLNN = calculatenumderiv(data5,0.002)
#plt.scatter(derivMLNN[:, 0], derivMLNN[:, 1],color='red', label='Deriv MLNN')

#plt.title('SCAN: H202 Force vs. H-O-O Angle')
#plt.xlabel('H-O-O Angle Scan Step')
#plt.ylabel('Force')
#plt.legend(bbox_to_anchor=(0.5, 0.95), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
#plt.show()

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
