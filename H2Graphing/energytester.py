__author__ = 'jujuman-home'

import numpy as np
import matplotlib.pyplot as plt
import math
import statistics as stat


# -----------------------
# readfile into np array
# -----------------------
def getfltsfromfile(filein, cols):
    # Open and read the data file
    infile = open(filein, 'r')

    infile_s = []

    for line in infile:
        row = line.strip().split(",")
        infile_s.append(row)

    # Truncate and convert to numpy array
    nparray = np.array(infile_s)
    data = nparray[:-1, cols]
    data = np.array(data, dtype=float)
    return data


# ----------------------------------------
# Rewrite Data file with outliers removed
# ----------------------------------------
def removeoutliers(filein, fileout, ex, es, en, ecol):
    # Open and read the data file
    infile = open(filein, 'r')
    otfile = open(fileout, 'w')

    vall = ex - float(en) * es
    valh = ex + float(en) * es

    cnt = int(0)
    for line in infile:
        row = line.strip().split(",")
        E = float(row[ecol])
        if vall < E < valh:
            otfile.write(line)

        cnt += 0

    otfile.close()


# ------------------------------------------
# Calculate the RMSD given two structures
# ------------------------------------------
def calculatermsd(compdata, actdata):

    n = int(compdata.shape[0])/3

    sum = float(0.0)
    for i in range(0, n):
        sum += (compdata[i*3] - actdata[i*3])**2\
            + (compdata[i*3+1] - actdata[i*3+1])**2\
            + (compdata[i*3+2] - actdata[i*3+2])**2

    rtn = math.sqrt(sum/float(n))

    return rtn


# -----------------------
#   Make index array
# -----------------------
def computermsdofdata(filein):

    infile = open(filein, 'r')

    infile_s = []

    for line in infile:
        row = line.strip().split(",")
        infile_s.append(row)

    nparray = np.array(infile_s)

    N = int(nparray[0][0])
    K = int(nparray.shape[0])
    colstart = N + 1
    colend = colstart + N * 3

    compdatastr = nparray[1, colstart:colend]
    compdataflt = np.array(compdatastr, dtype=float)

    datatmp = []
    for i in range(0, K-1):
        actdatastr = nparray[i, colstart:colend]
        actdataflt = np.array(actdatastr, dtype=float)

        datatmp.append(calculatermsd(compdataflt, actdataflt))

    data = np.array(datatmp, dtype=float)

    return data


# ------------------------------------------
# Calculate the RMSD given two structures
# ------------------------------------------
def countoutlier(data, mean, stddev, ndev):

    n = int(data.shape[0])

    vall = mean - float(ndev) * stddev
    valh = mean + float(ndev) * stddev

    cnt = int(0)
    for i in range(0, n):
        if data[i] > valh or data[i] < vall:
            cnt += 1
            print(data[i])

    return cnt

# -----------------------
#     MAIN PROGRAM
# -----------------------

data1 = getfltsfromfile('/home/jujuman-home/ServerAccess/MordorScratch/Research/ANN-Test-Data/H2O2_3/trainData.dat', [16])
data2 = computermsdofdata('/home/jujuman-home/ServerAccess/MordorScratch/Research/ANN-Test-Data/H2O2_3/trainData.dat')

X = stat.mean(data1[:, 0])
S = stat.stdev(data1[:, 0], X)

print(X)
print(S)

print(S * 2.5)
print(countoutlier(data1, X, S, 2.5))

#removeoutliers('/home/jujuman-home/ServerAccess/MordorScratch/Research/ANN-Test-Data/H2O2_3/trainData.dat'
#               , '/home/jujuman-home/ServerAccess/MordorScratch/Research/ANN-Test-Data/H2O2_3/trainDataFIXED.dat', X, S, 2.5, 25)

data4 = getfltsfromfile('/home/jujuman-home/ServerAccess/MordorScratch/Research/ANN-Test-Data/H2O2_2/trainData.dat', [16])
data5 = computermsdofdata('/home/jujuman-home/ServerAccess/MordorScratch/Research/ANN-Test-Data/H2O2_2/trainData.dat')

data3 = np.sort(data1, axis=None)
data6 = np.sort(data4, axis=None)

print(data3[0])
print(data6[0])

plt.scatter(data2, data1, label='UB3LYP/6-31g*', color='blue')
plt.scatter(data5, data4, label='UB3LYP/6-31g*', color='red')

plt.title('SCAN: H202 Energy vs. H-O-O Angle')
plt.xlabel('RMSD')
plt.ylabel('Energy (Hartrees)')
plt.legend(bbox_to_anchor=(0.5, 0.95), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()
