__author__ = 'jujuman'

import numpy as np
#import statsmodels.api as sm

# -----------------------
# readfile into np array
# -----------------------
def getfltsfromfile(file, delim, cols):
    # Open and read the data file
    infile = open(file, 'r')

    infile_s = []

    for line in infile:
        row = line.strip().split(delim)
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

def calculateelementdiff(data1):
    C = np.shape(data1)[0]
    data = np.zeros((C-1, 2))
    for i in range(0, C-1):
        data[i,0] = i
        data[i,1] = (data1[i] - data1[i+1])

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
