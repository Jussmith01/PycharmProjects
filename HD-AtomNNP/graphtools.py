__author__ = 'jujuman'

import numpy as np
#import statsmodels.api as sm
import re

hatokcal = 627.509469

convert = hatokcal  # Ha to Kcal/mol

def readxyz (file):
    xyz = []
    typ = []
    Na = []

    fd = open(file, 'r').read()

    rb = re.compile('(\d+?)\s*?\n((?:[A-Z][a-z]?\s+?\S+?\s+?\S+?\s+?\S+?\s)+)')
    ra = re.compile('([A-Z][a-z]?)\s+?(\S+?)\s+?(\S+?)\s+?(\S+)')

    s = rb.findall(fd)

    for i in s:
        Na.append(int(i[0]))
        atm = ra.findall(i[1])

        ntyp = []
        nxyz = []
        for j in range(0, int(i[0])):
            ntyp.append(atm[j][0])
            nxyz.append(float(atm[j][1]))
            nxyz.append(float(atm[j][2]))
            nxyz.append(float(atm[j][3]))

        xyz.append(nxyz)
        typ.append(ntyp)

    return xyz,typ,Na

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

def getfltsfromfileprob(file, delim, col1, col2, prob):
    # Open and read the data file
    infile = open(file, 'r')

    infile_s = []

    for line in infile:
        if np.random.binomial(1,prob):
            row = line.strip().split(delim)
            infile_s.append(row)

    # Truncate and convert to numpy array
    nparray = np.array(infile_s)
    data1 = nparray[:, col1]
    data2 = nparray[:, col2]

    data1 = np.array(data1, dtype='f8')
    data2 = np.array(data2, dtype='f8')
    return data1,data2

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

def calculateelementdiff2D(data1):
    C = np.shape(data1)[0]
    x = np.zeros((C*C))
    y = np.zeros((C*C))
    z = np.zeros((C*C))
    d = np.zeros((C*(C-1)/2))
    cnt = 0
    for i in range(0, C):
        for j in range(0, C):
            x[i+j*C] = i
            y[i+j*C] = j
            z[i+j*C] = (data1[i] - data1[j])
            if i > j:
                d[cnt] = z[i+j*C]
                cnt += 1

    return x,y,z,d

def calculatecompareelementdiff2D(data1,data2):
    C = np.shape(data1)[0]
    x = np.zeros(C*C)
    y = np.zeros(C*C)
    z = np.zeros(C*C)
    d = np.zeros(C*(C+1)/2)
    cnt = 0
    for i in range(0, C):
        for j in range(0, C):

            if i > j:
                x[i+j*C] = i
                y[i+j*C] = j
                z[i+j*C] = (data1[i] - data1[j])

            if j > i:
                x[i+j*C] = i
                y[i+j*C] = j
                z[i+j*C] = -(data2[i] - data2[j])

            if j == i:
                x[i+j*C] = i
                y[i+j*C] = j
                z[i+j*C] = 0.0

            if i > j:
                d[cnt] = z[i+j*C]
                cnt += 1

    return x,y,z,d

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
