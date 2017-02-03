__author__ = 'jujuman'

import numpy as np
#import statsmodels.api as sm
import re
import os.path
import time as tm
import pandas as pd

hatokcal = 627.509469

convert = hatokcal  # Ha to Kcal/mol

def convertatomicnumber(X):
    X = int(X)
    if X == 1:
        return 'H'
    elif X == 6:
        return 'C'
    elif X == 7:
        return 'N'
    elif X == 8:
        return 'O'

def readxyz (file):
    xyz = []
    typ = []
    Na  = []

    fd = open(file, 'r').read()

    #rb = re.compile('\s*?\n?\s*?(\d+?)\s*?\n((?:\s*?[A-Z][a-z]?.+(?:\n|))+)')
    rb = re.compile('(\d+)[\s\S]+?(?=[A-Z])((?:\s*?[A-Z][a-z]?\s+[-+]?\d+?\.\d+?\s+?[-+]?\d+?\.\d+?\s+?[-+]?\d+?\.\d+?\s.+(?:\n|))+)')
    ra = re.compile('([A-Z][a-z]?)\s+?(\S+?)\s+?(\S+?)\s+?(\S+)')

    s = rb.findall(fd)

    for i in s:
        Na.append(int(i[0]))
        atm = ra.findall(i[1])

        print(atm)

        ntyp = []
        nxyz = []
        for j in range(0, int(i[0])):
            ntyp.append(atm[j][0])
            nxyz.append(float(atm[j][1]))
            nxyz.append(float(atm[j][2]))
            nxyz.append(float(atm[j][3]))

        xyz.append(nxyz)
        typ.append(ntyp)

    xyz = np.asarray(xyz,dtype=np.float32)
    xyz = xyz.reshape((xyz.shape[0],len(typ[0]),3))

    return xyz,typ[0],Na

def writexyzfile (fn,xyz,typ):
    f = open(fn, 'w')
    f.write('\n')
    N = len(typ)

    f.write(str(N) + '\n')
    for i in range(N):
        x = xyz[i,0]
        y = xyz[i,1]
        z = xyz[i,2]
        f.write(typ[i] + ' ' + "{:.7f}".format(x) + ' ' + "{:.7f}".format(y) + ' ' + "{:.7f}".format(z) + '\n')
    f.write('\n')
    f.close()

def readncdatwforce (file,N = 0):
    xyz = []
    typ = []
    frc = []
    Eact = []

    readf = False

    if os.path.isfile(file):
        readf = True

        fd = open(file, 'r')

        fd.readline()
        fd.readline()

        types=fd.readline().split(",")

        Na = int(types[0])
        typ = types[1:Na+1]

        cnt = 0

        for i in fd.readlines():
            cnt += 1
            sd = i.strip().split(",")
            #print(sd)
            xyz.append(list(map(float, sd[0:3*Na])))
            Eact.append(float(sd[3*Na]))
            frc.append(list(map(float,sd[3*Na+1:2*3*Na+1])))
            if cnt >= N and N > 0:
                break


    return xyz,frc,typ,Eact,readf

def readncdat (file,type = np.float):

    if os.path.isfile(file):

        fd = open(file, 'r')

        fd.readline()
        Nconf = int(fd.readline())

        types=fd.readline().split(",")
        Na = int(types[0])
        spc = np.asarray(types[1:Na+1])

        if Nconf > 0:
            #data = np.asarray(pd.read_csv(fd, dtype=type,header=None),dtype=type)
            data = np.loadtxt(fd, delimiter=',',usecols=range(Na*3+1),dtype=type)
        else:
            data = np.empty((0,Na*3+1))

        # Pandas sucks!! Have to copy memory else weird stuff happens
        #print(data.shape)
        #print('NS: ', Nconf,' ',Na,' ',3)

        if Nconf != data.shape[0]:
            print('Error: shapes dont match for file: ',file)
            print(Nconf, '!=', data.shape[0])
            exit(1)

        xyz    = data[:Nconf,0:3*Na].reshape(Nconf,Na,3).copy()
        energy = data[:Nconf,3*Na].flatten()

        return xyz,spc,energy
    else :
        raise (FileNotFoundError("File not found: " + file))

def writencdat (file,xyz,spec,energy,comment=''):
    f = open(file,'w')
    f.write(comment+'\n')
    f.write(str(xyz.shape[0])+'\n')
    f.write(str(len(spec)) + ',')
    for i in spec:
        f.write(i + ',')
    f.write('\n')

    for m,e in zip(xyz,energy):
        for a in m:
            for c in a:
                f.write("{:.7f}".format(c) + ',')
        f.write("{:.7f}".format(e) + ',\n')
    f.close()



def readg09trajdat (file,type):
    xyz = []
    typ = []
    Enr = []

    fd = open(file, 'r').read()

    # Get energies
    rE = re.compile('SCF Done:\s+?E\(\S+\)\s+?=\s+?([+,-]?\d+?\.\d+?E?[+,-]?\d+?)\s')
    s = rE.findall(fd)

    for i in s:
        Enr.append(float(i))

    # Get coords
    rB = re.compile('Input orientation:([\S\s]+?)(?=Distance|Rotational)')
    b = rB.findall(fd)

    for i in b:
        rX = re.compile('\s+?\d+?\s+?(\d+?)\s+?\d+?\s+?([+,-]?\d+?\.\d+?)\s+?([+,-]?\d+?\.\d+?)\s+?([+,-]?\d+?\.\d+?)\s')
        c = rX.findall(i)

        t_xyz = []
        t_typ = []

        for j in c:
            t_typ.append(convertatomicnumber(j[0]))
            t_xyz.append(float(j[1]))
            t_xyz.append(float(j[2]))
            t_xyz.append(float(j[3]))

        typ.append(t_typ)
        xyz.append(t_xyz)

    #typ.pop(len(typ)-1)
    #xyz.pop(len(xyz)-1)
    #typ.pop(0)
    #xyz.pop(0)

    Enr.pop(0)
    #typ.pop(0)
    #xyz.pop(0)

    #print (len(typ[0]))
    xyz = np.asarray(xyz,dtype=type)
    xyz = xyz.reshape((xyz.shape[0],len(typ[0]),3))

    Enr = np.asarray(Enr,dtype=type)
    return xyz,typ,Enr

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

def generatedmat(crds,Na):
    dmat = []

    for i in range(0,Na):
        for j in range(i+1, Na):
            dmat.append(((crds[i*3] - crds[j*3])**2+(crds[i*3+1] - crds[j*3+1])**2+(crds[i*3+2] - crds[j*3+2])**2)**0.5)

    return np.array(dmat)

def generatedmat(crds,Na):
    dmat = []

    for i in range(0,Na):
        for j in range(i+1, Na):
            dmat.append(((crds[i*3] - crds[j*3])**2+(crds[i*3+1] - crds[j*3+1])**2+(crds[i*3+2] - crds[j*3+2])**2)**0.5)

    return dmat

def generatedmats(crds,Na):
    print(crds.shape[0])
    dmat = np.zeros([crds.shape[0],int((Na*(Na+1))/2)],np.float)

    for a in dmat:
        count = 0
        for i in range(0,Na):
            for j in range(i+1, Na):
                a[count] = ((crds[i*3] - crds[j*3])**2+(crds[i*3+1] - crds[j*3+1])**2+(crds[i*3+2] - crds[j*3+2])**2)**0.5
                count += 1

    return dmat

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
    d = np.zeros(int(C*(C-1)/2))
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
