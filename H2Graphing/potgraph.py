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
        row = line.strip().split(",")
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
#user = 'jujuman'
dir = '/Research/ANN-Test-Data/GDB-11-B3LYP-6-31gd/dnntsgdb11_03/trajdata/'

set = 3
mol = 0

cmap = mpl.cm.jet
for i in range(0,8):
    file = 'gdb11_s0' + str(set) + '-' + str(mol) + '_train.dat_trajdata' + str(i)
    #file = 'fixmolecule-0_train.dat'

    data1 = getfltsfromfile('/home/' + user + dir + file, [0])
    data2 = getfltsfromfile('/home/' + user + dir + file, [2])
    data3 = getfltsfromfile('/home/' + user + dir + file, [3])

    #data1 = data1 - data2
    data1 = data1

    color = i/float(8)
    plt.plot(data1, data3, color=cmap(color), label='Thread '+ str(i),linewidth=1)
    plt.scatter(data1, data3, color=cmap(color), label='Thread '+ str(i),linewidth=3)
    #plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.)
    #plt.show()

plt.show()


for i in range(0,8):
    file = 'gdb11_s0' + str(set) + '-' + str(mol) + '_valid.dat_trajdata' + str(i)
    #file = 'fixmolecule-0_train.dat'

    data1 = getfltsfromfile('/home/' + user + dir + file, [0])
    data2 = getfltsfromfile('/home/' + user + dir + file, [2])
    data3 = getfltsfromfile('/home/' + user + dir + file, [3])

    #data1 = data1 - data2
    data1 = data1

    color = i/float(8)
    plt.plot(data1, data3, color=cmap(color), label='Thread '+ str(i),linewidth=1)
    plt.scatter(data1, data3, color=cmap(color), label='Thread '+ str(i),linewidth=3)
    #plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.)
    #plt.show()

plt.show()

#dir = '/Research/ANN-Test-Data/GDB-11-M062X-6-311Gdd/rddata_test/'
#for i in range(0,8):
#    file = 'fixmolecule-0_train.dat_thread' + str(i)
    #file = 'fixmolecule-0_train.dat'

#    data1 = getfltsfromfile('/home/' + user + dir + file, [2])
#    data2 = getfltsfromfile('/home/' + user + dir + file, [5])
#    data3 = getfltsfromfile('/home/' + user + dir + file, [6])

#    data1 = data1 - data2

#    color = i/float(8)
#    plt.scatter(data1, data3, color=cmap(0.9), label='Thread '+ str(i),linewidth=3)

#data4 = getfltsfromfile('/home/' + user + dir + file2, [1])
#data5 = getfltsfromfile('/home/' + user + dir + file3, [1])

#dir = '/Research/ANN-Test-Data/GDB-11/train4/'
#data6 = getfltsfromfile('/home/' + user + dir + 'graph.dat', [2])
#AM1 = getfltsfromfile('/home/' + user + dir + file, [1])

#dir = '/Research/ANN-Test-Data/GDB-11/train3/'
#data5 = getfltsfromfile('/hom-7.558904531543e/' + user + dir + file, [2])
#PM6 = getfltsfromfile('/home/' + user + dir + file, [1])

#dir = '/Research/ANN-Test-Data/GDB-11/train4/'
#data6 = getfltsfromfile('/home/' + user + dir + file, [2])

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

#print(AM1.shape[0])
#AM1 = AM1 - AM1[AM1.shape[0]-1]
#PM6 = PM6 - PM6[AM1.shape[0]-1]
#data2 = data2 - data2[data2.shape[0]-1]
#data3 = data3 - data3[data2.shape[0]-1]
#data4 = data4 - data4[data2.shape[0]-1]
#data5 = data5 - data5[data2.shape[0]-1]
#data6 = data6 - data6[data2.shape[0]-1]
# --------------
# Setup 2D Plot
# --------------
#plt.plot(data1, AM1, color='black', label='AM1',linewidth=4)
#plt.plot(data1, PM6, color='grey', label='PM6',linewidth=4)
#plt.plot(data1, data4, color='grey', label='AM1',linewidth=3)
#plt.plot(data1, data5, color='black', label='PM6',linewidth=4)
#plt.scatter(data1, data3, color='blue', label='B3LYP/6-31g*',linewidth=3)
#plt.plot(data1, data3, color='red', label='ANN - GDB-3',linewidth=3)
#plt.plot(data1, data6, color='orange', label='ANN - GDB-4',linewidth=3)
#plt.scatter(data1, data6, color='green', label='ANN - up to GDB-4',linewidth=4)
#plt.plot(data1, data3, color='red', label='ANN - up to GDB-3',linewidth=4)
#plt.plot(data1, data4, color='green', label='ANN - up to GDB-4',linewidth=4)

#plt.plot(data1, data4, color='green', label='ANN - only $H_3 C H_2 C H_2 N$',linewidth=4)
#plt.plot(data1, data4, color='orange', label='ANN - up to GDB-3',linewidth=4)
#plt.scatter(data1, data6, color='black', label='ANN GDB-3',linewidth=4)

plt.title('Formic Acid Reaction Scan H-O1 -> H-O2')
#plt.title('SCAN: Formic Acid Energy vs. H-O-H Angle')
plt.xlabel('Reaction Coordinate (Angstroms)')
plt.ylabel('Energy (Hartrees)')
plt.legend(bbox_to_anchor=(0.7, 0.3), loc=2, borderaxespad=0.)


# -----
# PLOT
# -----
#plt.show()
