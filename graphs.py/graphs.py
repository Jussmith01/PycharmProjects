import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Open and read the data file
#infile = open('/home/jujuman-home/Gits/ForcePredictionNetwork/g09DNNTSData/functionTrainingSetBuilder/functionTrainingSetBuilder/bin/Release/traindata.dat', 'r')
#infile = open('/home/jujuman-home/Gits/ForcePredictionNetwork/g09DNNTSData/H20631gd-UNI/tdatatestUHFH2O.dat', 'r')
infile = open('/home/jujuman-home/Gits/ForcePredictionNetwork/g09DNNTSData/H20631gd-UNI/GPU1/graph_expvact.dat', 'r')

infile_s = []

for line in infile:
    row = line.strip().split(",")
    infile_s.append(row)

# Truncate and convert to numpy array
indata_f = np.array(infile_s)

#data = indata_f[:-1, [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]]
#data = indata_f[:-1, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]]
data = indata_f[:-1, [0, 1]]
data = np.array(data, dtype=float)

infile2 = open('/home/jujuman-home/Gits/ForcePredictionNetwork/g09DNNTSData/H20631gd-UNI/GPU1/testx.dat', 'r')

infile2_s = []

for line in infile2:
    row = line.strip().split(",")
    infile2_s.append(row)

# Truncate and convert to numpy array
indata2_f = np.array(infile2_s)

#data = indata_f[:-1, [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]]
#data = indata_f[:-1, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]]
data2 = indata2_f[:-1, [0, 9]]
data2 = np.array(data2, dtype=float)

#r01 = data[:, [0, 1, 2]] - data[:, [3, 4, 5]]
#print data[0, [0, 1, 2]]
#print data[0, [3, 4, 5]]
#print np.linalg.norm(r01[0,:]) * r01[0,:]
#r02 = data[:, [0, 1, 2]] - data[:, [6, 7, 8]]
#r12 = data[:, [3, 4, 5]] - data[:, [6, 7, 8]]

#f0 = data[:, [6, 7, 8]] - data[:, [9, 10, 11]]
#print f0[0, :]
#print np.dot(np.linalg.norm(r01[0,:]) * r01[0,:], np.linalg.norm(f0[0,:]) * f0[0,:])

#f1 = data[:, [12, 13, 14]]
#f2 = data[:, [15, 16, 17]]

#d01 = np.divide(r01[:, 0], np.power(np.sqrt(np.power(r01[:, 0],2) + np.power(r01[:, 1],2) + np.power(r01[:, 2],2)), 3))
#d02 = r02[:, 0] * (1.0 / np.power(np.sqrt(r02[:, 0] * r02[:, 0] + r02[:, 1] * r02[:, 1] + r02[:, 2] * r02[:, 2]), 3))
#d03 = r12[:, 0] * (1.0 / np.power(np.sqrt(r12[:, 0] * r12[:, 0] + r12[:, 1] * r12[:, 1] + r12[:, 2] * r12[:, 2]), 3))

#lr01 = np.sqrt(r01[:, 0] * r01[:, 0] + r01[:, 1] * r01[:, 1] + r01[:, 2] * r01[:, 2])
#lf01 = np.sqrt(f0[:, 0] * f0[:, 0] + f0[:, 1] * f0[:, 1] + f0[:, 2] * f0[:, 2])

#np.shape(d01)

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(d01, d02, f0[:, 0])
#ax.set_xlabel("d01")
#ax.set_ylabel("d02")
#ax.set_zlabel("f0")
#plt.scatter(lr01,lf01)
plt.scatter(data[:, 0],data[:, 1],color='red')
plt.scatter(data2[:, 0],data2[:, 1])
plt.xlabel("AM1/ML-AM1")
plt.ylabel("6-31g*")
plt.show()
