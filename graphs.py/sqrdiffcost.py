import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Open and read the data file
infile = open('/home/jujuman-home/Gits/ForcePredictionNetwork/g09DNNTSData/H20631gd-UNI/tdatatestUHFH2O.dat', 'r')

infile_s = []

for line in infile:
    row = line.strip().split(",")
    infile_s.append(row)

# Truncate and convert to numpy array
indata_f = np.array(infile_s)

data = indata_f[:-1, [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]]

data = np.array(data, dtype=float)



#d01 = np.divide(r01[:, 0], np.power(np.sqrt(np.power(r01[:, 0],2) + np.power(r01[:, 1],2) + np.power(r01[:, 2],2)), 3))
#d02 = r02[:, 0] * (1.0 / np.power(np.sqrt(r02[:, 0] * r02[:, 0] + r02[:, 1] * r02[:, 1] + r02[:, 2] * r02[:, 2]), 3))
#d03 = r12[:, 0] * (1.0 / np.power(np.sqrt(r12[:, 0] * r12[:, 0] + r12[:, 1] * r12[:, 1] + r12[:, 2] * r12[:, 2]), 3))

#lr01 = np.sqrt(r01[:, 0] * r01[:, 0] + r01[:, 1] * r01[:, 1] + r01[:, 2] * r01[:, 2])
#lf01 = np.sqrt(f0[:, 0] * f0[:, 0] + f0[:, 1] * f0[:, 1] + f0[:, 2] * f0[:, 2])

