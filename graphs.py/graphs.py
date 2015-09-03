__author__ = 'dustintracy'
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Open and read the data file
#infile = open('/home/jujuman/Gits/ForcePredictionNetwork/g09DNNTSData/TestCase6/RawData2/graph_expvact.dat', 'r')
infile = open('/home/jujuman/Gits/ForcePredictionNetwork/g09DNNTSData/UNIH22.5Network/GPU0/testDataz.dat', 'r')

data = []

for line in infile:
    row = line.strip().split(",")
    data.append(row)

# Truncate and convert to numpy array
data_array = np.array(data)

data = data_array[:-1, [0, 2]]
data = np.array(data, dtype=float)
#data = abs(data)
#data_n = data[:, [0, 1]]
#data_n[:, 1] = data[:, 2] - data[:, 1]

d01 = data[:, 0]
d02 = data[:, 1]
#f1 = data[:, 2]

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='2d')
#ax.scatter(d01,d02)
#ax.set_xlabel("d01")
#ax.set_ylabel("d02")
#ax.set_zlabel("f1")
plt.scatter(d01,d02)
plt.xlabel("d01")
plt.ylabel("d02")
plt.show()
