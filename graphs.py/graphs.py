import numpy as np
import statsmodels.api as sm
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
datall = getfltsfromfile('/home/jujuman/Gits/ForcePredictionNetwork/g09DNNTSData/H2OSpherTest/GPU1/tdatatest.dat', [13, 14, 15])
datahl = getfltsfromfile('/home/jujuman/Gits/ForcePredictionNetwork/g09DNNTSData/H2OSpherTest/GPU1/tdatatest.dat', [29, 30, 31])

data1 = makedatalinear(datall)
data2 = makedatalinear(datahl)

calculatemeansqrdiff(data1, data2)

# Do stats data
resultsom = sm.OLS(data2, sm.add_constant(data1)).fit()
print (resultsom.summary())

plt.scatter(data1, data2)
X_plotom = np.linspace(-0.15, 0.55, 100)
plt.plot(X_plotom, X_plotom * resultsom.params[1] + resultsom.params[0])

# ------------------
# ML network vs Act
# ------------------
# Open and read the data file
dataml = getfltsfromfile('/home/jujuman/Gits/ForcePredictionNetwork/g09DNNTSData/H2OSpherTest/GPU1/graph_expvact.dat', [0, 1])

calculatemeansqrdiff(dataml[:, 0], dataml[:, 1])

# Do stats data
resultsml = sm.OLS(dataml[:, 1], sm.add_constant(dataml[:, 0])).fit()
print (resultsml.summary())

plt.scatter(dataml[:, 0], dataml[:, 1],color='green')
X_plotml = np.linspace(-0.15, 0.55, 100)
plt.plot(X_plotml, X_plotml * resultsml.params[1] + resultsml.params[0])
plt.xlabel("Blue = AM1, Green = MLNN")
plt.ylabel("HF/6-31g*")
plt.title("Scatter plot of AM1/MLNN vs HF/6-31g*")

# -----
# PLOT
# -----
plt.show()
