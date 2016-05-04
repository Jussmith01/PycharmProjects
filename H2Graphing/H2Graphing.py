import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

# Color map
cmap = mpl.cm.nipy_spectral

data1 = np.array([[1,2,3,4,5],[0.05288,0.07957,0.05016,0.06075,0.06860]],dtype=float)
data2 = np.array([[1,2,3,4,5],[0.22605,0.42500,0.85969,0.80775,1.04419]],dtype=float)

plt.plot(data1[0,:], data1[1,:], color='green',linewidth=3)
plt.scatter(data1[0,:], data1[1,:], color='green', label='Modified SF',linewidth=6)

plt.plot(data2[0,:], data2[1,:], color='red',linewidth=3)
plt.scatter(data2[0,:], data2[1,:], color='red', label='Original SF',linewidth=6)

#plt.title(r'$log_{10}(Cost)$ with increaing GDB set size')
plt.xlabel('GDB-N training set size')
plt.ylabel('RMSE (eV)')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()

# ------------
# AM1 vs Act
# ------------


data1 = np.array([[6,5,4,3,2,1],[0.000122966,0.000468814,0.0083554,0.31195,2641.23,21575.49]],dtype=float)
data2 = np.array([[6,5,4,3,2,1],[0.000304857,0.000222569,0.0041405,1.0874,1840.69,18718.00]],dtype=float)
data3 = np.array([[6,5,4,3,2,1],[6.95115e-05,0.000446427,0.0008363,2.1571,1618.07,19275.83]],dtype=float)
data4 = np.array([[6,5,4,3,2,1],[9.39637e-05,0.000113011,0.0045689,1.0759,666.94,20479.38]],dtype=float)
data5 = np.array([[6,5,4,3,2,1],[7.9669e-05,0.000169404,0.0066724,4.4002,3392.95,26263.40]],dtype=float)
data6 = np.array([[6,5,4,3,2,1],[0.000137245,0.000140665,0.0032928,2.3660,4448.22,26127.42]],dtype=float)
data7 = np.array([[6,5,4,3,2,1],[8.61805e-05,0.000449303,0.0112065,4.4691,2208.78,24690.80]],dtype=float)
data8 = np.array([[6,5,4,3,2,1],[0.000249144,0.000333319,0.0265647,3.7400,2142.00,21745.52]],dtype=float)

data1[1,:] = np.log10(data1[1,:])
data2[1,:] = np.log10(data2[1,:])
data3[1,:] = np.log10(data3[1,:])
data4[1,:] = np.log10(data4[1,:])
data5[1,:] = np.log10(data5[1,:])
data6[1,:] = np.log10(data6[1,:])
data7[1,:] = np.log10(data7[1,:])
data8[1,:] = np.log10(data8[1,:])

# --------------
# Setup 2D Plot
# --------------

plt.plot(data1[0,:], data1[1,:], color=cmap(1/9.0),linewidth=1)
plt.scatter(data1[0,:], data1[1,:], color=cmap(1/9.0), label='$C_7 N_3 H_{15}$',linewidth=4)

plt.plot(data2[0,:], data2[1,:], color=cmap(2/9.0),linewidth=1)
plt.scatter(data2[0,:], data2[1,:], color=cmap(2/9.0), label='$C_6 N_2 O_{2} H_{10}$',linewidth=4)

plt.plot(data3[0,:], data3[1,:], color=cmap(3/9.0),linewidth=1)
plt.scatter(data3[0,:], data3[1,:], color=cmap(3/9.0), label='$C_7 N_2 O_{1} H_{14}$',linewidth=4)

plt.plot(data4[0,:], data4[1,:], color=cmap(4/9.0),linewidth=1)
plt.scatter(data4[0,:], data4[1,:], color=cmap(4/9.0), label='$C_8 O_{2} H_{10}$',linewidth=4)

plt.plot(data5[0,:], data5[1,:], color=cmap(5/9.0),linewidth=1)
plt.scatter(data5[0,:], data5[1,:], color=cmap(5/9.0), label='$C_6 N_2 O_{2} H_{8}$',linewidth=4)

plt.plot(data6[0,:], data6[1,:], color=cmap(6/9.0),linewidth=1)
plt.scatter(data6[0,:], data6[1,:], color=cmap(6/9.0), label='$C_5 N_3 O_{2} H_{7}$',linewidth=4)

plt.plot(data7[0,:], data7[1,:], color=cmap(7/9.0),linewidth=1)
plt.scatter(data7[0,:], data7[1,:], color=cmap(7/9.0), label='$C_7 N_2 O_{1} H_{10}$',linewidth=4)

plt.plot(data8[0,:], data8[1,:], color=cmap(8/9.0),linewidth=1)
plt.scatter(data8[0,:], data8[1,:], color=cmap(8/9.0), label='$C_6 N_2 O_{2} H_{8}$',linewidth=4)

#plt.title(r'$log_{10}(Cost)$ with increaing GDB set size')
plt.xlabel('GDB-N training set size')
plt.ylabel('$log_{10}(Cost)$')
plt.legend(bbox_to_anchor=(0.05, 0.6), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()
