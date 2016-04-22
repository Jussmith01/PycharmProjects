import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# ------------
# AM1 vs Act
# ------------


data1 = np.array([[5,4,3,2,1],[0.0022979,0.70203,21.298,2641.23,21575.49]],dtype=float)
data2 = np.array([[5,4,3,2,1],[0.0053911,0.90914,11.237,1840.69,18718.00]],dtype=float)
data3 = np.array([[5,4,3,2,1],[0.0009889,1.53928,17.519,1618.07,19275.83]],dtype=float)
data4 = np.array([[5,4,3,2,1],[0.0026380,1.25865,18.024,666.94,20479.38]],dtype=float)
data5 = np.array([[5,4,3,2,1],[0.0040166,2.68068,21.360,3392.95,26263.40]],dtype=float)
data6 = np.array([[5,4,3,2,1],[0.0018244,0.46006,9.0621,4448.22,26127.42]],dtype=float)
data7 = np.array([[5,4,3,2,1],[0.0063093,1.27487,18.018,2208.78,24690.80]],dtype=float)
data8 = np.array([[5,4,3,2,1],[0.0070196,1.13667,19.006,2142.00,21745.52]],dtype=float)

data1[1,:] = np.log10(data1[1,:])
data2[1,:] = np.log10(data2[1,:])
data3[1,:] = np.log10(data3[1,:])
data4[1,:] = np.log10(data4[1,:])
data5[1,:] = np.log10(data5[1,:])
data6[1,:] = np.log10(data6[1,:])
data7[1,:] = np.log10(data7[1,:])
data8[1,:] = np.log10(data8[1,:])

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

# Color map
cmap = mpl.cm.nipy_spectral

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
