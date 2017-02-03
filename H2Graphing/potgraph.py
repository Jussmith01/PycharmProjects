
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
import graphtools as gt

# -----------------------
cmap = mpl.cm.brg
# ------------
# AM1 vs Act
# ------------
user = os.environ['USER']

dir = '/home/jujuman/Python/PycharmProjects/HD-AtomNNP/'
file = 'temp.dat'

data1 = gt.getfltsfromfile(dir + file, ' ', [0])
data2 = gt.getfltsfromfile(dir + file, ' ', [1])

print('Datasize: ' + str(data1.shape[0]))

font = {'family' : 'Bitstream Vera Sans',
            'weight' : 'normal',
            'size'   : 14}

plt.rc('font', **font)

#print(data2)

plt.plot(data1, data2, color='black', label='ANI-c08e',linewidth=2)
plt.plot([0,data1.max()], [300.0,300.0], color='red',linewidth=2)

#plt.scatter(data1[:,0], data2[:,1], color='black',linewidth=4)

#plt.title("H-Gly-Pro-Hyp-Gly-Ala-Gly-OH")
#plt.title("Temperature per step of Langevin MD (T=300K, ts=0.5fs, fc=0.01AU)\nPolypeptide Chain: H-Gly-Pro-Hyp-Gly-Ala-Gly-OH")
plt.title("Temperature per step of Langevin MD (T=300K, ts=0.25fs, fc=0.05AU)\n152 water box")
plt.xlabel('Step')
plt.ylabel('T(K)')
plt.legend(bbox_to_anchor=(0.025, 0.975), loc=2, borderaxespad=0.)


# -----
# PLOT
# -----
plt.show()
