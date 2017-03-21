import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

exp = np.loadtxt('/home/jujuman/Desktop/exp_295.dat',delimiter=' ')
ani = np.loadtxt('/home/jujuman/Desktop/rdf_OO_scaled.dat',delimiter=' ')

print(exp[:,0])

plt.plot(exp[:,0],exp[:,1],color='black',label='Exp.',linewidth=4)
plt.plot(ani[:,0],ani[:,1],'r--',color='red',label='ANI-1.1',linewidth=4)

plt.title('O-O radial distribution function')

plt.xlabel('Distance $\AA$')
plt.ylabel('g(r)')

plt.legend(bbox_to_anchor=(0.6, 0.98), loc=2, borderaxespad=0., fontsize=14)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)


plt.show()