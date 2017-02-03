__author__ = 'jujuman'

import numpy as np
from scipy import stats as st
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
import hdnntools as hdt

def set_polar_grid(axes):
    # Add polar grid
    circle0 = plt.Circle((0, 0), 0.5, color='black', fill=False)
    circle1 = plt.Circle((0, 0), 1.0, color='black', fill=False)
    circle2 = plt.Circle((0, 0), 2.0, color='black', fill=False)
    circle3 = plt.Circle((0, 0), 2.0, color='black', fill=False)

    axes.add_artist(circle0)
    axes.add_artist(circle1)
    axes.add_artist(circle2)
    axes.add_artist(circle3)

    axes.plot([-3.0, 3.0], [0.0, 0.0], color='black')
    axes.plot([0.0, 0.0], [-3.0, 3.0], color='black')

# Plotting kde using matplotlib and scipy
def kde_estimate(x, xmin, xmax, y, ymin, ymax):
    # Peform the kernel density estimate
    xx, yy = np.mgrid[xmin:xmax:200j, ymin:ymax:200j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)

    return(xx, yy, f)

def nparray_compare(t,center,atom1,atom2):

    if np.array_equal(t, np.array([center, atom1, atom2])):
        return True
    elif np.array_equal(t, np.array([center, atom2, atom1])):
        return True
    else:
        return False

print('Loading data...')
p = '04'
data = np.load('data/angular/GDB-' + p + '_data.npz')['arr_0']
spec = np.load('data/angular/GDB-' + p + '_spec.npz')['arr_0']

an1 = 6
an2 = 6
an3 = 6

print(spec.shape)

n_atoms = spec.shape[0]
fraction = 0.05
randindices = np.random.permutation(n_atoms)[:int(n_atoms*fraction)]

#print('randint: ',randindices)

plot_data = []
for i in randindices:
    a = data[i]
    t = spec[i]
    #print(a, ' : ', t)
    if nparray_compare(t,an1,an2,an3):
        plot_data.append(a)

pld = np.vstack(plot_data)

#print ('Plot data: ',pld)

# Creating subplots and axes dynamically
axes = plt.subplot(111)
#fig.set_size_inches(10.0, 5.0, forward=True)

# pH - 4.00
da = 0.5 * (pld[:,0] + pld[:,1])
#an = 2.0 * pld[:,2]

x = da * np.cos(pld[:,2])
y = da * np.sin(pld[:,2])

xmin, xmax = -3.5,3.5
ymin, ymax = -3.5,3.5

xx, yy, f = kde_estimate(x, xmin, xmax, y, ymin, ymax)

#axes.set_xlim(xmin, xmax)
#axes.set_ylim(xmin, xmax)

# Contourf plot
cfset = axes.contourf(xx, yy, f, 50, cmap='Reds')
axes.scatter(x,da,marker='.',color='blue',linewidths=1)
#set_polar_grid(axes)

plt.title(hdt.convertatomicnumber(an2)+"-"+
          hdt.convertatomicnumber(an1)+"-"+
          hdt.convertatomicnumber(an3)+
          " polar plot of angle and average distance",y=1.08)

plt.show()