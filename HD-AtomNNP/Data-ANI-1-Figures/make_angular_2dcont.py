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
def kde_estimate(x, xmin, xmax, y, ymin, ymax,resolution):
    # Peform the kernel density estimate
    xx, yy = np.mgrid[xmin:xmax:resolution, ymin:ymax:resolution]
    print(xx.shape)
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

def make_polar_plot(axes,data,spec,an,color='black',fraction=1.0):
    n_atoms = spec.shape[0]
    randindices = np.random.permutation(n_atoms)[:int(n_atoms * fraction)]

    plot_data = []
    for i in randindices:
        a = data[i]
        t = spec[i]
        # print(a, ' : ', t)
        if nparray_compare(t, an[0], an[1], an[2]):
            plot_data.append(a)

    # print ('Plot data: ',pld)

    pld = np.vstack(plot_data)

    # pH - 4.00
    da = 0.5 * (pld[:, 0] + pld[:, 1])
    an = 2.0 * pld[:, 2]

    #x = da * np.cos(an)
    #y = da * np.sin(an)

    #xmin, xmax = -3.5, 3.5
    #ymin, ymax = -3.5, 3.5

    #xx, yy, f = kde_estimate(an, 0, 2*3.14159, da, 0.0, da.max(), 200j)

    # axes.set_xlim(xmin, xmax)
    # axes.set_ylim(xmin, xmax)

    # Contourf plot
    #cfset = axes.contourf(xx, yy, f, 100, cmap='jet')
    #plt.colorbar(cfset)
    axes.scatter(an, da, marker='.', color=color, linewidths=1)
    # set_polar_grid(axes)

print('Loading data...')
P = ['04','05','06','07']
an = [6,6,6]

# Creating subplots and axes dynamically
axes = plt.subplot(111, projection='polar')

for p in P:
    print('loading ',p,'...')

    data = np.load('data/angular/GDB-' + p + '_data.npz')['arr_0']
    spec = np.load('data/angular/GDB-' + p + '_spec.npz')['arr_0']

    data2 = np.load('data/minimized/angular/GDB-' + p + '_data.npz')['arr_0']
    spec2 = np.load('data/minimized/angular/GDB-' + p + '_spec.npz')['arr_0']

    make_polar_plot(axes,data,spec,an,'blue',0.2)
    make_polar_plot(axes,data2,spec2,an,'red',1.0)

plt.title(hdt.convertatomicnumber(an[1])+"-"+
          hdt.convertatomicnumber(an[0])+"-"+
          hdt.convertatomicnumber(an[2])+
          " 2d angle and average distance representation",y=1.00)

plt.show()