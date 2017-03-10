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

def make_polar_cont(axes,data,spec,an1,an2,fraction=1.0,saturation=1.0):
    n_atoms = spec.shape[0]
    randindices = np.random.permutation(n_atoms)[:int(n_atoms * fraction)]

    plot_data1 = []
    plot_data2 = []
    for i in randindices:
        a = data[i]
        t = spec[i]
        # print(a, ' : ', t)
        if nparray_compare(t, an1[0], an1[1], an1[2]):
            plot_data1.append(a)
        elif nparray_compare(t, an2[0], an2[1], an2[2]):
            plot_data2.append(a)
    # print ('Plot data: ',pld)

    pld1 = np.vstack(plot_data1)
    pld2 = np.vstack(plot_data2)

    # Angle 1 data
    da1 = 0.5 * (pld1[:, 0] + pld1[:, 1])
    an1 = pld1[:, 2]

    # Angle 2 data
    da2 = 0.5 * (pld2[:, 0] + pld2[:, 1])
    an2 = 2.0 * 3.14159 - pld2[:, 2]

    # Normalize data
    #if da1.max() > da2.max():



    #print(da1)
    #print(da2)

    #print(an1)
    #print(an2)


    # Combine data
    da = np.concatenate([da1,da2])
    an = np.concatenate([an1,an2])

    # Kernel desity estimate
    xx, yy, f = kde_estimate(an, 0, 2.0 * 3.14159, da, 0.0, da.max(), 200j)

    # Saturate the plot
    maxi = f.max()
    f[f > saturation*maxi] = saturation*maxi

    # Contourf plot
    cfset = axes.contourf(xx, yy, f, 100, cmap='Blues',)
    #plt.colorbar(cfset)


def make_polar_plot(axes,data,spec,Z1,Z2,color1='black',color2='black',fraction=1.0):
    n_atoms = spec.shape[0]
    randindices = np.random.permutation(n_atoms)[:int(n_atoms * fraction)]

    plot_data1 = []
    plot_data2 = []
    for i in randindices:
        a = data[i]
        t = spec[i]
        # print(a, ' : ', t)
        if nparray_compare(t, Z1[0], Z1[1], Z1[2]):
            plot_data1.append(a)
        elif nparray_compare(t, Z2[0], Z2[1], Z2[2]):
            plot_data2.append(a)

    # print ('Plot data: ',pld)

    pld1 = np.vstack(plot_data1)
    pld2 = np.vstack(plot_data2)

    # pH - 4.00
    da1 = 0.5 * (pld1[:, 0] + pld1[:, 1])
    an1 = pld1[:, 2]

    da2 = 0.5 * (pld2[:, 0] + pld2[:, 1])
    an2 = 2.0 * 3.14159 - pld2[:, 2]

    label1 = hdt.convertatomicnumber(an1[1]) + "-"
    label2 = hdt.convertatomicnumber(an1[1])

    #label1 = hdt.convertatomicnumber(an1[1]) + "-" + hdt.convertatomicnumber(an1[0]) + "-" + hdt.convertatomicnumber(an1[2])
    #label2 = hdt.convertatomicnumber(an2[1]) + "-" + hdt.convertatomicnumber(an2[0]) + "-" + hdt.convertatomicnumber(an2[2])

    axes.scatter(an1, da1, marker='.', label=label1, color=color1, linewidths=1)
    axes.scatter(an2, da2, marker='.', label=label2, color=color2, linewidths=1)

print('Loading data...')
P = ['05']
an1 = [1,6,8]
an2 = [1,6,7]

# Creating subplots and axes dynamically
axes = plt.subplot(111, projection='polar')

dir = '/home/jujuman/Python/PycharmProjects/HD-AtomNNP/Data-ANI-1-Figures/'

for p in P:
    print('loading ',p,'...')

    data = np.load(dir + 'data/angular/GDB-' + p + '_data.npz')['arr_0']
    spec = np.load(dir + 'data/angular/GDB-' + p + '_spec.npz')['arr_0']

    data2 = np.load(dir + 'data/minimized/angular/GDB-' + p + '_data.npz')['arr_0']
    spec2 = np.load(dir + 'data/minimized/angular/GDB-' + p + '_spec.npz')['arr_0']

    make_polar_cont(axes,data,spec,an1,an2,0.1,0.05)
    #make_polar_plot(axes, data, spec, an, 'blue', 0.2)

    make_polar_plot(axes,data2,spec2,an1,an2,'red','green',1.0)


plt.legend(bbox_to_anchor=(0.9, 0.975), loc=2, borderaxespad=0.)
plt.title("Polar angle and average distance representation",y=1.10)

plt.show()