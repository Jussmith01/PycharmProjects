__author__ = 'jujuman'

import numpy as np
from scipy import stats as st
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
import hdnntools as hdt
#import seaborn as sns
from scipy.interpolate import spline

def get_distance_data(data, spec, an1, an2, fraction=1.0):
    n_atoms = spec.shape[0]
    randindices = np.random.permutation(n_atoms)[:int(n_atoms * fraction)]

    # print('randint: ',randindices)

    plot_data = []
    for i in randindices:
        a = data[i]
        t = spec[i]
        # print(a, ' : ', t)
        c, d = nparray_compare(a, t, an1, an2)
        if c:
            plot_data.append(d)

    return np.array(plot_data)


def nparray_compare(d,t,center,atom):

    if t[0] == center and t[1] == atom:
        return True,d[0]
    elif  t[0] == center and t[2] == atom:
        return True,d[1]
    else:
        return False,0.00


print('Loading data...')

an1 = 6
an2 = 6

P = ['04','05','06']

pld_1l = []
pld_2l = []

for p in P:
    data1 = np.load('data/angular/GDB-' + p + '_data.npz')['arr_0']
    spec1 = np.load('data/angular/GDB-' + p + '_spec.npz')['arr_0']

    data2 = np.load('data/minimized/angular/GDB-' + p + '_data.npz')['arr_0']
    spec2 = np.load('data/minimized/angular/GDB-' + p + '_spec.npz')['arr_0']

    pld_1l.append(get_distance_data(data1,spec1,an1,an2,0.2))
    pld_2l.append(get_distance_data(data2,spec2,an1,an2,1.0))

pld1 = np.concatenate(pld_1l)
pld2 = np.concatenate(pld_2l)

print("Max 1: ", pld1.max())
print("Max 2: ", pld2.max())

'''
ax1 = sns.distplot(pld1,bins=400,hist=False,color='blue',norm_hist=True)

#print(ax1.pro)

ax2 = sns.distplot(pld2,bins=400,kde=False,color='red',norm_hist=True)

print(ax2.properties())

sns.plt.show()
'''

fig, axes = plt.subplots(nrows=1, ncols=1)

axes.set_title("CC Distances")
axes.set_ylabel('Normalized distant count')
axes.set_xlabel('Distance ($\AA$)')


x1, y1, p1 = axes.hist(pld1, 800, color='blue', normed=True, linewidth=2, range=[1.1, 4.6])

print(x1)
print(y1)

x1 = x1 / x1.max()
y1 = y1 - (y1[1] - y1[0])/2.0
y1 = y1[1:]

x_smooth = np.linspace(y1.min(), y1.max(), 2000)
y_smooth = spline(y1, x1, x_smooth)

axes.plot(x_smooth,y_smooth, label='NMS', linewidth=3)

high = float(max([r.get_height() for r in p1]))
for r in p1:
    #r.set_height(r.get_height() / high)
    r.set_height(0.0)
    axes.add_patch(r)
axes.set_ylim(0, 1)


color2 = 'black'
x2, y2, p2 = axes.hist(pld2, 800, color=color2, edgecolor="none", label='EQL',alpha=0.5, normed=True, linewidth=2, range=[1.1, 4.6])

high = float(max([r.get_height() for r in p2]))
for r in p2:
    r.set_height(r.get_height() / high)
#    axes.add_patch(r)

axes.set_ylim(0, 1)
axes.set_xlim(1.1, 4.6)

plt.legend(bbox_to_anchor=(0.7, 0.8), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()
