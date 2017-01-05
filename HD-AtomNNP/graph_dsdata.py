import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.ticker import LogLocator, FormatStrFormatter
import matplotlib
#------------------Training-------------------------
data = np.array([[1.493359, 2.065477, 2.103107, 3.205436, 860000.00],
                 [1.560947, 2.066859, 2.128276, 3.163914, 860000.00],
                 [1.437753, 2.020403, 2.094926, 3.018101, 860000.00],
                 [1.599113, 2.062258, 2.144063, 3.113774, 860000.00],
                 [1.390251, 1.730940, 1.800485, 2.680970, 1720000.00],
                 [1.289571, 1.684550, 1.770406, 2.826720, 1720000.00],
                 [1.441021, 1.795275, 1.828275, 2.805933, 1720000.00],
                 [1.183598, 1.428861, 1.453163, 2.283762, 4300000.00],
                 [1.173134, 1.418987, 1.447598, 2.407515, 4300000.00],
                 [1.151351, 1.400439, 1.436869, 2.457207, 4300000.00],
                 [1.197017, 1.418429, 1.455871, 2.373005, 4300000.00],
                 [1.170853, 1.323417, 1.335303, 2.220478, 8600000.00],
                 [1.195701, 1.332740, 1.360572, 2.215108, 8600000.00],
                 [1.087846, 1.200135, 1.214918, 2.056411, 12900000.00],
                 [1.157926, 1.278807, 1.275528, 1.914648, 17200000.00]])

#data = np.log10(data)

x1 = np.mean(data[0:4, 4])
x2 = np.mean(data[4:7, 4])
x3 = np.mean(data[7:11,4])
x4 = np.mean(data[11:13,4])
x5 = np.mean(data[13, 4])
x6 = np.mean(data[14, 4])

x = np.array([x1,x2,x3,x4,x5,x6])
print(x)

tr_y = np.array([np.mean(data[0:4, 0])
                ,np.mean(data[4:7, 0])
                ,np.mean(data[7:11,0])
                ,np.mean(data[11:13,0])
                ,np.mean(data[13, 0])
                ,np.mean(data[14, 0])])
print(tr_y)

vd_y = np.array([np.mean(data[0:4, 1])
                ,np.mean(data[4:7, 1])
                ,np.mean(data[7:11,1])
                ,np.mean(data[11:13,1])
                ,np.mean(data[13, 1])
                ,np.mean(data[14, 1])])
print(vd_y)

te_y = np.array([np.mean(data[0:4, 2])
                ,np.mean(data[4:7, 2])
                ,np.mean(data[7:11,2])
                ,np.mean(data[11:13,2])
                ,np.mean(data[13, 2])
                ,np.mean(data[14, 2])])
print(te_y)

te10_y = np.array([np.mean(data[0:4, 3])
                  ,np.mean(data[4:7, 3])
                  ,np.mean(data[7:11,3])
                  ,np.mean(data[11:13,3])
                  ,np.mean(data[13, 3])
                  ,np.mean(data[14, 3])])
print(te10_y)

#data = np.log10(data[])

# First illustrate basic pyplot interface, using defaults where possible.
fig1, ax1 = plt.subplots()

ax1.loglog(x, tr_y, color='blue', label='Training', linewidth=2)
ax1.loglog(data[:,4], data[:,0], marker='o', color='blue', linewidth=0)

ax1.loglog(x, vd_y, color='red', label='Validation', linewidth=2)
ax1.loglog(data[:,4], data[:,1], marker='o', color='red', linewidth=0)

ax1.loglog(x, te_y, color='green', label='Test', linewidth=2)
ax1.loglog(data[:,4], data[:,2], marker='o', color='green', linewidth=0)

ax1.loglog(x, te10_y, color='orange', label='GDB-10 test', linewidth=2)
ax1.loglog(data[:,4], data[:,3], marker='o', color='orange', linewidth=0)

ax1.set_xticks([1000000, 2000000, 3000000, 6000000, 10000000, 20000000])
ax1.set_yticks([1.0, 1.2, 1.5, 2.0, 2.4, 3.0])

ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

# display 5 minor ticks between each major tick
minorLocator = LogLocator(subs=np.linspace(2,10,6,endpoint=False))
# format the labels (if they're the x values)
majorFormatter = FormatStrFormatter('%1.1f')

# for no labels use default NullFormatter
ax1.yaxis.set_minor_locator(minorLocator)

# or if you want to see some constrained floats
# this may make things busy!
ax1.yaxis.set_minor_formatter(majorFormatter)

def myticks(x,pos):
    if x == 0: return "$0$"

    exponent = int(np.log10(x))
    coeff = x/10**exponent

    return r"${:2.0f} \times 10^{{ {:2d} }}$".format(coeff,exponent)

ax1.xaxis.set_major_formatter(mtick.FuncFormatter(myticks))

ax1.set_ylabel('Total energy RMSE (kcal/mol)')
ax1.set_xlabel('Number of data points')

ax1.legend(bbox_to_anchor=(0.68, 0.98), loc=2, borderaxespad=0., fontsize=14)

plt.title("ANI Performance with increasing data set size")

plt.grid(True,which="major",ls="-")

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

plt.show()