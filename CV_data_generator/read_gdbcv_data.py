import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def perc_better(a1,a2):
    N = 0
    for i,j in zip(a1,a2):
        if i == j:
            continue
        elif i < j:
            N = N + 1
        else:
            N = N - 1
    return 100.0*(N / float(a1.shape[0]))

def plot_pie(grid, std, rstd, Na, title):
    pb = perc_better(std, rstd)

    std = std/Na

    f1 = std[(std <= 0.023)].shape[0]
    f2 = std[(std > 0.023) & (std <= 0.05)].shape[0]
    f3 = std[(std > 0.05)].shape[0]

    explode = (0.05, 0.05, 0.05)
    labels = 's <= 0.023', '0.023 < s <= 0.05', 's > 0.05'
    fracs = [f1, f2, f3]

    plt.subplot(grid, aspect=1)
    patches, texts, autotexts = plt.pie(fracs, explode=explode, labels=labels,
                                        autopct='%1.1f%%', shadow=True, startangle=90)
    #for t in texts:
    #    t.set_size('smaller')
    #for t in autotexts:
    #    t.set_size('small')
    plt.title(title + ' (Improv.: ' + "{:.1f}".format(pb) +'%;' + str(int((pb/100.0)*rstd.shape[0])) + ')', bbox={'facecolor': '0.8', 'pad': 5})

wkdir = '/home/jujuman/Research/CrossValidation/GDB-06-High-sdev/'

data = np.loadtxt(wkdir + 'gdb-06-cvsdev_c08f.dat', delimiter=':', dtype=str)
std1 = np.asarray([np.float32(str(i).split("=")[-1].split(" ")[0]) for i in data[:,3]])
Na1  = np.asarray([np.float32(str(i).split("(")[-1].split(")")[0]) for i in data[:,0]])

data = np.loadtxt(wkdir + 'gdb-06-cvsdev_c08f09bad.dat', delimiter=':', dtype=str)
std2 = np.asarray([np.float32(str(i).split("=")[-1].split(" ")[0]) for i in data[:,3]])
Na2  = np.asarray([np.float32(str(i).split("(")[-1].split(")")[0]) for i in data[:,0]])

data = np.loadtxt(wkdir + 'gdb-06-cvsdev_c08f09dd.dat', delimiter=':', dtype=str)
std3 = np.asarray([np.float32(str(i).split("=")[-1].split(" ")[0]) for i in data[:,3]])
Na3  = np.asarray([np.float32(str(i).split("(")[-1].split(")")[0]) for i in data[:,0]])

data = np.loadtxt(wkdir + 'gdb-06-cvsdev_c08f09div.dat', delimiter=':', dtype=str)
std4 = np.asarray([np.float32(str(i).split("=")[-1].split(" ")[0]) for i in data[:,3]])
Na4  = np.asarray([np.float32(str(i).split("(")[-1].split(")")[0]) for i in data[:,0]])

the_grid = GridSpec(2, 2)

plot_pie(the_grid[0, 0],std1,std1, Na1,'Original')
plot_pie(the_grid[0, 1],std2,std1, Na2,'Worst 2500')
plot_pie(the_grid[1, 0],std3,std1, Na3,'Diverse dimers')
plot_pie(the_grid[1, 1],std4,std1, Na4,'Diverse 2500')

plt.suptitle("Comparison of strategies to sample molecules from GDB-06 using cross validation\n(Total molecules: " + str(std1.shape[0]) + ") [s=std. dev.; kcal/mol/atom]",fontsize=16)
plt.show()
