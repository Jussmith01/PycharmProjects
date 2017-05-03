import numpy as np
import matplotlib.pyplot as plt; plt.rcdefaults()
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

def plot_pie(grid, std, rstd, Na, title, stride, shift):
    pb = perc_better(std, rstd)

    std = std/Na

    f1 = std[(std <= 0.023)].shape[0]
    f2 = std[(std > 0.023) & (std <= 0.05)].shape[0]
    f3 = std[(std > 0.05)].shape[0]

    labels = ['s <= 0.023', '0.023 < s <= 0.05', 's > 0.05']
    fracs = [f1, f2, f3]

    x_pos = np.arange(len(labels))
    print(x_pos)
    print(fracs)

    from matplotlib import cm
    cs = cm.Paired([0.2,0.9,0.42])

    plt.subplot(grid,aspect=0.002)
    plt.bar(x_pos, fracs, align='center', alpha=0.5)
    plt.xticks(x_pos, labels)

    plt.title(title + ' (Improv.: ' + "{:.1f}".format(pb) +'%;' + str(int((pb/100.0)*rstd.shape[0])) + ')', bbox={'facecolor': '0.8', 'pad': 5})

def plot_bar(grid, s, name_list, width, color, index, expr_list, label):

    fracs = []
    for e in expr_list:
        fracs.append( s[eval(e)].shape[0] )
    x_pos = np.arange(len(fracs))
    fracs = np.array(fracs)

    plt.subplot(grid)
    rects = plt.bar(x_pos+index*width, fracs, width, color=color, label=label, align='center', alpha=0.5)
    plt.xticks(x_pos,expr_list)
    #plt.title(title+' '+expression+'')
    return rects, fracs.min(),fracs.max()

N = '09'

wkdir = '/home/jujuman/Research/CrossValidation/GDB-' + N + '-High-sdev/'

data = np.loadtxt(wkdir + 'gdb-' + N + '-cvsdev_c08f.dat', delimiter=':', dtype=str)
Na1  = np.asarray([np.float32(str(i).split("(")[-1].split(")")[0]) for i in data[:,0]])
std1 = np.asarray([np.float32(str(i).split("=")[-1].split(" ")[0]) for i in data[:,3]])/Na1

data = np.loadtxt(wkdir + 'gdb-' + N + '-cvsdev_c08f09bad.dat', delimiter=':', dtype=str)
Na2  = np.asarray([np.float32(str(i).split("(")[-1].split(")")[0]) for i in data[:,0]])
std2 = np.asarray([np.float32(str(i).split("=")[-1].split(" ")[0]) for i in data[:,3]])/Na2

data = np.loadtxt(wkdir + 'gdb-' + N + '-cvsdev_c08f09dd.dat', delimiter=':', dtype=str)
Na3  = np.asarray([np.float32(str(i).split("(")[-1].split(")")[0]) for i in data[:,0]])
std3 = np.asarray([np.float32(str(i).split("=")[-1].split(" ")[0]) for i in data[:,3]])/Na3

data = np.loadtxt(wkdir + 'gdb-' + N + '-cvsdev_c08f09div.dat', delimiter=':', dtype=str)
Na4  = np.asarray([np.float32(str(i).split("(")[-1].split(")")[0]) for i in data[:,0]])
std4 = np.asarray([np.float32(str(i).split("=")[-1].split(" ")[0]) for i in data[:,3]])/Na4

the_grid = GridSpec(1, 1)

name_list = ['Original', 'Worst 2500', 'Diverse dimers', 'Diverse 2500']
expr_list = ['s <= 0.023', '(s > 0.023) & (s <= 0.05)', 's > 0.05']
data_list = [std1, std2, std3, std4]

rect1, min1, max1 = plot_bar(the_grid[0, 0],data_list[0], name_list, 0.2, 'red', -1.5, expr_list, 'Original')
rect2, min2, max2 = plot_bar(the_grid[0, 0],data_list[1], name_list, 0.2, 'green', -0.5, expr_list, '2500 Worst')
rect3, min3, max3 = plot_bar(the_grid[0, 0],data_list[3], name_list, 0.2, 'blue', 0.5, expr_list, '2500 Diverse')
rect4, min4, max4 = plot_bar(the_grid[0, 0],data_list[2], name_list, 0.2, 'yellow', 1.5, expr_list, 'Diverse dimers')

mins = np.array([min1,min2,min3,min4,])
maxs = np.array([max1,max2,max3,max4,])

dy = np.abs(mins.min() - maxs.max())
plt.ylim([mins.min() - dy * 0.1, maxs.max() + dy * 0.1])

def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = 100.0 * (rect.get_height() / std1.shape[0])
        plt.text(rect.get_x() + rect.get_width()/2., 1.01*rect.get_height(),
                "{:.1f}".format(height) + '%',
                ha='center', va='bottom')

autolabel(rect1)
autolabel(rect2)
autolabel(rect3)
autolabel(rect4)

plt.legend(loc='upper right', shadow=True)

plt.suptitle('Comparison of strategies to sample molecules from GDB-' + N + ' using cross validation\n(Total molecules: ' + str(std1.shape[0]) + ') [s=std. dev.; kcal/mol/atom]',fontsize=16)
plt.show()
