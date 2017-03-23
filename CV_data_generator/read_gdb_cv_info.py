import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

f = open('/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/CrossValidation/GDB-09-High-sdev/gdb-09-1.0sdev.dat','r')

sigma = []
for i in f.readlines():
    data = i.split(":")

    mid = int(data[0].split('(')[0])

    sd = float(data[3].split('=')[1])
    sigma.append(sd)
    print(mid,' ',sd)

x = np.sort(np.array(sigma))
x1 = x[:x.shape[0]-2500]
x2 = x[x.shape[0]-2500:]

print(x)

f, axarr = plt.subplots(2)

n, bins, patches = axarr[0].hist(x1, 250, normed=0, facecolor='green', alpha=0.75)
n, bins, patches = axarr[1].hist(x2, 250, normed=0, facecolor='green', alpha=0.75)

axarr[0].grid(True)
axarr[1].grid(True)

axarr[0].set_xlabel('Standard deviation (kcal/mol)')
axarr[0].set_ylabel('Count')

axarr[1].set_xlabel('Standard deviation (kcal/mol)')
axarr[1].set_ylabel('Count')

print('Total: ', x.shape[0], ' Split: ', x.shape[0]-2500)

#plt.suptitle(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])

plt.show()
