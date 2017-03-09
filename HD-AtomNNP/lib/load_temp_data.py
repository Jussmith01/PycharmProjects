import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('/home/jujuman/Python/PycharmProjects/HD-AtomNNP/ANI_ASE_Programs/temp.dat',dtype=np.float32,delimiter=' ')
data = np.loadtxt('/home/jujuman/Research/CrossValidation/MD_CV/md-peptide-cv.dat',dtype=np.float32,delimiter=' ')

print(data)

#t = (data[:,0] * 0.25) / 1000.0
t = data[:,0]

plt.plot(t,data[:,1],color='blue',label='E0')
plt.plot(t,data[:,2],color='red',label='E1')
plt.plot(t,data[:,3],color='green',label='E2')
plt.plot(t,data[:,4],color='orange',label='E3')
plt.plot(t,data[:,5],color='black',label='E4')

plt.ylabel('E (eV)')
plt.xlabel('t (ps)')
plt.legend(bbox_to_anchor=(0.2, 0.98), loc=2, borderaxespad=0., fontsize=14)

plt.show()

plt.plot(t,data[:,6],color='black',label='E4')
plt.show()

