__author__ = 'jujuman'

from matplotlib import pyplot as plt
import numpy as np

#El = float(5.0e-4)
#Tl = 1.0e-13
#Th = 1.0e-12

#Eh = np.linspace(5.5e-4, 9.0e-2, 70, endpoint=True)
#dE = (Eh/El - 1.0)
#print(dE)
#C = Tl/Th - 1.0
#print(C)
#Z = dE*C
#print(Z)
#P = np.exp(Z)
#print(P)
#plt.plot(dE, P)

#plt.show()

# ------------------------------------------
# Calculate the RMSD given two structures
# ------------------------------------------
def computedataset(x1,x2,pts,eta,Rc,plt,scolor,slabel):

    X = np.linspace(x1, x2, pts, endpoint=True)
    F = np.exp(-eta*(X-0.25)**2.0) * (0.25 * (np.cos((np.pi * X)/Rc) + 1.0))
    plt.plot(X, F, label=slabel, color=scolor)


computedataset(0.0, 10.0, 1000, 0.01, 6.0, plt, 'blue', 'eta = 0.1')
computedataset(0.0, 10.0, 1000, 0.05, 6.0, plt, 'green', 'eta = 0.1')
computedataset(0.0, 10.0, 1000, 0.10, 6.0, plt, 'red', 'eta = 0.2')
computedataset(0.0, 10.0, 1000, 0.30, 6.0, plt, 'orange', 'eta = 0.3')
computedataset(0.0, 10.0, 1000, 0.40, 6.0, plt, 'yellow', 'eta = 0.4')
computedataset(0.0, 10.0, 1000, 5.0, 6.0, plt, 'black', 'eta = 0.5')

plt.show()
