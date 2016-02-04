__author__ = 'jujuman'

from matplotlib import pyplot as plt
import numpy as np

El = float(5.0e-4)
Tl = 1.0e-13
Th = 1.0e-12

Eh = np.linspace(5.5e-4, 9.0e-2, 70, endpoint=True)
dE = (Eh/El - 1.0)
print(dE)
C = Tl/Th - 1.0
print(C)
Z = dE*C
print(Z)
P = np.exp(Z)
print(P)
plt.plot(dE, P)

plt.show()
