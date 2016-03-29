__author__ = 'jujuman'

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
import random

def length(x,y,z):
    length = np.sqrt(x*x + y*y + z*z)
    return length


def normalizev3(vec):
    N = int(np.shape(vec)[0]) / int(3)

    for i in range(0,int(N)):
        len = length(vec[i*3],vec[i*3+1],vec[i*3+2])
        vec[i*3] = vec[i*3] / len
        vec[i*3+1] = vec[i*3+1] / len
        vec[i*3+2] = vec[i*3+2] / len

    return vec


#******* PARAMETERS *******

N = 200
M = 1
Na = 12
S = 0.00001
Rmax = 0.3

pos = np.array([[-1.360317, -0.297041, 0.322606, 0.052895, 0.306677, 0.270264, 0.087893, 1.355359, -0.767814, 1.058859, -0.720139, -0.036504, 1.885093, -1.568059, -0.273193, -1.622531, -0.714577, -0.653242, -2.090527, 0.480493, 0.577576, -1.419225,-1.085437, 1.078976, 0.287058, 0.698697, 1.278092, -0.520308, 2.120572, -0.476570, 1.029213, 1.740674, -0.829174, 2.611907, -2.317230, -0.491018]])

#pos = np.array([[0.0, 0.5, 0.0, 0.0, -0.5, 0.0, 0.5, 0.0, 0.0]])

#**************************

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

D = 3
cmap = mpl.cm.jet
for l in range(0,M):
    poswk = pos
    for i in range(0,N):
        found = False
        while not found:
            rand = np.random.random((1, Na * 2))

            T = 2.0 * np.pi * rand[0,0:Na]
            Z = 2.0 * rand[0,Na:Na+Na] - 1.0

            upos = np.zeros([D * Na],dtype=np.float)
            for m in range(0,Na):
                upos[3 * m]     = np.sqrt(1.0 - Z[m]*Z[m]) * np.cos(T[m])
                upos[3 * m + 1] = np.sqrt(1.0 - Z[m]*Z[m]) * np.sin(T[m])
                upos[3 * m + 2] = Z[m]

            upos = S * upos # scale random vector to step size
            tpos = [upos+poswk[i]] # Take a random walk

            found = True
            for k in range(0,Na):
                xm = tpos[0][D*k] - poswk[0][D*k]
                ym = tpos[0][D*k+1] - poswk[0][D*k+1]
                zm = tpos[0][D*k+2] - poswk[0][D*k+2]

                #print('LEN:',length(upos[D*k],upos[D*k+1],upos[D*k+2]))

                if np.sqrt( xm*xm + ym*ym + zm*zm ) > Rmax:
                    found = False

        #for n in range(0,Na):
        #    print(' ',tpos[0][D*n],' ',tpos[0][D*n+1],' ',tpos[0][D*n+2])
        #print('\n')

        poswk = np.append(poswk,tpos,axis=0)

    color = random.uniform(0.0,1.0)
    for j in range(0,Na):
        ax.scatter(poswk[:-1,D*j], poswk[:-1,D*j+1], poswk[:-1,D*j+2], color=cmap(color), label='RUN'+str(l),linewidth=3)

ax.set_xlim(-1.0,1.0)
ax.set_ylim(-1.0,1.0)
ax.set_zlim(-1.0,1.0)

plt.show()