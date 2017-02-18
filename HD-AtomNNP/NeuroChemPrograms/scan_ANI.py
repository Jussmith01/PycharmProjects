__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import numpy as np
import matplotlib.pyplot as plt
import re


# ------------------------------------------
#          Radial Function Cos
# ------------------------------------------
def cutoffcos(X,Rc,Ymax):
    Xt = X.copy()
    for i in range(0,Xt.shape[0]):
        if Xt[i] > Rc:
            Xt[i] = Rc

    return Ymax * (0.5 * (np.cos((np.pi * Xt)/Rc) + 1.0))

# ------------------------------------------
#          Radial Function Cos
# ------------------------------------------
def cutoffxsqr(X,Rc,slope):
    Xt = X.copy()
    for i in range(0,Xt.shape[0]):
        if Xt[i] > Rc:
            Xt[i] = Rc

    return ((Xt - Rc) * (Xt - Rc)) * slope

# ------------------------------------------
#          Radial Function Cos
# ------------------------------------------
def morse(X, De, a, Re):
    return De * np.power ( 1.0 - np.exp( -a * ( X - Re ) ), 2.0 )


def readxyz (file):
    xyz = []
    typ = []
    Na  = []

    fd = open(file, 'r').read()

    #rb = re.compile('\s*?\n?\s*?(\d+?)\s*?\n((?:\s*?[A-Z][a-z]?.+(?:\n|))+)')
    rb = re.compile('(\d+)[\s\S]+?(?=[A-Z])((?:\s*?[A-Z][a-z]?\s+[-+]?\d+?\.\d+?\s+?[-+]?\d+?\.\d+?\s+?[-+]?\d+?\.\d+?\s.+(?:\n|))+)')
    ra = re.compile('([A-Z][a-z]?)\s+?(\S+?)\s+?(\S+?)\s+?(\S+)')

    s = rb.findall(fd)

    for i in s:
        Na.append(int(i[0]))
        atm = ra.findall(i[1])

        print(atm)

        ntyp = []
        nxyz = []
        for j in range(0, int(i[0])):
            ntyp.append(atm[j][0])
            nxyz.append(float(atm[j][1]))
            nxyz.append(float(atm[j][2]))
            nxyz.append(float(atm[j][3]))

        xyz.append(nxyz)
        typ.append(ntyp)

    xyz = np.asarray(xyz,dtype=np.float32)
    xyz = xyz.reshape((xyz.shape[0],len(typ[0]),3))

    return xyz,typ[0],Na

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

xyz,typ,Na = readxyz('/home/jujuman/Scratch/Research/GDB-11-wB97X-6-31gd/dnnts_dissociation/scans_cc_bonds_dft/single/methane.xyz')
xyz = xyz[0]

xyz[2][2] = xyz[2][2] - xyz[0][2]
xyz[3][2] = xyz[3][2] - xyz[0][2]
xyz[4][2] = xyz[4][2] - xyz[0][2]
xyz[0][2] = 0

xyz[5][2] = xyz[5][2] - xyz[1][2]
xyz[6][2] = xyz[6][2] - xyz[1][2]
xyz[7][2] = xyz[7][2] - xyz[1][2]
xyz[1][2] = 0

print (xyz)

# Construct pyNeuroChem class
mol = pync.molecule(cnstfile, saefile, nnfdir, 0)

E1 = []

min = 100
max = 700
inc= 0.001

x = np.arange(0.2, 5.0, inc)

for i in x:

    xyzt = xyz.copy()

    xyzt[0][2] = i
    xyzt[2][2] = xyz[2][2] + i
    xyzt[3][2] = xyz[3][2] + i
    xyzt[4][2] = xyz[4][2] + i

    xyzt[1][2] = -i
    xyzt[5][2] = xyz[5][2] - i
    xyzt[6][2] = xyz[6][2] - i
    xyzt[7][2] = xyz[7][2] - i

    mol.setMolecule(coords=xyzt,types=typ)
    E = mol.energy()

    print(xyzt, ' ', E)

    E1.append(E[0])

E1 = np.array(E1)
E1 = E1 - E1.min()

x = 2 * x

plt.plot (x, E1, color='red',  label='ANI',  linewidth=2)

y = cutoffxsqr(x,1.1,10.0)
plt.plot (x, y, color='blue',  label='SRC',  linewidth=2)

plt.plot (x, E1 + y, color='orange',  label='SRC+ANI',  linewidth=2)

#plt.plot (x, morse(x,1.0,1.9,1.53), color='grey',  label='SRC+ANI',  linewidth=2)

plt.title("C-C Dissociation")

plt.ylabel('E (kcal/mol)')
plt.xlabel('Distance $\AA$')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.,fontsize=16)

plt.show()