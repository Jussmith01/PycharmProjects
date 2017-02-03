from os import listdir
import numpy as np
import re
import hdnntools as hdt
import matplotlib.pyplot as plt
import matplotlib as mpl

def get_rcbd_input_crds (file):
    xyz = []
    typ = []
    Na  = []

    fd = open(file, 'r').read()

    ra = re.compile('\$coordinates\s*?\n([\s\S]+?)\&')
    rb = re.compile('([A-Z][a-z]?)\s+?[A-Z][a-z]?\s+?(\S+?)\s+?(\S+?)\s+?(\S+?)\s+?')

    s = ra.findall(fd)

    for i in s:
        atm = rb.findall(i)

        ntyp = []
        nxyz = []
        for j in atm:
            ntyp.append(j[0])
            nxyz.append(float(j[1]))
            nxyz.append(float(j[2]))
            nxyz.append(float(j[3]))

        xyz.append(nxyz)
        typ.append(ntyp)

    xyz = np.asarray(xyz,dtype=np.float32)
    xyz = xyz.reshape((xyz.shape[0],len(typ[0]),3))

    return xyz,typ

def compute_dist(xyz,spc,a1,a2):
    dists = []

    for i in range(0,xyz.shape[0]):
        for j in range(i+1, xyz.shape[0]):
            if (spc[i] == a1 and spc[j] == a2) or (spc[i] == a2 and spc[j] == a1):
                dists.append(np.linalg.norm(xyz[i] - xyz[j]))

    return dists

comment = 'Optimized at wb97x/6-31g*'

dtdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/gdb-ani-inputs/'
files = listdir(dtdir)

S = ['01','02','03','04','05','06','07','08']

for set in S:
    for i in files:
        if 's' + set in i:
            print('File: ', i)
            xyz, spc = get_rcbd_input_crds(dtdir+i)

            energy = np.zeros(xyz.shape[0],dtype=np.float32)

            name = 'data/minimized/GDB-' + set + '/' + i.split('.')[0] + ".dat"
            hdt.writencdat(name, xyz, spc[0], energy, comment)

    #dist.extend(compute_dist(xyz[0],spc[0],'C','C'))

'''
f= open("ccdist.dat",'w')
for i in dist:
    f.write("{:.7f}".format(i)+'\n')
f.close()

fig, axes = plt.subplots(nrows=1, ncols=1)

axes.set_title("All Minimized CC distances")
axes.set_ylabel('Normalized distant count')
axes.set_xlabel('Distance ($\AA$)')

axes.hist(dist, 2000, color='blue', normed=True, label='CC', linewidth=2, range=[0.5, 5.0])

# -----
# PLOT
# -----
plt.show()
'''