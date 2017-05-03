import numpy as np
import re

def readxyz (file):
    xyz = []
    typ = []
    Na  = []

    fd = open(file, 'r').read()

    rb = re.compile('((?:[A-Z][a-z]? +?[-+]?\d+?\.\S+? +?[-+]?\d+?\.\S+? +?[-+]?\d+?\.\S+?\s*?(?:\n|$))+)')
    ra = re.compile('([A-Z][a-z]?) +?([-+]?\d+?\.\S+?) +?([-+]?\d+?\.\S+?) +?([-+]?\d+?\.\S+?)\s*?(?:\n|$)')

    s = rb.findall(fd)
    Nc = len(s)
    for i in s:

        c = ra.findall(i)
        Na = len(c)
        for j in c:
            typ.append(j[0])
            xyz.append(j[1])
            xyz.append(j[2])
            xyz.append(j[3])

    xyz = np.asarray(xyz,dtype=np.float32)
    print(xyz.shape)
    print(Nc)
    print(Na)
    xyz = xyz.reshape((Nc,Na,3))

    return xyz,typ[0:Na],Na

def writexyzfile (fn,xyz,typ):
    f = open(fn, 'w')
    N = len(typ)
    for m in xyz:
        f.write(str(N)+'\n')
        f.write('       Comment line\n')
        #print(m)
        for i in range(N):
            x = m[i,0]
            y = m[i,1]
            z = m[i,2]
            f.write(typ[i] + ' ' + "{:.7f}".format(float(x)) + ' ' + "{:.7f}".format(y) + ' ' + "{:.7f}".format(z) + '\n')
    f.close()

file_old = '/home/jujuman/Downloads/mdcrd.xyz'
file_new = '/home/jujuman/Python/PycharmProjects/HD-AtomNNP/ANI_ASE_Programs/mdcrd_new.xyz'

data = readxyz(file_old)
writexyzfile(file_new,data[0],data[1])