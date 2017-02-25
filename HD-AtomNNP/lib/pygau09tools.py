import numpy as np
import re

def convertatomicnumber(X):
    X = int(X)
    if X == 1:
        return 'H'
    elif X == 6:
        return 'C'
    elif X == 7:
        return 'N'
    elif X == 8:
        return 'O'

def read_charge (file, type='mulliken'):

    fd = open(file, 'r').read()

    if type == 'mulliken':
        ra = re.compile('(?:Mulliken atomic charges:)\n +?\d\n([\s\S]+?)(?:Sum)')
    elif type == 'esp':
        ra = re.compile('(?:Fitting point charges to electrostatic potential)\n.+?\n.+?\n +?\d\n([\s\S]+?)(?:Charges)')

    rb = re.compile(' +?\d+? +?([A-Z][a-z]?) +?([+-]?\d+?\.\S+?)\n')

    species = []
    charges = []
    block = ra.findall(fd)
    for b in block:
        lines = rb.findall(b)
        for l in lines:
            species.append(l[0])
            charges.append(l[1])

    return charges, species

def read_irc (file):
    f = open(file,'r').read()

    ty = []
    ls = []
    en = []
    cd = []

    for i in range(0, 10):
        ls.append([])

    r1 = re.compile('SCF Done:\s+?E\(.+?\)\s+?=\s+?(.+?)A.U.')
    r2 = re.compile('Input orientation:.+\n.+\n.+\n.+\n.+\n([\S\s]+?)(?:---)')
    r3 = re.compile('\d+?\s+?(\d+?)\s+?\d+?\s+?(\S+?)\s+?(\S+?)\s+?(\S+)')

    s1 = r1.findall(f)
    s2 = r2.findall(f)

    for i in s1:
        en.append(i.strip())

    Nc = 0
    for i in s2:

        sm = r3.findall(i.strip())
        Na = len(sm)

        ty = []
        xyz = []

        for j in sm:
            ty.append(convertatomicnumber(j[0]))
            xyz.append(float(j[1]))
            xyz.append(float(j[2]))
            xyz.append(float(j[3]))

        cd.append(xyz)
        Nc += 1

    en = np.asarray(en, dtype=np.float32)
    cd = np.asarray(cd, dtype=np.float32).reshape(Nc, Na, 3)[:-1]

    return [en, cd, ty[0:Na]]

def read_scan (file):
    f = open(file,'r').read()

    ty = []
    ls = []
    en = []
    cd = []

    for i in range(0, 10):
        ls.append([])

    r0 = re.compile('(?:Optimization completed\.)([\s\S]+?)(?:\*\*\*\*\* Axes restored to original set \*\*\*\*\*)')
    r1 = re.compile('SCF Done:\s+?E\(.+?\)\s+?=\s+?(.+?)A.U.')
    r2 = re.compile('Input orientation:.+\n.+\n.+\n.+\n.+\n([\S\s]+?)(?:---)')
    r3 = re.compile('\d+?\s+?(\d+?)\s+?\d+?\s+?(\S+?)\s+?(\S+?)\s+?(\S+)')

    s0 = r0.findall(f)

    for i in s0:
        en.append(r1.findall(i)[0])
        B = r2.findall(i)

        rx = r3.findall(B[0])
        for i in rx:
            ty.append(convertatomicnumber(i[0]))
            cd.append(np.asarray([i[1],i[2],i[3]], dtype=np.float32))

    en = np.asarray(en, dtype=np.float32)
    cd = np.concatenate(cd).reshape(en.shape[0],int(len(cd)/en.shape[0]),3)

    return [en, cd, ty[0:cd.shape[1]]]