import re
import numpy as np

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

f = open('/home/jujuman/Downloads/l.out','r').read()

ty = []
ls = []
en = []
cd = []


for i in range(0,10):
	ls.append([])

r1 = re.compile('SCF Done:\s+?E\(.+?\)\s+?=\s+?(.+?)A.U.')
r2 = re.compile('Standard orientation:.+\n.+\n.+\n.+\n.+\n([\S\s]+?)(?:---)')
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

print (ty)

en = np.asarray(en)[1:]
cd = np.asarray(cd).reshape(Nc,Na,3)[1:Nc-1]

print (cd.shape,' ',en.shape)

f = open("ircdata.dat",'w')
f.write("comment\n")
f.write(str(Nc-2)+'\n')
f.write(str(Na)+',')
for j in ty: f.write(j+',')
f.write('\n')
mol = 0
for i in cd:

    for j in i:
        for k in j:
            f.write(k+',')
    f.write(en[i] + ',')
    mol += 1
f.close()