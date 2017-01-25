import numpy as np
import graphtools as gt
from os import listdir

dtdirs = ["/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_01/testdata/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_02/testdata/",
         ]



namelist = ["_train.dat","_valid.dat","_test.dat"]

for d in dtdirs:

    tmp = listdir(d)
    files = {("_".join(f.split("_")[:-1])) for f in tmp}
    for i in files:
        allarr = []
        for k in namelist:
            try:
                readarrays = gt.readncdat(d+i+k)
                shapesarr = [x.shape for x in readarrays]
            except FileNotFoundError:
                readarrays = [np.zeros((0,*x[1:])) for x in shapesarr]
            allarr.append(readarrays)

        xyz,typ,E = zip(*allarr)

        xyz = np.concatenate(xyz)
        E = np.concatenate(E)
        typ = typ[0]

        print(xyz)

