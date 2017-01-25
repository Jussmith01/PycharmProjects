import numpy as np
import graphtools as gt
from os import listdir
import os
import pandas as pd
import itertools

path = "store_test.h5"

dtdirs = ["/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_01/testdata/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_02/testdata/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_03/testdata/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_04/testdata/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_05/testdata/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_06/testdata/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_07/testdata/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_08/testdata/",
         ]



namelist = ["_train.dat","_valid.dat","_test.dat"]

if os.path.exists(path):
    os.remove(path)
#open an HDF5 for compressed storage.
#Note that if the path exists, it will open whatever is there.
store = pd.HDFStore(path,complib='zlib',complevel=9)

for d in dtdirs:

    tmp = listdir(d)
    files = {("_".join(f.split("_")[:-1])) for f in tmp}
    for i in files:
        print(d+i)
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

        #mi = pd.MultiIndex.from_product(
        #    [np.arange(xyz.shape[0]), np.arange(xyz.shape[1])],
        #    names=["conform"])

        #print (xyz.shape)

        clabel = [string + str(num) for num,string in itertools.product(range(xyz.shape[1]),["x", "y", "z"])]

        df_xyz = pd.DataFrame(xyz.reshape(xyz.shape[0], xyz.shape[2] * xyz.shape[1]), columns=clabel)
        df_typ = pd.DataFrame(typ, columns=["typ"])
        df_E   = pd.DataFrame(E, columns=["E"])

        gn = i.replace("-","m")
        store.put(gn+"/xyz",df_xyz,complevel=9)
        store.put(gn+"/typ",df_typ,complevel=9)
        store.put(gn+"/E",  df_E,  complevel=9)

store.close()


'''
# opening file
store = pd.HDFStore(path, complib='zlib', complevel=9)
print(store)
# HDFStore iterates over the names of its contents:

for x in store.get_node(""):
    print("Name:", x._v_name)

    for i in x._v_children:
        print(i)

        # getting an item from the hdf5 file:
        z = store.select(x._v_name + "/" + i)

        #print(np.asarray(z))

        # print converting to numpy array because numpy has pretty printing
        maxlen = min(len(z), 100)
        print("Shape:", z.shape)
        print("Value:\n", z[:maxlen])
store.close()
'''

