import numpy as np
import hdnntools as gt
from os import listdir
import os
import pandas as pd
import itertools
import time as tm

path = "/home/jujuman/Research/test_data2.h5"

dtdirs = [#"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_dissociation/scans_cc_bonds_dft/single/data/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_01/testdata/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_02/testdata/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_03/testdata/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_04/testdata/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_05/testdata/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_06/testdata/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_07/testdata/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_08/testdata/",
         ]

namelist = ["_train.dat","_valid.dat","_test.dat"]
#namelist = ["_test.dat"]

if os.path.exists(path):
    os.remove(path)
#open an HDF5 for compressed storage.
#Note that if the path exists, it will open whatever is there.
store = pd.HDFStore(path,complib='blosc',complevel=8)

totaltime = 0.0
for d in dtdirs:

    tmp = listdir(d)
    files = list({("_".join(f.split("_")[:-1])) for f in tmp})
    files = sorted(files, key=lambda x: int(x.split('-')[1].split('.')[0]))
    gn = files[0].split("-")[0]

    fcounter = 0
    for n,i in enumerate(files):
        allarr = []

        print(d+i)

        nc = 0
        for k in namelist:
            try:
                _timeloop = tm.time()
                readarrays = gt.readncdat(d+i+k)
                _timeloop2 = (tm.time() - _timeloop)
                totaltime += _timeloop2

                shapesarr = [x.shape for x in readarrays]
                typ = readarrays[1]
            except FileNotFoundError:
                readarrays = [np.zeros((0,*x[1:])) for x in shapesarr]

            ncsub, nat, ndim = readarrays[0].shape
            nc += ncsub
            readarrays[0] = readarrays[0].reshape(ncsub*nat,ndim)

            allarr.append(readarrays)

        xyz, typ, E = zip(*allarr)

        xyz = np.concatenate(xyz).reshape((nc,3 * nat))
        xyz = np.array(xyz,dtype=np.float32)

        E = np.concatenate(E).reshape(nc,1)

        print("Build xyz data frames...")
        cols = [["x" + str(z), "y" + str(z), "z" + str(z)] for z in list(range(nat))]
        cols = [item for sublist in cols for item in sublist]
        cols = [('coordinates',l) for l in cols]
        #cols.append(('energy','E'))
        cols = pd.MultiIndex.from_tuples(cols)  # Notice these are un-named

        # Combine data
        #data = np.append(xyz,E,1)

        df_xyz = pd.DataFrame(xyz, columns=cols)
        df_xyz['energy'] = E
        df_xyz.index.name = 'conformer'

        #print(df_xyz.dtypes)

        print("Store xyzs...")
        store_loc = gn + "/mol" + str(n)
        store.put(store_loc, df_xyz, complib='blosc', complevel=0, format='table')
        store.get_storer(store_loc).attrs.species = typ[0]

        fcounter = fcounter + 1

    print('Total load function time: ' + "{:.4f}".format(totaltime) + 's')

store.close()


'''
# opening file
store = pd.HDFStore(path, complib='blosc', complevel=8)
print(store)
# HDFStore iterates over the names of its contents:


for x in store.get_node(""):
    print("Name:", x._v_name)
    for i in x._v_children:
        store_loc = x._v_name + "/" + i
        data = store.select(store_loc)

        #xyz = np.asarray(xyz)
        xyz = np.asarray(data['coordinates'])
        print(xyz.shape)
        eng = np.asarray(data['energy']).flatten()
        species = store.get_storer(store_loc).attrs.species
        print(xyz)
        print(eng)
        print(species)


    for i in x._v_children:
        print(i)

        # getting an item from the hdf5 file:
        z = store.select(x._v_name + "/" + i,'molecule > 0')

        #print(np.asarray(z))

        # print converting to numpy array because numpy has pretty printing
        maxlen = min(len(z), 10)
        print("Shape:", z.shape)
        print("Value:\n", z[:maxlen])

store.close()
'''

