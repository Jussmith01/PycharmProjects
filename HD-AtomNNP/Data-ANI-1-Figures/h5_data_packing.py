import numpy as np
import hdnntools as gt
from os import listdir
import os
import pandas as pd
import itertools
import time as tm

path = "test_data.h5"

dtdirs = [#"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_dissociation/scans_cc_bonds_dft/single/data/",
          "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_01/testdata/",
          #"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_02/testdata/",
          #"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_03/testdata/",
          #"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_04/testdata/",
          #"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_05/testdata/",
          #"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_06/testdata/",
          #"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_07/testdata/",
          #"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_08/testdata/",
         ]

#namelist = ["_train.dat"]
namelist = ["_test.dat"]

if os.path.exists(path):
    os.remove(path)
#open an HDF5 for compressed storage.
#Note that if the path exists, it will open whatever is there.
store = pd.HDFStore(path,complib='blosc',complevel=0)

totaltime = 0.0

for d in dtdirs:

    tmp = listdir(d)
    files = list({("_".join(f.split("_")[:-1])) for f in tmp})

    files = sorted(files, key=lambda x: int(x.split('-')[1].split('.')[0]))

    print(files)

    gn = files[0].split("-")[0]

    Natoms = []
    Nconfs = []
    typarr = []

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
                #print('Computation complete. Time: ' + "{:.4f}".format(_timeloop2) + 'ms')

                shapesarr = [x.shape for x in readarrays]
                typ = readarrays[1]
            except FileNotFoundError:
                readarrays = [np.zeros((0,*x[1:])) for x in shapesarr]

            ncsub, nat, ndim = readarrays[0].shape
            nc += ncsub
            readarrays[0] = readarrays[0].reshape(ncsub*nat,ndim)

            allarr.append(readarrays)

        xyz, typ, E = zip(*allarr)

        xyz = np.array(xyz)
        xyz = xyz.reshape((xyz.shape[0]*xyz.shape[1],xyz.shape[2]))

        E =   np.array(E)
        typ = np.array(typarr)

        print("Build data frames...")
        xyz_idx =  pd.MultiIndex.from_product([list([fcounter]),list(np.arange(nc)), list(np.arange(nat))], names=["molecule", "conformer", "atom"])
        df_xyz = pd.DataFrame(xyz, index=xyz_idx, columns=["x", "y", "z"])

        if fcounter == 0:
            store.put(gn + "/xyz", df_xyz, complib='blosc', complevel=0, format='table')
        else:
            store.append(gn + "/xyz",df_xyz)



        #df_typ = pd.DataFrame(typ, columns=["species"])
        #df_E   = pd.DataFrame(E  , columns=["energy"])

        print("Store xyzs...")
        #store.put(gn + "/species", df_typ, complib='blosc', complevel=0, format='table')
        #store.put(gn + "/energy",  df_E,   complib='blosc', complevel=0, format='table')
        #else:
        #    s_xyz = store.select("/"+gn+"/xyz")
        #    s_spc = store.select("/"+gn+"/species")
        #    s_enr = store.select("/"+gn+"/energy")

#            print("s_xyz: ", s_xyz)
        fcounter = fcounter + 1
        print(store)

    print('Total load function time: ' + "{:.4f}".format(totaltime) + 's')

    '''
    molnum = [int(molecule.split('-')[-1]) for molecule in files]

    print ("Tuple 2...")
    typetuples = []
    for i,molecule in enumerate(molnum):
        typetuples.extend([(molecule,atomnum) for atomnum in range(Natoms[i])])
    #print(typetuples)
    mitype = pd.MultiIndex.from_tuples(typetuples, names=["molecule","atom"])
    del typetuples

    print ("Tuple 3...")
    etuples = []
    for i,molecule in enumerate(molnum):
        etuples.extend([(molecule,conform) for conform in range(Nconfs[i])])
    mienergy = pd.MultiIndex.from_tuples(etuples,names=["molecule","conformer"])
    del etuples

    def mkconf(i): return itertools.product(range(Nconfs[i]), range(Natoms[i]))
    def mkmol(k): return itertools.starmap(lambda i, j: (molnum[k], i, j), mkconf(k))
    def mkdata(): return itertools.chain.from_iterable(mkmol(i) for i in range(len(files)))
    x = list(mkdata())
    mi = pd.MultiIndex.from_tuples(x, names=["molecule","conformer","atom"])
    #print("Length: ", len(x), "Sample:", x[:100], sep='\n')

    print("Shape Coords: ",xyz.shape)

    print(mi)

    print("Build data frame...")
    df_xyz = pd.DataFrame(xyz, index=mi , columns=["x", "y", "z"])
    df_typ = pd.DataFrame(typ, index=mitype, columns=["species"])
    df_E   = pd.DataFrame(E,   index=mienergy, columns=["energy"])

    del xyz,typ,E
    print(df_xyz)

    gn = files[0].split("-")[0]
    print("Store xyzs...")
    store.put(gn+"/xyz",    df_xyz,complib='blosc',complevel=0,format='table')
    print("Store species...")
    store.put(gn+"/species",df_typ,complib='blosc',complevel=0,format='table')
    print("Store energies...")
    store.put(gn+"/energy", df_E,  complib='blosc',complevel=0,format='table')
    '''
store.close()



# opening file
store = pd.HDFStore(path, complib='blosc', complevel=0)
print(store)
# HDFStore iterates over the names of its contents:


for x in store.get_node(""):
    print("Name:", x._v_name)

    xyz = store.select(x._v_name + "/" + x.xyz._v_name, 'molecule == 2 & conformer == 1')
    #energy = store.select(x._v_name + "/" + x.energy._v_name, 'molecule == 1  & conformer == 1')
    #species = store.select(x._v_name + "/" + x.species._v_name, 'molecule == 1')

    print (np.asarray(xyz))
    #print (np.asarray(energy).flatten())
    #print (np.asarray(species).flatten())

    '''
    for i in x._v_children:
        print(i)

        # getting an item from the hdf5 file:
        z = store.select(x._v_name + "/" + i,'molecule > 0')

        #print(np.asarray(z))

        # print converting to numpy array because numpy has pretty printing
        maxlen = min(len(z), 10)
        print("Shape:", z.shape)
        print("Value:\n", z[:maxlen])
    '''
store.close()

