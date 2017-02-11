import numpy as np
import hdnntools as gt
from os import listdir
import os
import pandas as pd
import itertools
import time as tm

path = "ethane_CC_disso_data.h5"

dtdirs = ["/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_dissociation/scans_cc_bonds_dft/double/data/",
          #"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_01/testdata/",
          #"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_02/testdata/",
          #"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_03/testdata/",
          #"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_04/testdata/",
          #"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_05/testdata/",
          #"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_06/testdata/",
          #"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_07/testdata/",
          #"/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_08/testdata/",
         ]

namelist = ["_train.dat","_valid.dat","_test.dat"]

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

    allarr = []
    Natoms = []
    Nconfs = []
    typarr = []
    for i in files:
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

        Natoms.append(nat)
        typarr.append(typ)
        Nconfs.append(nc)

    xyz, typ, E = zip(*allarr)

    print('Total load function time: ' + "{:.4f}".format(totaltime) + 's')

    #print(*[thing.shape for thing in xyz],sep='\n')

    #print(typarr[0])
    #print(xyz[0])
    #print(E[0])

    print ("Contruct...")

    xyz = np.concatenate(xyz)
    E = np.concatenate(E)
    typ = np.concatenate(typarr)

    molnum = [int(molecule.split('-')[-1]) for molecule in files]

    '''
    print ("Tuple 1...")
    mituples = []
    for i,molecule in enumerate(molnum):
        mituples.extend([(molecule,conform,atom) for conform,atom in itertools.product(range(Nconfs[i]),range(Natoms[i]))])
    #mituples = [(int(molecule.split("-")[-1]),conform,atom) for molecule,conform,atom in itertools.product(files,Nconfs,Natoms)]
    #print(len(mituples))
    mi = pd.MultiIndex.from_tuples(mituples, names=["molecule","conformer","atom"])
    del mituples
    '''

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
    #pd.MultiIndex.from_product(
    #[np.arange(xyz.shape[0]), np.arange(xyz.shape[1])],
    #names=["conform"])

    #print (xyz.shape)

    #clabel = [string + str(num) for num,string in itertools.product(range(xyz.shape[1]),["x", "y", "z"])]

    print(files)

    def mkconf(i): return itertools.product(range(Nconfs[i]), range(Natoms[i]))
    def mkmol(k): return itertools.starmap(lambda i, j: (molnum[k], i, j), mkconf(k))
    def mkdata(): return itertools.chain.from_iterable(mkmol(i) for i in range(len(files)))
    x = list(mkdata())
    mi = pd.MultiIndex.from_tuples(x, names=["molecule","conformer","atom"])
    #print("Length: ", len(x), "Sample:", x[:100], sep='\n')

    print("Shape Coords: ",xyz.shape)

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

store.close()



# opening file
store = pd.HDFStore(path, complib='zlib', complevel=9)
print(store)
# HDFStore iterates over the names of its contents:

for x in store.get_node(""):
    print("Name:", x._v_name)

    xyz = store.select(x._v_name + "/" + x.xyz._v_name, 'molecule == 1 & conformer == 1')
    energy = store.select(x._v_name + "/" + x.energy._v_name, 'molecule == 1  & conformer == 1')
    species = store.select(x._v_name + "/" + x.species._v_name, 'molecule == 1')

    print (np.asarray(xyz))
    print (np.asarray(energy).flatten())
    print (np.asarray(species).flatten())
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

