import numpy as np
import pyanitools as pyt
from pyNeuroChem import cachegenerator as cg

saef   = "/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk/sae_6-31gd.dat"
h5file = "/home/jujuman/Research/ANI-DATASET/ani_data_c08e_gdb09divaug.h5"

store_dir = "/home/jujuman/Research/CrossValidation/GDB-09-Retrain-DIV/cache-c08e-"

adl = pyt.anidataloader(h5file)

#adl.split_load(10)
N = 10

train_idx = [[2, 3, 4, 5, 6, 7, 8, 9],
             [0, 1, 4, 5, 6, 7, 8, 9],
             [0, 1, 2, 3, 6, 7, 8, 9],
             [0, 1, 2, 3, 4, 5, 8, 9],
             [0, 1, 2, 3, 4, 5, 6, 7]
             ]

valid_idx = [[0]
            ,[2]
            ,[4]
            ,[6]
            ,[8]
             ]

cachet = [cg('_train', saef, store_dir + str(r) + '/') for r in range(5)]
cachev = [cg('_valid', saef, store_dir + str(r) + '/') for r in range(5)]

for data in adl.getnextdata():
    # Print file
    print('Processing file: ', data['parent'], '/', data['child'])

    # Extract the data
    xyz = data['coordinates']
    erg = data['energies']
    spc = data['species']

    xyz = np.array_split(xyz, N)
    erg = np.array_split(erg, N)

    for i,(t,v) in enumerate(zip(cachet, cachev)):
        xyz_t = np.array(np.concatenate([xyz[j] for j in train_idx[i]]), order='C', dtype=np.float32)
        erg_t = np.array(np.concatenate([erg[j] for j in train_idx[i]]), order='C', dtype=np.float64)

        xyz_v = np.array(np.concatenate([xyz[j] for j in valid_idx[i]]), order='C', dtype=np.float32)
        erg_v = np.array(np.concatenate([erg[j] for j in valid_idx[i]]), order='C', dtype=np.float64)

        print(spc)

        t.insertdata(xyz_t, erg_t, list(spc))
        v.insertdata(xyz_v, erg_v, list(spc))

for t,v in zip(cachet, cachev):
    t.makemetadata()
    v.makemetadata()

adl.cleanup()
