import numpy as np

# Store test split
import pyanitools as pyt
from pyNeuroChem import cachegenerator as cg

# Set the HDF5 file containing the data
hdf5files = [('/home/jujuman/Research/ANI-DATASET/RXN1_TNET/rxn_h5_data/ani_benz_rxn_test_1.h5','dbmrxn1'),
             ('/home/jujuman/Research/ANI-DATASET/RXN1_TNET/rxn_h5_data/ani_benz_rxn_test_2.h5','dbmrxn2'),
             ('/home/jujuman/Research/ANI-DATASET/RXN1_TNET/rxn_h5_data/ani_benz_rxn_test_3.h5','dbmrxn3'),
             ('/home/jujuman/Research/ANI-DATASET/RXN1_TNET/rxn_h5_data/ani_benz_rxn_test_4.h5','dbmrxn4'),
             ('/home/jujuman/Research/ANI-DATASET/RXN1_TNET/rxn_h5_data/ani_benz_rxn_test_5.h5','dbmrxn5'),
             ('/home/jujuman/Research/ANI-DATASET/RXN1_TNET/rxn_h5_data/ani_benz_rxn_test_6.h5','dbmrxn6'),
             #'/home/jujuman/Research/ANI-DATASET/h5data/ani-gdb-c08f.h5',
	     ('/home/jujuman/Research/ANI-DATASET/h5data/gdb9-2500-bad_new.h5','gdb9-2500-bad'),
	     ('/home/jujuman/Research/ANI-DATASET/h5data/gdb9-2500-div_new.h5','gdb9-2500-div'),
	     ('/home/jujuman/Research/ANI-DATASET/h5data/gdb9-2500-div-dim_35.h5','gdb9-2500-div-dim'),
	     ('/home/jujuman/Research/ANI-DATASET/h5data/sf_new.h5','sf_data'),
             ('/home/jujuman/Research/ANI-DATASET/h5data/ani-gdb-c08f.h5','gdbc08f'),]

#hdf5file = '/home/jujuman/Research/ANI-DATASET/ani-1_data_c03.h5'
storecac = '/home/jujuman/Research/SingleNetworkTest/cachec08f09augdbm6/'
saef   = "/home/jujuman/Research/SingleNetworkTest/sae_wb97x-631gd_HCNOFS.dat"
path = "/home/jujuman/Research/SingleNetworkTest/cachec08f09augdbm6/testset/testset.h5"


# Declare data cache
cachet = cg('_train', saef, storecac)
cachev = cg('_valid', saef, storecac)

# Declare test cache
dpack = pyt.datapacker(path)

for f in hdf5files:
    # Construct the data loader class
    print(f)
    adl = pyt.anidataloader(f[0])

    print(adl.get_group_list())

    # Loop over data in set
    dc = 0
    for i,data in enumerate(adl):

        xyz = np.array_split(data['coordinates'], 10)
        eng = np.array_split(data['energies'], 10)
        spc = data['species']
        nme = data['parent']

        #print('Parent: ', nme, eng)
        dc = dc + np.concatenate(eng[0:8]).shape[0]

        # Prepare and store the training and validation data
        cachet.insertdata(np.concatenate(xyz[0:8]), np.array(np.concatenate(eng[0:8]), dtype=np.float64), list(spc))
        cachev.insertdata(xyz[8], np.array(eng[8], dtype=np.float64), list(spc))

        Na = len(list(spc))

        # Prepare and store the test data set
        if xyz[9].shape[0] != 0:
            #print(xyz[9].shape)
            t_xyz = xyz[9].reshape(xyz[9].shape[0], Na*3)
            dpack.store_data(f[1] + '/mol' + str(i), coordinates=t_xyz, energies=np.array(eng[9]), species=spc)
    print('Count: ',dc)

    adl.cleanup()

# Make meta data file for caches
cachet.makemetadata()
cachev.makemetadata()

# Cleanup the disk
dpack.cleanup()
