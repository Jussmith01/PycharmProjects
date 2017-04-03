import numpy as np
import hdnntools as hdn

# Store test split
import pyanitools as pyt
from pyNeuroChem import cachegenerator as cg
import pyanitools as pya
import re

# Set the HDF5 file containing the data
hdf5file = '/home/jujuman/Research/ANI-DATASET/Roman-DATA/gdb9-2500-div.h5'
#hdf5file = '/home/jujuman/Research/ANI-DATASET/ani-1_data_c03.h5'
storecac = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/cache01/'
saef   = "/home/jujuman/Research/GDB-11-wB97X-6-31gd/sae_6-31gd.dat"
path = "/home/jujuman/Research/GDB-11-wB97X-6-31gd/cache01/testset/c01-testset.h5"

# Construct the data loader class
adl = pya.anidataloader(hdf5file)

# Declare data cache
cachet = cg('_train', saef, storecac)
cachev = cg('_valid', saef, storecac)

# Declare test cache
dpack = pyt.datapacker(path)

# Loop over data in set
for data in adl.getnextdata():

    xyz = np.array_split(data['coordinates'], 10)
    eng = np.array_split(data['energies'], 10)
    spc = data['species']

    # Prepare and store the training and validation data
    cachet.insertdata(np.concatenate(xyz[0:8])), np.array(np.concatenate(eng[0:8]), dtype=np.float64), list(spc))
    cachev.insertdata(xyz[8], np.array(eng[8], dtype=np.float64), list(spc))

    # Prepare and store the test data set
    if xyz[9].shape[0] != 0:
        t_xyz = xyz[9].reshape(xyz[9].shape[0],xyz[9].shape[1]*xyz[9].shape[2])
        dpack.store_data(loc, t_xyz, np.array(eng[9]), list(spc))

# Make meta data file for caches
cachet.makemetadata()
cachev.makemetadata()

# Cleanup the disk
dpack.cleanup()
adl.cleanup()