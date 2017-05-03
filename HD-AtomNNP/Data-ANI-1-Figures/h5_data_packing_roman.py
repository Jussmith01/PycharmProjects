import numpy as np
import pyanitools as pyt
import hdnntools as hdt

def check_for_outsider (okayl, chckl):
    for i in chckl:
        if i not in okayl:
            return False
    return True

dst = "/home/jujuman/Research/ANI-DATASET/h5data/gdb9-2500-bad_new.h5"
src = "/home/jujuman/Research/ANI-DATASET/GDB-09-Data/gdb9-2500-bad.h5"

#open an HDF5 for compressed storage.
#Note that if the path exists, it will open whatever is there.
dpack = pyt.datapacker(dst)
aload = pyt.anidataloader(src)

at = ['H',
      'C',
      'N',
      'O',
      #'F',
      #'S',
	]

for id,data in enumerate(aload.get_roman_data()):

    xyz = np.asarray(data['coordinates'], dtype=np.float32)
    erg = np.asarray(data['energies'], dtype=np.float64)
    spc = [str(a.decode('ascii')) for a in data['species']]

    if check_for_outsider(at,spc):
        print('Packing ', id, ' ...')
        dpack.store_data("/sfdata/mol" + str(id),coordinates=xyz, energies=erg, species=spc)
    else:
        print('Invalid type (', id, ')...')

dpack.cleanup()

