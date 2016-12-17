from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

from sklearn.neural_network import MLPRegressor

from os import listdir
import time as tm

import numpy as np
import graphtools as gt
import coulomb_matrix as cm

cmdfile = '/home/jujuman/Research/CMatrixBaseline_data/test_cm_data.dat'

dtdir='/home/jujuman/Research/GDB-11-wB97X-6-31gd/traindata/'
#dtdir='/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_02/testdata/'
files = listdir(dtdir)

N = 30

print('Computing coulomb matrices...')
_t1b = tm.time()

cnt=0

for i in files:
    cnt += 1
    print('FILE: ' + str(cnt) + ' of ' + str(len(files)) + ' ' + i)

    # Set file
    file = dtdir + i

    # Get training molecules
    xyz_tr, typ_tr, Eact_tr, readf = gt.readncdat(file, np.float32)

    # Compute energy of atoms at infinite separation
    ise = cm.computerISE(typ_tr)

    Eact_tr = (Eact_tr - ise)

    Nm = Eact_tr.shape[0]
    cmat = cm.GenCMatData(xyz_tr,typ_tr,N)

    cmat[:, N * N ] = Eact_tr

    output_file = open(cmdfile + '_' + str(cnt), 'wb')
    cmat.tofile(output_file)
    output_file.close()


# End timer
_t1e = tm.time()
print('Computation complete in: ' + "{:.4f}".format((_t1e - _t1b)) + 's')
